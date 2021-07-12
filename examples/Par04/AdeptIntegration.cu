// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#include "AdeptIntegration.h"
#include "AdeptIntegration.cuh"

#include <AdePT/Atomic.h>
#include <AdePT/LoopNavigator.h>
#include <AdePT/MParray.h>

#include <CopCore/Global.h>
#include <CopCore/PhysicalConstants.h>
#include <CopCore/Ranluxpp.h>

#include <VecGeom/base/Config.h>
#include <VecGeom/base/Stopwatch.h>
#ifdef VECGEOM_ENABLE_CUDA
#include <VecGeom/backend/cuda/Interface.h>
#endif

#include <G4TransportationManager.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4HepEmData.hh>
#include <G4HepEmElectronInit.hh>
#include <G4HepEmGammaInit.hh>
#include <G4HepEmMatCutData.hh>
#include <G4HepEmMaterialInit.hh>
#include <G4HepEmParameters.hh>
#include <G4HepEmParametersInit.hh>

#include <iostream>
#include <iomanip>
#include <stdio.h>

__constant__ __device__ struct G4HepEmParameters g4HepEmPars;
__constant__ __device__ struct G4HepEmData g4HepEmData;

__constant__ __device__ int *MCIndex        = nullptr;
__constant__ __device__ double BzFieldValue = 0;

__constant__ __device__ int Zero = 0;

// A bundle of queues per particle type:
//  * Two for active particles, one for the current iteration and the second for the next.
//  * One for all particles that need to be relocated to the next volume.
struct ParticleQueues {
  adept::MParray *currentlyActive;
  adept::MParray *nextActive;
  adept::MParray *relocate;

  void SwapActive() { std::swap(currentlyActive, nextActive); }
};

struct ParticleType {
  Track *tracks;
  SlotManager *slotManager;
  ParticleQueues queues;
  cudaStream_t stream;
  cudaEvent_t event;

  enum {
    Electron = 0,
    Positron = 1,
    Gamma    = 2,

    NumParticleTypes,
  };
};

// A bundle of queues for the three particle types.
struct AllParticleQueues {
  ParticleQueues queues[ParticleType::NumParticleTypes];
};

// A data structure to transfer statistics after each iteration.
struct Stats {
  int inFlight[ParticleType::NumParticleTypes];
};

struct G4HepEmState {
  G4HepEmData data;
  G4HepEmParameters parameters;
};

struct GPUstate {
  ParticleType particles[ParticleType::NumParticleTypes];
  G4HepEmState *hepem_state{nullptr};
  // Create a stream to synchronize kernels of all particle types.
  cudaStream_t stream;
  AdeptIntegration::TrackData *toDevice_dev{nullptr};
  Stats *stats_dev{nullptr};
  Stats *stats{nullptr};
};

static GPUstate gpuState[AdeptIntegration::kMaxThreads];

static G4HepEmState *InitG4HepEm()
{
  G4HepEmState *state = new G4HepEmState;
  InitG4HepEmData(&state->data);
  InitHepEmParameters(&state->parameters);

  InitMaterialAndCoupleData(&state->data, &state->parameters);

  InitElectronData(&state->data, &state->parameters, true);
  InitElectronData(&state->data, &state->parameters, false);
  InitGammaData(&state->data, &state->parameters);

  G4HepEmMatCutData *cutData = state->data.fTheMatCutData;
  std::cout << "fNumG4MatCuts = " << cutData->fNumG4MatCuts << ", fNumMatCutData = " << cutData->fNumMatCutData
            << std::endl;

  // Copy to GPU.
  CopyG4HepEmDataToGPU(&state->data);
  COPCORE_CUDA_CHECK(cudaMemcpyToSymbol(g4HepEmPars, &state->parameters, sizeof(G4HepEmParameters)));

  // Create G4HepEmData with the device pointers.
  G4HepEmData dataOnDevice;
  dataOnDevice.fTheMatCutData   = state->data.fTheMatCutData_gpu;
  dataOnDevice.fTheMaterialData = state->data.fTheMaterialData_gpu;
  dataOnDevice.fTheElementData  = state->data.fTheElementData_gpu;
  dataOnDevice.fTheElectronData = state->data.fTheElectronData_gpu;
  dataOnDevice.fThePositronData = state->data.fThePositronData_gpu;
  dataOnDevice.fTheSBTableData  = state->data.fTheSBTableData_gpu;
  dataOnDevice.fTheGammaData    = state->data.fTheGammaData_gpu;
  // The other pointers should never be used.
  dataOnDevice.fTheMatCutData_gpu   = nullptr;
  dataOnDevice.fTheMaterialData_gpu = nullptr;
  dataOnDevice.fTheElementData_gpu  = nullptr;
  dataOnDevice.fTheElectronData_gpu = nullptr;
  dataOnDevice.fThePositronData_gpu = nullptr;
  dataOnDevice.fTheSBTableData_gpu  = nullptr;
  dataOnDevice.fTheGammaData_gpu    = nullptr;

  COPCORE_CUDA_CHECK(cudaMemcpyToSymbol(g4HepEmData, &dataOnDevice, sizeof(G4HepEmData)));

  return state;
}

static void FreeG4HepEm(G4HepEmState *state)
{
  FreeG4HepEmData(&state->data);
  delete state;
}

// Kernel to initialize the set of queues per particle type.
__global__ void InitParticleQueues(ParticleQueues queues, size_t Capacity)
{
  adept::MParray::MakeInstanceAt(Capacity, queues.currentlyActive);
  adept::MParray::MakeInstanceAt(Capacity, queues.nextActive);
  adept::MParray::MakeInstanceAt(Capacity, queues.relocate);
}

// Kernel function to initialize tracks comming from a Geant4 buffer
__global__ void InitTracks(AdeptIntegration::TrackData *trackinfo, int ntracks, int event, Secondaries secondaries,
                           const vecgeom::VPlacedVolume *world, GlobalScoring *globalScoring)
{
  constexpr unsigned long kMaxTracks = 10000;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < ntracks; i += blockDim.x * gridDim.x) {
    ParticleGenerator *generator = nullptr;
    switch (trackinfo[i].pdg) {
    case 11:
      generator = &secondaries.electrons;
      atomicAdd(&globalScoring->numElectrons, 1);
      break;
    case -11:
      generator = &secondaries.positrons;
      atomicAdd(&globalScoring->numPositrons, 1);
      break;
    case 22:
      generator = &secondaries.gammas;
      atomicAdd(&globalScoring->numGammas, 1);
    };
    assert(generator != nullptr && "Unsupported pdg type");

    Track &track = generator->NextTrack();
    track.rngState.SetSeed(314159265 * (kMaxTracks * event + i + 1));
    track.energy       = trackinfo[i].energy;
    track.numIALeft[0] = -1.0;
    track.numIALeft[1] = -1.0;
    track.numIALeft[2] = -1.0;

    track.pos = {trackinfo[i].position[0], trackinfo[i].position[1], trackinfo[i].position[2]};
    track.dir = {trackinfo[i].direction[0], trackinfo[i].direction[1], trackinfo[i].direction[2]};
    track.currentState.Clear();
    LoopNavigator::LocatePointIn(world, track.pos, track.currentState, true);
    // nextState is initialized as needed.
  }
}

// Finish iteration: clear queues and fill statistics.
__global__ void FinishIteration(AllParticleQueues all, Stats *stats)
{
  for (int i = 0; i < ParticleType::NumParticleTypes; i++) {
    all.queues[i].currentlyActive->clear();
    stats->inFlight[i] = all.queues[i].nextActive->size();
    all.queues[i].relocate->clear();
  }
}

// Deposit energy of particles still in flight.
__global__ void DepositEnergy(Track *allTracks, const adept::MParray *queue, GlobalScoring *globalScoring,
                              ScoringPerVolume *scoringPerVolume)
{
  int queueSize = queue->size();
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < queueSize; i += blockDim.x * gridDim.x) {
    const int slot      = (*queue)[i];
    Track &currentTrack = allTracks[slot];
    auto volume         = currentTrack.currentState.Top();
    if (volume == nullptr) {
      // The particle left the world, why wasn't it killed before?!
      continue;
    }
    int volumeID = volume->id();

    double energy = currentTrack.energy;
    atomicAdd(&globalScoring->energyDeposit, energy);
    atomicAdd(&scoringPerVolume->energyDeposit[volumeID], energy);
  }
}

__global__ void ClearQueue(adept::MParray *queue)
{
  queue->clear();
}

void AdeptIntegration::InitializeGPU(const vecgeom::cxx::VPlacedVolume *world, int max_batch)
{
  // Bridge method from CPU to initialize geometry and physics on device

  // Upload geometry, asserting that the GeoManager must be closed. Geometry creation is a user recipe.
  auto &cudaManager = vecgeom::cxx::CudaManager::Instance();
  cudaManager.LoadGeometry(world);
  auto world_dev = cudaManager.Synchronize();
  assert(world_dev && "VecGeom geometry upload to GPU failed");

  // Initialize physics
  gpuState[0].hepem_state = InitG4HepEm();

  // Initialize field
  auto field =
      (G4UniformMagField *)G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
  auto field_vect = field->GetConstantFieldValue();
  double bz       = field_vect[2];
  COPCORE_CUDA_CHECK(cudaMemcpyToSymbol(BzFieldValue, &bz, sizeof(double)));

  // Transfer MC indices.
  COPCORE_CUDA_CHECK(cudaMalloc(&fUserData[0].MCIndex_dev, sizeof(int) * NumVolumes));
  COPCORE_CUDA_CHECK(
      cudaMemcpy(fUserData[0].MCIndex_dev, fUserData[0].MCIndex, sizeof(int) * NumVolumes, cudaMemcpyHostToDevice));
  COPCORE_CUDA_CHECK(cudaMemcpyToSymbol(MCIndex, &fUserData[0].MCIndex_dev, sizeof(int *)));

  // Capacity of the different containers aka the maximum number of particles.
  constexpr int Capacity = 256 * 1024;
  std::cout << "INFO: batching " << max_batch << " particles for transport on the GPU" << std::endl;

  // Allocate structures to manage tracks of an implicit type:
  //  * memory to hold the actual Track elements,
  //  * objects to manage slots inside the memory,
  //  * queues of slots to remember active particle and those needing relocation,
  //  * a stream and an event for synchronization of kernels.
  constexpr size_t TracksSize  = sizeof(Track) * Capacity;
  constexpr size_t ManagerSize = sizeof(SlotManager);
  const size_t QueueSize       = adept::MParray::SizeOfInstance(Capacity);
  for (int tid = 0; tid < fNthreads; ++tid) {
    // Create a stream to synchronize kernels of all particle types.
    COPCORE_CUDA_CHECK(cudaStreamCreate(&gpuState[tid].stream));
    if (tid > 0) {
      fUserData[tid].MCIndex_dev = fUserData[0].MCIndex_dev;
      gpuState[tid].hepem_state  = gpuState[0].hepem_state;
    }
    for (int i = 0; i < ParticleType::NumParticleTypes; i++) {
      // Share hepem state between threads
      COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].particles[i].tracks, TracksSize));

      COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].particles[i].slotManager, ManagerSize));

      COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].particles[i].queues.currentlyActive, QueueSize));
      COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].particles[i].queues.nextActive, QueueSize));
      COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].particles[i].queues.relocate, QueueSize));
      InitParticleQueues<<<1, 1>>>(gpuState[tid].particles[i].queues, Capacity);

      COPCORE_CUDA_CHECK(cudaStreamCreate(&gpuState[tid].particles[i].stream));
      COPCORE_CUDA_CHECK(cudaEventCreate(&gpuState[tid].particles[i].event));
    }
  }
  COPCORE_CUDA_CHECK(cudaDeviceSynchronize());

  for (int tid = 0; tid < fNthreads; ++tid) {
    // Allocate memory to score charged track length and energy deposit per volume.
    COPCORE_CUDA_CHECK(cudaMalloc(&fUserData[tid].chargedTrackLength_dev, sizeof(double) * NumVolumes));
    COPCORE_CUDA_CHECK(cudaMemset(fUserData[tid].chargedTrackLength_dev, 0, sizeof(double) * NumVolumes));
    COPCORE_CUDA_CHECK(cudaMalloc(&fUserData[tid].energyDeposit_dev, sizeof(double) * NumVolumes));
    COPCORE_CUDA_CHECK(cudaMemset(fUserData[tid].energyDeposit_dev, 0, sizeof(double) * NumVolumes));

    // Allocate and initialize scoring and statistics.
    COPCORE_CUDA_CHECK(cudaMalloc(&fUserData[tid].globalScoring_dev, sizeof(GlobalScoring)));
    COPCORE_CUDA_CHECK(cudaMemset(fUserData[tid].globalScoring_dev, 0, sizeof(GlobalScoring)));

    ScoringPerVolume scoringPerVolume_devPtrs;
    scoringPerVolume_devPtrs.chargedTrackLength = fUserData[tid].chargedTrackLength_dev;
    scoringPerVolume_devPtrs.energyDeposit      = fUserData[tid].energyDeposit_dev;
    COPCORE_CUDA_CHECK(cudaMalloc(&fUserData[tid].scoringPerVolume_dev, sizeof(ScoringPerVolume)));
    COPCORE_CUDA_CHECK(cudaMemcpy(fUserData[tid].scoringPerVolume_dev, &scoringPerVolume_devPtrs,
                                  sizeof(ScoringPerVolume), cudaMemcpyHostToDevice));

    // initialize statistics
    COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].stats_dev, sizeof(Stats)));
    COPCORE_CUDA_CHECK(cudaMallocHost(&gpuState[tid].stats, sizeof(Stats)));

    // initialize buffer of tracks on device
    COPCORE_CUDA_CHECK(cudaMalloc(&gpuState[tid].toDevice_dev, max_batch * sizeof(TrackData)));
  }
}

void AdeptIntegration::FreeGPU()
{
  // Free resources.
  for (int tid = 0; tid < fNthreads; ++tid) {
    COPCORE_CUDA_CHECK(cudaFree(fUserData[tid].chargedTrackLength_dev));
    COPCORE_CUDA_CHECK(cudaFree(fUserData[tid].energyDeposit_dev));

    COPCORE_CUDA_CHECK(cudaFree(fUserData[tid].globalScoring_dev));
    COPCORE_CUDA_CHECK(cudaFree(fUserData[tid].scoringPerVolume_dev));
    COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].stats_dev));
    COPCORE_CUDA_CHECK(cudaFreeHost(gpuState[tid].stats));
    COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].toDevice_dev));

    COPCORE_CUDA_CHECK(cudaStreamDestroy(gpuState[tid].stream));

    for (int i = 0; i < ParticleType::NumParticleTypes; i++) {
      COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].particles[i].tracks));
      COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].particles[i].slotManager));

      COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].particles[i].queues.currentlyActive));
      COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].particles[i].queues.nextActive));
      COPCORE_CUDA_CHECK(cudaFree(gpuState[tid].particles[i].queues.relocate));

      COPCORE_CUDA_CHECK(cudaStreamDestroy(gpuState[tid].particles[i].stream));
      COPCORE_CUDA_CHECK(cudaEventDestroy(gpuState[tid].particles[i].event));
    }
  }

  COPCORE_CUDA_CHECK(cudaFree(fUserData[0].MCIndex_dev));
  FreeG4HepEm(gpuState[0].hepem_state);
}

void AdeptIntegration::ShowerGPU(int event, int tid, TrackBuffer const &buffer)
{
  // Capacity of the different containers aka the maximum number of particles.
  constexpr int Capacity                        = 256 * 1024;
  auto &cudaManager                             = vecgeom::cxx::CudaManager::Instance();
  const vecgeom::cuda::VPlacedVolume *world_dev = cudaManager.world_gpu();

  ParticleType &electrons = gpuState[tid].particles[ParticleType::Electron];
  ParticleType &positrons = gpuState[tid].particles[ParticleType::Positron];
  ParticleType &gammas    = gpuState[tid].particles[ParticleType::Gamma];

  vecgeom::Stopwatch timer;
  timer.Start();

  // copy buffer of tracks to device
  COPCORE_CUDA_CHECK(cudaMemcpy(gpuState[tid].toDevice_dev, buffer.toDevice.data(),
                                buffer.toDevice.size() * sizeof(TrackData), cudaMemcpyHostToDevice));

  // initialize slot manager
  SlotManager slotManagerInit(Capacity);
  for (int i = 0; i < ParticleType::NumParticleTypes; i++) {
    COPCORE_CUDA_CHECK(cudaMemcpy(gpuState[tid].particles[i].slotManager, &slotManagerInit, sizeof(SlotManager),
                                  cudaMemcpyHostToDevice));
  }

  std::cout << std::endl << "GPU transporting event " << event << " for CPU thread " << tid << " ..." << std::flush;

  // Initialize AdePT tracks using the track buffer copied from CPU
  constexpr int initThreads = 32;
  int initBlocks            = (buffer.toDevice.size() + initThreads - 1) / initThreads;
  Secondaries secondaries   = {
      .electrons = {electrons.tracks, electrons.slotManager, electrons.queues.nextActive},
      .positrons = {positrons.tracks, positrons.slotManager, positrons.queues.nextActive},
      .gammas    = {gammas.tracks, gammas.slotManager, gammas.queues.nextActive},
  };

  InitTracks<<<initBlocks, initThreads>>>(gpuState[tid].toDevice_dev, buffer.toDevice.size(), event, secondaries,
                                          world_dev, fUserData[tid].globalScoring_dev);

  COPCORE_CUDA_CHECK(cudaDeviceSynchronize());

  gpuState[tid].stats->inFlight[ParticleType::Electron] = fBuffer[tid].nelectrons;
  gpuState[tid].stats->inFlight[ParticleType::Positron] = fBuffer[tid].npositrons;
  gpuState[tid].stats->inFlight[ParticleType::Gamma]    = fBuffer[tid].ngammas;

  constexpr int MaxBlocks        = 1024;
  constexpr int TransportThreads = 32;
  constexpr int RelocateThreads  = 32;
  int transportBlocks, relocateBlocks;

  int inFlight;
  int iterNo = 0;

  do {

    // *** ELECTRONS ***
    int numElectrons = gpuState[tid].stats->inFlight[ParticleType::Electron];
    if (numElectrons > 0) {
      transportBlocks = (numElectrons + TransportThreads - 1) / TransportThreads;
      transportBlocks = std::min(transportBlocks, MaxBlocks);

      relocateBlocks = std::min(numElectrons, MaxBlocks);

      TransportElectrons<<<transportBlocks, TransportThreads, 0, electrons.stream>>>(
          electrons.tracks, electrons.queues.currentlyActive, secondaries, electrons.queues.nextActive,
          electrons.queues.relocate, fUserData[tid].globalScoring_dev, fUserData[tid].scoringPerVolume_dev);

      RelocateToNextVolume<<<relocateBlocks, RelocateThreads, 0, electrons.stream>>>(electrons.tracks,
                                                                                     electrons.queues.relocate);

      COPCORE_CUDA_CHECK(cudaEventRecord(electrons.event, electrons.stream));
      COPCORE_CUDA_CHECK(cudaStreamWaitEvent(gpuState[tid].stream, electrons.event, 0));
    }

    // *** POSITRONS ***
    int numPositrons = gpuState[tid].stats->inFlight[ParticleType::Positron];
    if (numPositrons > 0) {
      transportBlocks = (numPositrons + TransportThreads - 1) / TransportThreads;
      transportBlocks = std::min(transportBlocks, MaxBlocks);

      relocateBlocks = std::min(numPositrons, MaxBlocks);

      TransportPositrons<<<transportBlocks, TransportThreads, 0, positrons.stream>>>(
          positrons.tracks, positrons.queues.currentlyActive, secondaries, positrons.queues.nextActive,
          positrons.queues.relocate, fUserData[tid].globalScoring_dev, fUserData[tid].scoringPerVolume_dev);

      RelocateToNextVolume<<<relocateBlocks, RelocateThreads, 0, positrons.stream>>>(positrons.tracks,
                                                                                     positrons.queues.relocate);

      COPCORE_CUDA_CHECK(cudaEventRecord(positrons.event, positrons.stream));
      COPCORE_CUDA_CHECK(cudaStreamWaitEvent(gpuState[tid].stream, positrons.event, 0));
    }

    // *** GAMMAS ***
    int numGammas = gpuState[tid].stats->inFlight[ParticleType::Gamma];
    if (numGammas > 0) {
      transportBlocks = (numGammas + TransportThreads - 1) / TransportThreads;
      transportBlocks = std::min(transportBlocks, MaxBlocks);

      relocateBlocks = std::min(numGammas, MaxBlocks);

      TransportGammas<<<transportBlocks, TransportThreads, 0, gammas.stream>>>(
          gammas.tracks, gammas.queues.currentlyActive, secondaries, gammas.queues.nextActive, gammas.queues.relocate,
          fUserData[tid].globalScoring_dev, fUserData[tid].scoringPerVolume_dev);

      RelocateToNextVolume<<<relocateBlocks, RelocateThreads, 0, gammas.stream>>>(gammas.tracks,
                                                                                  gammas.queues.relocate);

      COPCORE_CUDA_CHECK(cudaEventRecord(gammas.event, gammas.stream));
      COPCORE_CUDA_CHECK(cudaStreamWaitEvent(gpuState[tid].stream, gammas.event, 0));
    }

    // *** END OF TRANSPORT ***

    // The events ensure synchronization before finishing this iteration and
    // copying the Stats back to the host.
    AllParticleQueues queues = {{electrons.queues, positrons.queues, gammas.queues}};
    FinishIteration<<<1, 1, 0, gpuState[tid].stream>>>(queues, gpuState[tid].stats_dev);
    COPCORE_CUDA_CHECK(cudaMemcpyAsync(gpuState[tid].stats, gpuState[tid].stats_dev, sizeof(Stats),
                                       cudaMemcpyDeviceToHost, gpuState[tid].stream));

    // Finally synchronize all kernels.
    COPCORE_CUDA_CHECK(cudaStreamSynchronize(gpuState[tid].stream));

    // Count the number of particles in flight.
    inFlight = 0;
    for (int i = 0; i < ParticleType::NumParticleTypes; i++) {
      inFlight += gpuState[tid].stats->inFlight[i];
    }

    // Swap the queues for the next iteration.
    electrons.queues.SwapActive();
    positrons.queues.SwapActive();
    gammas.queues.SwapActive();

    // Update the active queues that for next iteration
    secondaries.electrons.SetActiveQueue(electrons.queues.nextActive);
    secondaries.positrons.SetActiveQueue(positrons.queues.nextActive);
    secondaries.gammas.SetActiveQueue(gammas.queues.nextActive);

    iterNo++;
  } while (inFlight > 0 && iterNo < 1000);

  if (inFlight > 0) {
    std::cout << std::endl;
    std::cout << "WARN: Depositing energy of " << inFlight << " particles still in flight!" << std::endl;
    constexpr int DepositThreads = 32;

    for (int i = 0; i < ParticleType::NumParticleTypes; i++) {
      ParticleType &pType   = gpuState[tid].particles[i];
      int inFlightParticles = gpuState[tid].stats->inFlight[i];
      if (inFlightParticles == 0) {
        continue;
      }

      int depositBlocks = (inFlightParticles + DepositThreads - 1) / DepositThreads;
      depositBlocks     = std::min(depositBlocks, MaxBlocks);
      DepositEnergy<<<depositBlocks, DepositThreads, 0, gpuState[tid].stream>>>(
          pType.tracks, pType.queues.currentlyActive, fUserData[tid].globalScoring_dev,
          fUserData[tid].scoringPerVolume_dev);

      ClearQueue<<<1, 1, 0, gpuState[tid].stream>>>(pType.queues.currentlyActive);
    }
    COPCORE_CUDA_CHECK(cudaStreamSynchronize(gpuState[tid].stream));
    std::cout << " ... ";
  }

  std::cout << "done!" << std::endl;

  auto time = timer.Stop();
  std::cout << "Run time: " << time << "\n";

  // Transfer back scoring.
  COPCORE_CUDA_CHECK(cudaMemcpy(&fUserData[tid].globalScoring, fUserData[tid].globalScoring_dev, sizeof(GlobalScoring),
                                cudaMemcpyDeviceToHost));

  // Transfer back the scoring per volume (charged track length and energy deposit).
  COPCORE_CUDA_CHECK(cudaMemcpy(fUserData[tid].scoringPerVolume.chargedTrackLength,
                                fUserData[tid].chargedTrackLength_dev, sizeof(double) * NumVolumes,
                                cudaMemcpyDeviceToHost));
  COPCORE_CUDA_CHECK(cudaMemcpy(fUserData[tid].scoringPerVolume.energyDeposit, fUserData[tid].energyDeposit_dev,
                                sizeof(double) * NumVolumes, cudaMemcpyDeviceToHost));
  fUserData[tid].globalScoring.Print();
}
