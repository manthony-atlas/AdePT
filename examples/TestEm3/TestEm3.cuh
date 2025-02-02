// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#ifndef TESTEM3_CUH
#define TESTEM3_CUH

#include "TestEm3.h"

#include <AdePT/MParray.h>
#include <CopCore/SystemOfUnits.h>
#include <CopCore/Ranluxpp.h>

#include <G4HepEmData.hh>
#include <G4HepEmParameters.hh>
#include <G4HepEmRandomEngine.hh>

#include <VecGeom/base/Vector3D.h>
#include <VecGeom/navigation/NavStateIndex.h>


// A data structure to represent a particle track. The particle type is implicit
// by the queue and not stored in memory.
struct Track {
  RanluxppDouble rngState;
  double energy;
  double numIALeft[3];

  vecgeom::Vector3D<double> pos;
  vecgeom::Vector3D<double> dir;
  vecgeom::NavStateIndex currentState;
  vecgeom::NavStateIndex nextState;

  __host__ __device__ double Uniform() { return rngState.Rndm(); }

  __host__ __device__ void SwapStates()
  {
    auto state         = this->currentState;
    this->currentState = this->nextState;
    this->nextState    = state;
  }

  __host__ __device__ void InitAsSecondary(const Track &parent)
  {
    // Initialize a new PRNG state.
    this->rngState = parent.rngState;
    this->rngState.Skip(1 << 15);

    // The caller is responsible to set the energy.
    this->numIALeft[0] = -1.0;
    this->numIALeft[1] = -1.0;
    this->numIALeft[2] = -1.0;

    // A secondary inherits the position of its parent; the caller is responsible
    // to update the directions.
    this->pos           = parent.pos;
    this->currentState = parent.currentState;
    this->nextState    = parent.nextState;
  }
};

// Defined in TestEm3.cu
extern __constant__ __device__ int Zero;

class RanluxppDoubleEngine : public G4HepEmRandomEngine {
  // Wrapper functions to call into RanluxppDouble.
  static __host__ __device__ __attribute__((noinline))
  double FlatWrapper(void *object)
  {
    return ((RanluxppDouble *)object)->Rndm();
  }
  static __host__ __device__ __attribute__((noinline))
  void FlatArrayWrapper(void *object, const int size, double *vect)
  {
    for (int i = 0; i < size; i++) {
      vect[i] = ((RanluxppDouble *)object)->Rndm();
    }
  }

public:
  __host__ __device__ RanluxppDoubleEngine(RanluxppDouble *engine)
      : G4HepEmRandomEngine(/*object=*/engine, &FlatWrapper, &FlatArrayWrapper)
  {
#ifdef __CUDA_ARCH__
    // This is a hack: The compiler cannot see that we're going to call the
    // functions through their pointers, so it underestimates the number of
    // required registers. By including calls to the (non-inlinable) functions
    // we force the compiler to account for the register usage, even if this
    // particular set of calls are not executed at runtime.
    if (Zero) {
      FlatWrapper(engine);
      FlatArrayWrapper(engine, 0, nullptr);
    }
#endif
  }
};


// A data structure to manage slots in the track storage.
class SlotManager {
  adept::Atomic_t<int> fNextSlot;
  const int fMaxSlot;

public:
  __host__ __device__ SlotManager(int maxSlot) : fMaxSlot(maxSlot) { fNextSlot = 0; }

  __host__ __device__ int NextSlot()
  {
    int next = fNextSlot.fetch_add(1);
    if (next >= fMaxSlot) return -1;
    return next;
  }
};

// A bundle of pointers to generate particles of an implicit type.
class ParticleGenerator {
  Track *fTracks;
  SlotManager *fSlotManager;
  adept::MParray *fActiveQueue;

public:
  __host__ __device__ ParticleGenerator(Track *tracks, SlotManager *slotManager, adept::MParray *activeQueue)
    : fTracks(tracks), fSlotManager(slotManager), fActiveQueue(activeQueue) {}

  __host__ __device__ Track &NextTrack()
  {
    int slot = fSlotManager->NextSlot();
    if (slot == -1) {
      COPCORE_EXCEPTION("No slot available in ParticleGenerator::NextTrack");
    }
    fActiveQueue->push_back(slot);
    return fTracks[slot];
  }
};

// A bundle of generators for the three particle types.
struct Secondaries {
  ParticleGenerator electrons;
  ParticleGenerator positrons;
  ParticleGenerator gammas;
};


// Kernels in different TUs.
__global__ void RelocateToNextVolume(Track *allTracks, const adept::MParray *relocateQueue);

__global__ void TransportElectrons(
    Track *electrons, const adept::MParray *active, Secondaries secondaries, adept::MParray *activeQueue,
    adept::MParray *relocateQueue, GlobalScoring *globalScoring, ScoringPerVolume *scoringPerVolume);
__global__ void TransportPositrons(
    Track *positrons, const adept::MParray *active, Secondaries secondaries, adept::MParray *activeQueue,
    adept::MParray *relocateQueue, GlobalScoring *globalScoring, ScoringPerVolume *scoringPerVolume);

__global__ void TransportGammas(Track *gammas, const adept::MParray *active, Secondaries secondaries,
                                adept::MParray *activeQueue, adept::MParray *relocateQueue,
                                GlobalScoring *globalScoring, ScoringPerVolume *scoringPerVolume);

// Constant data structures from G4HepEm accessed by the kernels.
// (defined in TestEm3.cu)
extern __constant__ __device__ struct G4HepEmParameters g4HepEmPars;
extern __constant__ __device__ struct G4HepEmData g4HepEmData;

extern __constant__ __device__ int *MCIndex;

// constexpr float BzFieldValue = 0.1 * copcore::units::tesla;
constexpr double BzFieldValue = 0;
constexpr double kPush = 1.e-8 * copcore::units::cm;

#endif
