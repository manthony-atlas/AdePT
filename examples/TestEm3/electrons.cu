// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#include "TestEm3.cuh"

#include <AdePT/LoopNavigator.h>
#include <fieldPropagatorConstBz.h>

#include <CopCore/PhysicalConstants.h>

#include <G4HepEmElectronManager.hh>
#include <G4HepEmElectronTrack.hh>
#include <G4HepEmElectronInteractionBrem.hh>
#include <G4HepEmElectronInteractionIoni.hh>
#include <G4HepEmPositronInteractionAnnihilation.hh>
// Pull in implementation.
#include <G4HepEmRunUtils.icc>
#include <G4HepEmInteractionUtils.icc>
#include <G4HepEmElectronManager.icc>
#include <G4HepEmElectronInteractionBrem.icc>
#include <G4HepEmElectronInteractionIoni.icc>
#include <G4HepEmPositronInteractionAnnihilation.icc>

__device__ struct G4HepEmElectronManager electronManager;

// Compute the physics and geometry step limit, transport the electrons while
// applying the continuous effects and maybe a discrete process that could
// generate secondaries.
template <bool IsElectron>
static __device__ __forceinline__ void TransportElectrons(Track *electrons, 
							  const adept::MParray *active,
                                                          Secondaries &secondaries, adept::MParray *activeQueue,
                                                          adept::MParray *relocateQueue, GlobalScoring *globalScoring,
                                                          ScoringPerVolume *scoringPerVolume,
							  HitRecord *hitRecord,
							  int *eventNumber,
							  ScoringPerParticle *scoringPerParticle,
							  int *mystream_stride
							  )
{
  constexpr int Charge  = IsElectron ? -1 : 1;
  constexpr double Mass = copcore::units::kElectronMassC2;
  fieldPropagatorConstBz fieldPropagatorBz(BzFieldValue);

  int activeSize = active->size();
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < activeSize; i += blockDim.x * gridDim.x) {
    const int slot      = (*active)[i];
    Track &currentTrack = electrons[slot];
    auto volume         = currentTrack.currentState.Top();
    if (volume == nullptr) {
      // The particle left the world, kill it by not enqueuing into activeQueue.
      continue;
    }
    int volumeID   = volume->id();
    int theMCIndex = MCIndex[volumeID];

    // Init a track with the needed data to call into G4HepEm.
    G4HepEmElectronTrack elTrack;
    G4HepEmTrack *theTrack = elTrack.GetTrack();
    theTrack->SetEKin(currentTrack.energy);
    theTrack->SetMCIndex(theMCIndex);
    theTrack->SetCharge(Charge);

    // Sample the `number-of-interaction-left` and put it into the track.
    for (int ip = 0; ip < 3; ++ip) {
      double numIALeft = currentTrack.numIALeft[ip];
      if (numIALeft <= 0) {
        numIALeft                  = -std::log(currentTrack.Uniform());
        currentTrack.numIALeft[ip] = numIALeft;
      }
      theTrack->SetNumIALeft(numIALeft, ip);
    }

    // Call G4HepEm to compute the physics step limit.
    electronManager.HowFar(&g4HepEmData, &g4HepEmPars, &elTrack);

    // Get result into variables.
    double geometricalStepLengthFromPhysics = theTrack->GetGStepLength();
    // The phyiscal step length is the amount that the particle experiences
    // which might be longer than the geometrical step length due to MSC. As
    // long as we call PerformContinuous in the same kernel we don't need to
    // care, but we need to make this available when splitting the operations.
    // double physicalStepLength = elTrack.GetPStepLength();
    int winnerProcessIndex = theTrack->GetWinnerProcessIndex();
    // Leave the range and MFP inside the G4HepEmTrack. If we split kernels, we
    // also need to carry them over!

    // Check if there's a volume boundary in between.
    double geometryStepLength;
    if (BzFieldValue != 0) {
      geometryStepLength = fieldPropagatorBz.ComputeStepAndPropagatedState</*Relocate=*/false>(
          currentTrack.energy, Mass, Charge, geometricalStepLengthFromPhysics, currentTrack.pos, currentTrack.dir,
          currentTrack.currentState, currentTrack.nextState);
    } else {
      geometryStepLength =
          LoopNavigator::ComputeStepAndNextVolume(currentTrack.pos, currentTrack.dir, geometricalStepLengthFromPhysics,
                                                  currentTrack.currentState, currentTrack.nextState);
      currentTrack.pos += (geometryStepLength + kPush) * currentTrack.dir;
    }
    atomicAdd(&globalScoring->chargedSteps, 1);
    atomicAdd(&scoringPerVolume->chargedTrackLength[volumeID], geometryStepLength);

    if (currentTrack.nextState.IsOnBoundary()) {
      theTrack->SetGStepLength(geometryStepLength);
      theTrack->SetOnBoundary(true);
    }

    // Apply continuous effects.
    bool stopped = electronManager.PerformContinuous(&g4HepEmData, &g4HepEmPars, &elTrack);
    // Collect the changes.
    currentTrack.energy  = theTrack->GetEKin();
    double energyDeposit = theTrack->GetEnergyDeposit();
    atomicAdd(&globalScoring->energyDeposit, energyDeposit);
    atomicAdd(&scoringPerVolume->energyDeposit[volumeID], energyDeposit);
    //now that we know it hit something, read off where it hit
    
    
    // Save the `number-of-interaction-left` in our track.
    for (int ip = 0; ip < 3; ++ip) {
      double numIALeft           = theTrack->GetNumIALeft(ip);
      currentTrack.numIALeft[ip] = numIALeft;
    }

    if (stopped) {
      if (!IsElectron) {
        // Annihilate the stopped positron into two gammas heading to opposite
        // directions (isotropic).
        Track &gamma1 = secondaries.gammas.NextTrack();
        Track &gamma2 = secondaries.gammas.NextTrack();
        atomicAdd(&globalScoring->numGammas, 2);

        const double cost = 2 * currentTrack.Uniform() - 1;
        const double sint = sqrt(1 - cost * cost);
        const double phi  = k2Pi * currentTrack.Uniform();
        double sinPhi, cosPhi;
        sincos(phi, &sinPhi, &cosPhi);

        gamma1.InitAsSecondary(/*parent=*/currentTrack);
        gamma1.energy = copcore::units::kElectronMassC2;
        gamma1.dir.Set(sint * cosPhi, sint * sinPhi, cost);

        gamma2.InitAsSecondary(/*parent=*/currentTrack);
        gamma2.energy = copcore::units::kElectronMassC2;
        gamma2.dir    = -gamma1.dir;
      }
      // Particles are killed by not enqueuing them into the new activeQueue.
      continue;
    }

    if (currentTrack.nextState.IsOnBoundary()) {
      // For now, just count that we hit something.
      atomicAdd(&globalScoring->hits, 1);
      // additional per volume info
      atomicAdd(&scoringPerVolume->numHits[volumeID],1);
      //now fill our event level counter
      int evtnum=*eventNumber;      
      atomicAdd(&scoringPerParticle->numHits_per_particle[evtnum],1);            
 
      //now that we know it hit something get the position of said hit
      vecgeom::Vector3D<vecgeom::Precision> hit_location=currentTrack.pos;
      //retrieve the x,y,z coordinates from the hit location
      double hit_location_x_pos=hit_location[0];
      double hit_location_y_pos=hit_location[1];
      double hit_location_z_pos=hit_location[2];
      // let's see if we can also get the thread ID
      int streamID=0;
      int streamstride=*mystream_stride;
      if(IsElectron){
	streamID=0;
      }
      else{
	streamID=1;
      }
      int streamIndex=streamID*streamstride;
      
      int threadIndex=threadIdx.x+blockIdx.x*blockDim.x; 
      int threadID=threadIdx.x;
      int blockID=blockIdx.x;
      // de-ref ptr to get original val

      int threadStreamIndex=threadIdx.x+blockDim.x*(blockIdx.x+streamIndex);

      //now get the hit
      int hit_number=scoringPerParticle->numHits_per_particle[evtnum];
      int hit_stride=2e6; //

      int nStream=3; // currently there are 3 streams only in use for particle transport - find a better way to do it at some point
      // 4 indices mapped into a single one: particle/batch, stream, block,thread
      // basically this maps what might be otherwise 4-dimensional (hit,stream,block,thread) into a 1D array for CUDA
      // uses formula (X,Y,Z,W) -> dim4_index= x+ y*dim(X) + z*dim(Y)*dim(X)+w*dim(Z)*dim(Y)*dim(X).
      // effectively this procedure is a flattening "row-wise " in 4D, and involves modular maths of strides. This flattening also produces a unique index guaranteed by the stride definition.
      long int index_first_term=hit_number; // x
      long int index_second_term=streamID*hit_stride; // y*dim(X)
      long int index_third_term=blockIdx.x * nStream * hit_stride; // z*dim(Y)*dim(X)
      long int index_fourth_term=threadIdx.x* nStream* hit_stride* blockDim.x; // w*dim(Z)*dim(Y)*dim(X)
      long int threadStreamHitindex=index_first_term+index_second_term+index_third_term+index_fourth_term;
      
      
      
      printf("Event Number: %i, Thread ID %i, Block ID %i, Thread Index %i, Stream ID %i, Stream+ThreadID %i, Streamstride %i, Unique hit identifier %li \n",
	     evtnum,
	     threadID,
	     blockID,
	     threadIndex,
	     streamID,
	     threadStreamIndex,
	     streamstride,
	     threadStreamHitindex
	     );

      //debug to print the unique identifier for each thread and stream combination
      //      printf("Event No.: %i ThreadID %i, Stream %li",evtnum,threadIndex,streamID);
      //      printf("Event No.: %i, ThreadID %i, Stream %li, VolumeID %i, Hit location: X= %f Y= %f Z=%f \n",evtnum,streamID,threadIndex,volumeID,hit_location_x_pos,hit_location_y_pos,hit_location_z_pos);


      activeQueue->push_back(slot);
      relocateQueue->push_back(slot);

      // Move to the next boundary.
      currentTrack.SwapStates();
      continue;
    } else if (winnerProcessIndex < 0) {
      // No discrete process, move on.
      activeQueue->push_back(slot);
      continue;
    }

    // Reset number of interaction left for the winner discrete process.
    // (Will be resampled in the next iteration.)
    currentTrack.numIALeft[winnerProcessIndex] = -1.0;

    // Check if a delta interaction happens instead of the real discrete process.
    if (electronManager.CheckDelta(&g4HepEmData, theTrack, currentTrack.Uniform())) {
      // A delta interaction happened, move on.
      activeQueue->push_back(slot);
      continue;
    }

    // Perform the discrete interaction.
    RanluxppDoubleEngine rnge(&currentTrack.rngState);

    const double energy   = currentTrack.energy;
    const double theElCut = g4HepEmData.fTheMatCutData->fMatCutData[theMCIndex].fSecElProdCutE;

    switch (winnerProcessIndex) {
    case 0: {
      // Invoke ionization (for e-/e+):
      double deltaEkin = (IsElectron) ? SampleETransferMoller(theElCut, energy, &rnge)
                                      : SampleETransferBhabha(theElCut, energy, &rnge);

      double dirPrimary[] = {currentTrack.dir.x(), currentTrack.dir.y(), currentTrack.dir.z()};
      double dirSecondary[3];
      SampleDirectionsIoni(energy, deltaEkin, dirSecondary, dirPrimary, &rnge);

      Track &secondary = secondaries.electrons.NextTrack();
      atomicAdd(&globalScoring->numElectrons, 1);

      secondary.InitAsSecondary(/*parent=*/currentTrack);
      secondary.energy = deltaEkin;
      secondary.dir.Set(dirSecondary[0], dirSecondary[1], dirSecondary[2]);

      currentTrack.energy = energy - deltaEkin;
      currentTrack.dir.Set(dirPrimary[0], dirPrimary[1], dirPrimary[2]);
      // The current track continues to live.
      activeQueue->push_back(slot);
      break;
    }
    case 1: {
      // Invoke model for Bremsstrahlung: either SB- or Rel-Brem.
      double logEnergy = std::log(energy);
      double deltaEkin = energy < g4HepEmPars.fElectronBremModelLim
                             ? SampleETransferBremSB(&g4HepEmData, energy, logEnergy, theMCIndex, &rnge, IsElectron)
                             : SampleETransferBremRB(&g4HepEmData, energy, logEnergy, theMCIndex, &rnge, IsElectron);

      double dirPrimary[] = {currentTrack.dir.x(), currentTrack.dir.y(), currentTrack.dir.z()};
      double dirSecondary[3];
      SampleDirectionsBrem(energy, deltaEkin, dirSecondary, dirPrimary, &rnge);

      Track &gamma = secondaries.gammas.NextTrack();
      atomicAdd(&globalScoring->numGammas, 1);

      gamma.InitAsSecondary(/*parent=*/currentTrack);
      gamma.energy = deltaEkin;
      gamma.dir.Set(dirSecondary[0], dirSecondary[1], dirSecondary[2]);

      currentTrack.energy = energy - deltaEkin;
      currentTrack.dir.Set(dirPrimary[0], dirPrimary[1], dirPrimary[2]);
      // The current track continues to live.
      activeQueue->push_back(slot);
      break;
    }
    case 2: {
      // Invoke annihilation (in-flight) for e+
      double dirPrimary[] = {currentTrack.dir.x(), currentTrack.dir.y(), currentTrack.dir.z()};
      double theGamma1Ekin, theGamma2Ekin;
      double theGamma1Dir[3], theGamma2Dir[3];
      SampleEnergyAndDirectionsForAnnihilationInFlight(energy, dirPrimary, &theGamma1Ekin, theGamma1Dir, &theGamma2Ekin,
                                                       theGamma2Dir, &rnge);

      Track &gamma1 = secondaries.gammas.NextTrack();
      Track &gamma2 = secondaries.gammas.NextTrack();
      atomicAdd(&globalScoring->numGammas, 2);

      gamma1.InitAsSecondary(/*parent=*/currentTrack);
      gamma1.energy = theGamma1Ekin;
      gamma1.dir.Set(theGamma1Dir[0], theGamma1Dir[1], theGamma1Dir[2]);

      gamma2.InitAsSecondary(/*parent=*/currentTrack);
      gamma2.energy = theGamma2Ekin;
      gamma2.dir.Set(theGamma2Dir[0], theGamma2Dir[1], theGamma2Dir[2]);

      // The current track is killed by not enqueuing into the next activeQueue.
      break;
    }
    }
  }
}

// Instantiate kernels for electrons and positrons.
__global__ void TransportElectrons(Track *electrons,
				   const adept::MParray *active,
				   Secondaries secondaries,
                                   adept::MParray *activeQueue,
				   adept::MParray *relocateQueue,
                                   GlobalScoring *globalScoring,
				   ScoringPerVolume *scoringPerVolume,
				   HitRecord *hitRecord, 
				   int *eventNumber,
				   ScoringPerParticle *scoringPerParticle,
				   int *mystream_stride
				   )
{
  TransportElectrons</*IsElectron*/ true>(electrons, 
					  active,
					  secondaries,
					  activeQueue,
					  relocateQueue,
					  globalScoring,
                                          scoringPerVolume, 
					  hitRecord, 
					  eventNumber,
					  scoringPerParticle,
					  mystream_stride
					  );
}
__global__ void TransportPositrons(Track *positrons,
				   const adept::MParray *active,
				   Secondaries secondaries,
                                   adept::MParray *activeQueue,
				   adept::MParray *relocateQueue,
                                   GlobalScoring *globalScoring,
				   ScoringPerVolume *scoringPerVolume,
				   HitRecord *hitRecord, 
				   int *eventNumber,
				   ScoringPerParticle *scoringPerParticle,
				   int *mystream_stride
				   )
{
  TransportElectrons</*IsElectron*/ false>(positrons, 
					   active, 
					   secondaries, 
					   activeQueue, 
					   relocateQueue, 
					   globalScoring,
                                           scoringPerVolume,
					   hitRecord, 
					   eventNumber,
					   scoringPerParticle,
					   mystream_stride
					   );
}
