// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "Par04EMShowerMessenger.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4FastHit.hh"
#include "Randomize.hh"
#include "G4FastSimHitMaker.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#include "AdeptIntegration.h"

#include "Par04EMShowerModel.hh"

#include <VecGeom/base/Config.h>
#include <VecGeom/management/GeoManager.h>

Par04EMShowerModel::Par04EMShowerModel(G4String aModelName, G4Region *aEnvelope)
    : G4VFastSimulationModel(aModelName, aEnvelope), fMessenger(new Par04EMShowerMessenger(this)),
      fHitMaker(new G4FastSimHitMaker)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EMShowerModel::Par04EMShowerModel(G4String aModelName)
    : G4VFastSimulationModel(aModelName), fMessenger(new Par04EMShowerMessenger(this)), fHitMaker(new G4FastSimHitMaker)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EMShowerModel::~Par04EMShowerModel() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04EMShowerModel::IsApplicable(const G4ParticleDefinition &aParticleType)
{
  return &aParticleType == G4Electron::ElectronDefinition() || &aParticleType == G4Positron::PositronDefinition() ||
         &aParticleType == G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04EMShowerModel::ModelTrigger(const G4FastTrack &aFastTrack)
{
  /*
  // Check energy
  if(aFastTrack.GetPrimaryTrack()->GetKineticEnergy() < 1 * GeV)
  {
    return false;
  }
  // Check length of detector
  // Calculate depth of the detector along shower axis to verify if shower
  // will fit inside. Required max shower depth is defined by fLongMaxDepth, and
  // can be changed with UI command `/Par04/fastSim/longitudinalProfile/maxDepth
  G4double X0 = aFastTrack.GetPrimaryTrack()->GetMaterial()->GetRadlen();
  auto particleDirection     = aFastTrack.GetPrimaryTrackLocalDirection();
  auto particlePosition      = aFastTrack.GetPrimaryTrackLocalPosition();
  G4double detectorDepthInMM = aFastTrack.GetEnvelopeSolid()->DistanceToOut(
    particlePosition, particleDirection);
  G4double detectorDepthInX0 = detectorDepthInMM / X0;
 */
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EMShowerModel::DoIt(const G4FastTrack &aFastTrack, G4FastStep &aFastStep)
{
  // Remove particle from further processing by G4
  aFastStep.KillPrimaryTrack();
  aFastStep.SetPrimaryTrackPathLength(0.0);
  G4double energy = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();
  // No need to create any deposit, it will be handled by this model (and
  // G4FastSimHitMaker that will call the sensitive detector)
  aFastStep.SetTotalEnergyDeposited(0);
  auto particlePosition  = aFastTrack.GetPrimaryTrackLocalPosition();
  auto particleDirection = aFastTrack.GetPrimaryTrackLocalDirection();

  auto pdg = aFastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGEncoding();

  int tid = G4Threading::G4GetThreadId();
  if (tid < 0) tid = 0;

  std::cout << "Thread " << tid << " calling AdePT for particle " << pdg << " energy " << energy << " position "
            << particlePosition[0] << " " << particlePosition[1] << " " << particlePosition[2] << " direction "
            << particleDirection[0] << " " << particleDirection[1] << " " << particleDirection[2] << std::endl;

  AdeptIntegration::Instance().AddTrack(tid, pdg, energy, particlePosition[0], particlePosition[1], particlePosition[2],
                                        particleDirection[0], particleDirection[1], particleDirection[2]);

  // I need to pass the particle from Geant4 to AdePT and simulate the shower
  AdeptIntegration::Instance().Shower(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID(), tid);

  // Create energy deposit in the detector

  for (auto id = 0; id != NumVolumes; id++) {
    // std::cout << " ID " << id << " Charged-TrakL " <<
    // AdeptIntegration::Instance().fUserData[tid].scoringPerVolume.chargedTrackLength[id] / copcore::units::mm
    //          << " mm; Energy-Dep " << AdeptIntegration::Instance().fUserData[tid].scoringPerVolume.energyDeposit[id]
    //          / copcore::units::MeV << " MeV" << std::endl;
    fHitMaker->make(
        G4FastHit(G4ThreeVector(id, 0, 0),
                  AdeptIntegration::Instance().fUserData[tid].scoringPerVolume.energyDeposit[id] / copcore::units::MeV),
        aFastTrack);
  }

  //      generatedHits++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EMShowerModel::Print() const
{
  G4cout << "Par04EMShowerModel: " << G4endl;
}

void Par04EMShowerModel::Initialize()
{
  // This is supposed to set the max batching for Adept to allocate properly the memory
  AdeptIntegration::Instance().SetMaxBatch(25);
  AdeptIntegration::Instance().Initialize();
}
