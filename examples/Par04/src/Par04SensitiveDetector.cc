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
#include "Par04SensitiveDetector.hh"
#include "Par04Hit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"

Par04SensitiveDetector::Par04SensitiveDetector(G4String aName) : G4VSensitiveDetector(aName)
{
  collectionName.insert("hits");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04SensitiveDetector::Par04SensitiveDetector(G4String aName, G4int aNumLayers, G4int aNumAbsorbers)
    : G4VSensitiveDetector(aName), fNumLayers(aNumLayers), fNumAbsorbers(aNumAbsorbers)
{
  collectionName.insert("hits");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04SensitiveDetector::~Par04SensitiveDetector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04SensitiveDetector::Initialize(G4HCofThisEvent *aHCE)
{
  fHitsCollection = new Par04HitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHitCollectionID < 0) {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  aHCE->AddHitsCollection(fHitCollectionID, fHitsCollection);

  // fill calorimeter hits with zero energy deposition
  for (G4int iz = 0; iz < fNumLayers; iz++) {
    for (G4int iz = 0; iz < fNumAbsorbers; iz++) {
      Par04Hit *hit = new Par04Hit();
      fHitsCollection->insert(hit);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return true;

  G4TouchableHistory *aTouchable = (G4TouchableHistory *)(aStep->GetPreStepPoint()->GetTouchable());

  auto hit = RetrieveAndSetupHit(aTouchable);

  // Add energy deposit from G4Step
  hit->AddEdep(edep);

  // Fill time information from G4Step
  // If it's already filled, choose hit with earliest global time
  if (hit->GetTime() == -1 || hit->GetTime() > aStep->GetTrack()->GetGlobalTime())
    hit->SetTime(aStep->GetTrack()->GetGlobalTime());

  // Set hit type to full simulation (only if hit is not already marked as fast
  // sim)
  if (hit->GetType() != 1) hit->SetType(0);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04SensitiveDetector::ProcessHits(const G4FastHit *aHit, const G4FastTrack *aTrack,
                                           G4TouchableHistory *aTouchable)
{
  G4double edep = aHit->GetEnergy();
  if (edep == 0.) return true;

  auto hit = RetrieveAndSetupHit(aTouchable);

  // Add energy deposit from G4FastHit
  hit->AddEdep(edep);

  // Fill time information from G4FastTrack
  // If it's already filled, choose hit with earliest global time
  if (hit->GetTime() == -1 || hit->GetTime() > aTrack->GetPrimaryTrack()->GetGlobalTime()) {
    hit->SetTime(aTrack->GetPrimaryTrack()->GetGlobalTime());
  }

  // Set hit type to fast simulation (even if hit was already marked as full
  // sim, overwrite it)
  hit->SetType(1);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04Hit *Par04SensitiveDetector::RetrieveAndSetupHit(G4TouchableHistory *aTouchable)
{
  G4int layer = aTouchable->GetCopyNumber(1); // layer
  G4int abso  = aTouchable->GetCopyNumber(0); // absorber

  std::size_t hitID = 2 * layer + (abso - 1);

  if (hitID >= fHitsCollection->entries()) {
    G4Exception("Par04SensitiveDetector::RetrieveAndSetupHit()", "InvalidSetup", FatalException,
                "Size of hit collection in Par04SensitiveDetector is smaller than the "
                "number of layers created in Par04DetectorConstruction!");
  }
  Par04Hit *hit = (*fHitsCollection)[hitID];

  if (hit->GetZid() < 0) {
    hit->SetZid(layer);
    hit->SetAbsoid(abso);
    hit->SetLogV(aTouchable->GetVolume(0)->GetLogicalVolume());
    G4AffineTransform transform = aTouchable->GetHistory()->GetTopTransform();
    hit->SetRot(transform.NetRotation());
    transform.Invert();
    hit->SetPos(transform.NetTranslation());
  }
  return hit;
}