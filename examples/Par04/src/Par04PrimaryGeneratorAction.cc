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
#include "Par04PrimaryGeneratorAction.hh"
#include "Par04PrimaryGeneratorMessenger.hh"
#include "Par04DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04PrimaryGeneratorAction::Par04PrimaryGeneratorAction(Par04DetectorConstruction *det)
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0), fDetector(det), fRndmBeam(0.), fGunMessenger(0)
{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);
  SetDefaultKinematic();

  // create a messenger for this class
  fGunMessenger = new Par04PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04PrimaryGeneratorAction::~Par04PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04PrimaryGeneratorAction::GeneratePrimaries(G4Event *aEvent)
{
  // this function is called at the begining of event
  //
  // randomize the beam, if requested.
  if (fRndmBeam > 0.) {
    G4ThreeVector oldPosition = fParticleGun->GetParticlePosition();
    G4double rbeam            = 0.5 * (fDetector->GetCalorSizeYZ()) * fRndmBeam;
    G4double x0               = oldPosition.x();
    G4double y0               = oldPosition.y() + (2 * G4UniformRand() - 1.) * rbeam;
    G4double z0               = oldPosition.z() + (2 * G4UniformRand() - 1.) * rbeam;
    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    fParticleGun->GeneratePrimaryVertex(aEvent);
    fParticleGun->SetParticlePosition(oldPosition);
  } else
    fParticleGun->GeneratePrimaryVertex(aEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04PrimaryGeneratorAction::SetDefaultKinematic()
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition *particle = particleTable->FindParticle(particleName = "e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  fParticleGun->SetParticleEnergy(1. * GeV);
  G4double position = -0.25 * (fDetector->GetWorldSizeX() + fDetector->GetCalorThickness());
  fParticleGun->SetParticlePosition(G4ThreeVector(position, 0. * cm, 0. * cm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Par04PrimaryGeneratorAction::Print() const
{
  std::cout << "=== Gun shooting " << fParticleGun->GetParticleDefinition()->GetParticleName() << " with energy "
            << fParticleGun->GetParticleEnergy() / GeV << "[GeV] from: " << fParticleGun->GetParticlePosition() / mm
            << " [mm]\n";
}
