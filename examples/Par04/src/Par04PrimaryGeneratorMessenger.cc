// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0
/// \brief Implementation of the PAr04PrimaryGeneratorMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Par04PrimaryGeneratorMessenger.hh"

#include "Par04PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04PrimaryGeneratorMessenger::Par04PrimaryGeneratorMessenger(Par04PrimaryGeneratorAction *Gun)
    : G4UImessenger(), fAction(Gun), fGunDir(0), fDefaultCmd(0), fRndmCmd(0)
{
  fGunDir = new G4UIdirectory("/Par04/gun/");
  fGunDir->SetGuidance("gun control");

  fDefaultCmd = new G4UIcmdWithoutParameter("/Par04/gun/setDefault", this);
  fDefaultCmd->SetGuidance("set/reset kinematic defined in PrimaryGenerator");
  fDefaultCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPrintCmd = new G4UIcmdWithoutParameter("/Par04/gun/print", this);
  fPrintCmd->SetGuidance("print gun kinematics in PrimaryGenerator");
  fPrintCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fRndmCmd = new G4UIcmdWithADouble("/Par04/gun/rndm", this);
  fRndmCmd->SetGuidance("random lateral extension on the beam");
  fRndmCmd->SetGuidance("in fraction of 0.5*sizeYZ");
  fRndmCmd->SetParameterName("rBeam", false);
  fRndmCmd->SetRange("rBeam>=0.&&rBeam<=1.");
  fRndmCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04PrimaryGeneratorMessenger::~Par04PrimaryGeneratorMessenger()
{
  delete fDefaultCmd;
  delete fPrintCmd;
  delete fRndmCmd;
  delete fGunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if (command == fDefaultCmd) {
    fAction->SetDefaultKinematic();
  }

  if (command == fPrintCmd) {
    fAction->Print();
  }

  if (command == fRndmCmd) {
    fAction->SetRndmBeam(fRndmCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
