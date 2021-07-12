// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0
//
/// \brief Definition of the Par04PrimaryGeneratorMessenger class
//
// $Id: PrimaryGeneratorMessenger.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Par04PrimaryGeneratorMessenger_h
#define Par04PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Par04PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Par04PrimaryGeneratorMessenger : public G4UImessenger {
public:
  Par04PrimaryGeneratorMessenger(Par04PrimaryGeneratorAction *);
  ~Par04PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand *, G4String);

private:
  Par04PrimaryGeneratorAction *fAction;

  G4UIdirectory *fGunDir;
  G4UIcmdWithoutParameter *fDefaultCmd;
  G4UIcmdWithoutParameter *fPrintCmd;
  G4UIcmdWithADouble *fRndmCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
