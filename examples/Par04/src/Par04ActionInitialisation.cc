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
#include "Par04ActionInitialisation.hh"
#include "Par04PrimaryGeneratorAction.hh"
#include "Par04EventAction.hh"
#include "Par04RunAction.hh"

Par04ActionInitialisation::Par04ActionInitialisation(Par04DetectorConstruction *aDetector)
    : G4VUserActionInitialization(), fDetector(aDetector)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ActionInitialisation::~Par04ActionInitialisation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ActionInitialisation::BuildForMaster() const
{
  SetUserAction(new Par04RunAction(fDetector));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ActionInitialisation::Build() const
{
  SetUserAction(new Par04PrimaryGeneratorAction(fDetector));
  SetUserAction(new Par04EventAction(fDetector));
  SetUserAction(new Par04RunAction(fDetector));
}
