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
#ifndef PAR03EMSHOWERMODEL_HH
#define PAR03EMSHOWERMODEL_HH

#include "G4VFastSimulationModel.hh"
//#include <AdePT/ArgParser.h>
#include <CopCore/SystemOfUnits.h>

class Par04EMShowerMessenger;
class G4FastSimHitMaker;

/**
 * @brief Example fast simulation model for EM showers.
 *
 * Parametrisation of electrons, positrons, and gammas. It is triggered if those
 * particles enter the detector so that there is sufficient length for the
 * shower development (max depth, controlled by the UI command).
 *
 * Using AdePT
 */

class Par04EMShowerModel : public G4VFastSimulationModel {
public:
  Par04EMShowerModel(G4String, G4Region *);
  Par04EMShowerModel(G4String);
  ~Par04EMShowerModel();

  /// There are no kinematics constraints. True is returned.
  virtual G4bool ModelTrigger(const G4FastTrack &) final;
  /// Model is applicable to electrons, positrons, and photons.
  virtual G4bool IsApplicable(const G4ParticleDefinition &) final;

  /// Take particle out of the full simulation (kill it at the entrance
  /// depositing all the energy). Simulate the full shower using AdePT library.
  virtual void DoIt(const G4FastTrack &, G4FastStep &) final;

  /// Print current settings.
  void Print() const;

  /// Initialized VecGeom, etc.

  void Initialize();

private:
  /// Messenger for configuration
  Par04EMShowerMessenger *fMessenger;
  /// Helper class for creation of hits within the sensitive detector
  std::unique_ptr<G4FastSimHitMaker> fHitMaker;

  G4double ProductionCut = 0.7 * copcore::units::mm;

  double CalorSizeYZ       = 40 * copcore::units::cm;
  int NbOfLayers           = 50;
  int NbOfAbsorbers        = 2;
  double GapThickness      = 2.3 * copcore::units::mm;
  double AbsorberThickness = 5.7 * copcore::units::mm;

  double LayerThickness = GapThickness + AbsorberThickness;
  double CalorThickness = NbOfLayers * LayerThickness;

  double WorldSizeX  = 1.2 * CalorThickness;
  double WorldSizeYZ = 1.2 * CalorSizeYZ;

  int NumVolumes = 1 + 1 + NbOfLayers * (1 + NbOfAbsorbers);
  int MCIndex[100];
};
#endif /* PAR03EMSHOWERMODEL_HH */