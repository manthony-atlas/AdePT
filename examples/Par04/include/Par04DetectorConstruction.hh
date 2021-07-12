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
#ifndef PAR03DETECTORCONSTRUCTION_H
#define PAR03DETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class G4UniformMagField;
class G4FieldManager;
class Par04PrimaryGeneratorAction;

class Par04DetectorMessenger;
const G4int kMaxAbsor = 10; // 0 + 9

/**
 * @brief Detector construction.
 *
 * Sensitive detector Par04SensitiveDetector is attached to the
 * Absorber volumes.
 * Region for the detector is created as an envelope of the fast simulation.
 *
 */

class Par04DetectorConstruction : public G4VUserDetectorConstruction {
public:
  Par04DetectorConstruction();
  virtual ~Par04DetectorConstruction();

  virtual G4VPhysicalVolume *Construct() final;
  void CreateVecGeomWorld();
  virtual void ConstructSDandField() final;

  void SetNbOfLayers(G4int);
  void SetNbOfAbsor(G4int);
  void SetAbsorMaterial(G4int, const G4String &);
  void SetAbsorThickness(G4int, G4double);

  void SetWorldMaterial(const G4String &);
  void SetCalorSizeYZ(G4double);

  G4double GetWorldSizeX() { return fWorldSizeX; }
  G4double GetWorldSizeYZ() { return fWorldSizeYZ; }

  G4double GetCalorThickness() { return fCalorThickness; }
  G4double GetCalorSizeYZ() { return fCalorSizeYZ; }

  G4int GetNbOfLayers() { return fNbOfLayers; }

  G4int GetNbOfAbsor() { return fNbOfAbsor; }
  G4Material *GetAbsorMaterial(G4int i) { return fAbsorMaterial[i]; }
  G4double GetAbsorThickness(G4int i) { return fAbsorThickness[i]; }

  const G4VPhysicalVolume *GetphysiWorld() { return fPhysiWorld; }
  const G4Material *GetWorldMaterial() { return fDefaultMaterial; }
  const G4VPhysicalVolume *GetAbsorber(G4int i) { return fPhysiAbsor[i]; }

  // Set uniform magnetic field
  inline void SetMagField(const G4ThreeVector &fv) { fMagFieldVector = fv; }
  void SetPrimaryGenerator(Par04PrimaryGeneratorAction *pg) { fPrimaryGenerator = pg; }

  // Print detector information
  void Print() const;

private:
  /// Messenger that allows to modify geometry
  Par04DetectorMessenger *fDetectorMessenger;

  G4int fNbOfAbsor = 2;
  G4Material *fAbsorMaterial[kMaxAbsor];
  G4double fAbsorThickness[kMaxAbsor];

  /// Number of layers = slices along z axis
  G4int fNbOfLayers = 50;
  G4double fLayerThickness;

  G4double fCalorSizeYZ = 40 * cm;
  G4double fCalorThickness;

  G4Material *fDefaultMaterial;
  G4double fWorldSizeYZ;
  G4double fWorldSizeX;

  G4Box *fSolidWorld;
  G4LogicalVolume *fLogicWorld;
  G4VPhysicalVolume *fPhysiWorld;

  G4Box *fSolidCalor;
  G4LogicalVolume *fLogicCalor;
  G4VPhysicalVolume *fPhysiCalor;

  G4Box *fSolidLayer;
  G4LogicalVolume *fLogicLayer;
  G4VPhysicalVolume *fPhysiLayer;

  G4Box *fSolidAbsor[kMaxAbsor];
  G4LogicalVolume *fLogicAbsor[kMaxAbsor];
  G4VPhysicalVolume *fPhysiAbsor[kMaxAbsor];

  // field related members
  G4ThreeVector fMagFieldVector;
  Par04PrimaryGeneratorAction *fPrimaryGenerator;
};

#endif /* PAR03DETECTORCONSTRUCTION_H */
