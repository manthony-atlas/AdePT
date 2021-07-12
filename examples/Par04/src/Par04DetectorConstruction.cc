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
#include "Par04DetectorConstruction.hh"
#include "Par04DetectorMessenger.hh"
#include "Par04SensitiveDetector.hh"
#include "Par04EMShowerModel.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4SDManager.hh"

#include "G4UnitsTable.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include <G4ProductionCutsTable.hh>

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include <VecGeom/base/Config.h>
#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/volumes/UnplacedBox.h>

#include "Scoring.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorConstruction::Par04DetectorConstruction() : G4VUserDetectorConstruction()
{
  fDetectorMessenger = new Par04DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorConstruction::~Par04DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *Par04DetectorConstruction::Construct()
{
  // Compute derived parameters of the calorimeter
  fLayerThickness = 0.;
  for (G4int iAbs = 1; iAbs <= fNbOfAbsor; iAbs++) {
    fLayerThickness += fAbsorThickness[iAbs];
  }
  fCalorThickness = fNbOfLayers * fLayerThickness;
  fWorldSizeX     = 1.2 * fCalorThickness;
  fWorldSizeYZ    = 1.2 * fCalorSizeYZ;
  // update the primary generator
  if (fPrimaryGenerator) {
    fPrimaryGenerator->SetDefaultKinematic();
  }

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //--------- Material definition ---------
  const char *WorldMaterial    = "G4_Galactic";
  const char *GapMaterial      = "G4_Pb";
  const char *AbsorberMaterial = "G4_lAr";

  G4Material *default_mat  = G4NistManager::Instance()->FindOrBuildMaterial(WorldMaterial);
  G4Material *gap_mat      = G4NistManager::Instance()->FindOrBuildMaterial(GapMaterial);
  G4Material *absorber_mat = G4NistManager::Instance()->FindOrBuildMaterial(AbsorberMaterial);

  fSolidWorld = new G4Box("World", fWorldSizeX / 2., fWorldSizeYZ / 2., fWorldSizeYZ / 2.);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, default_mat, "World");
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0);

  fSolidCalor = new G4Box("Calorimeter", fCalorThickness / 2., fCalorSizeYZ / 2., fCalorSizeYZ / 2.);
  fLogicCalor = new G4LogicalVolume(fSolidCalor, default_mat, "Calorimeter");
  fPhysiCalor = new G4PVPlacement(0, G4ThreeVector(), fLogicCalor, "Calorimeter", fLogicWorld, false, 0);

  //
  // Layers
  //

  fSolidLayer = new G4Box("Layer", fLayerThickness / 2, fCalorSizeYZ / 2, fCalorSizeYZ / 2);
  fLogicLayer = new G4LogicalVolume(fSolidLayer, default_mat, "Layer");
  fPhysiLayer = new G4PVReplica("Layer", fLogicLayer, fLogicCalor, kXAxis, fNbOfLayers, fLayerThickness);

  //
  // Absorbers
  //
  G4double xfront = -0.5 * fLayerThickness;
  for (G4int k = 1; k <= fNbOfAbsor; k++) {
    fSolidAbsor[k] = new G4Box("Absorber", // its name
                               fAbsorThickness[k] / 2, fCalorSizeYZ / 2, fCalorSizeYZ / 2);

    fLogicAbsor[k] = new G4LogicalVolume(fSolidAbsor[k],    // its solid
                                         fAbsorMaterial[k], // its material
                                         fAbsorMaterial[k]->GetName());

    G4double xcenter = xfront + 0.5 * fAbsorThickness[k];
    xfront += fAbsorThickness[k];
    fPhysiAbsor[k] = new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicAbsor[k], fAbsorMaterial[k]->GetName(),
                                       fLogicLayer, false,
                                       k); // copy number
  }

  // Region for fast simulation
  auto detectorRegion = new G4Region("DetectorRegion");
  detectorRegion->AddRootLogicalVolume(fLogicCalor);
  detectorRegion->UsedInMassGeometry(true);

  G4ProductionCuts *productionCuts = new G4ProductionCuts();
  productionCuts->SetProductionCut(ProductionCut);
  //
  detectorRegion->SetProductionCuts(productionCuts);
  //

  G4ProductionCutsTable *theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  theCoupleTable->UpdateCoupleTable(fPhysiWorld);

  // Print();
  CreateVecGeomWorld();
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::ConstructSDandField()
{
  if (fMagFieldVector.mag() > 0.0) {
    // Apply a global uniform magnetic field along the Z axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.
    auto uniformMagField     = new G4UniformMagField(fMagFieldVector);
    G4FieldManager *fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(uniformMagField);
    fieldMgr->CreateChordFinder(uniformMagField);
    G4cout << G4endl << " *** SETTING MAGNETIC FIELD : fieldValue = " << fMagFieldVector / kilogauss
           << " [kilogauss] *** " << G4endl << G4endl;

  } else {
    G4cout << G4endl << " *** NO MAGNETIC FIELD SET  *** " << G4endl << G4endl;
  }

  Par04SensitiveDetector *caloSD = new Par04SensitiveDetector("sensitiveDetector", fNbOfLayers, fNbOfAbsor);
  G4SDManager::GetSDMpointer()->AddNewDetector(caloSD);
  for (G4int k = 1; k <= fNbOfAbsor; k++) {
    SetSensitiveDetector(fLogicAbsor[k], caloSD);
  }

  auto detectorRegion             = G4RegionStore::GetInstance()->GetRegion("DetectorRegion");
  Par04EMShowerModel *showermodel = new Par04EMShowerModel("model", detectorRegion);
  showermodel->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::Print() const
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << fNbOfLayers << " layers of:";
  for (G4int i = 1; i <= fNbOfAbsor; i++) {
    G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() << ": " << std::setw(6)
           << G4BestUnit(fAbsorThickness[i], "Length");
  }
  G4cout << "\n-------------------------------------------------------------\n";

  // G4cout << "\n" << fDefaultMaterial << G4endl;
  for (G4int j = 1; j <= fNbOfAbsor; j++)
    G4cout << "\n" << fAbsorMaterial[j] << G4endl;

  G4cout << "\n-------------------------------------------------------------\n";
}

void Par04DetectorConstruction::CreateVecGeomWorld()
{
  auto worldSolid = new vecgeom::UnplacedBox(0.5 * WorldSizeX, 0.5 * WorldSizeYZ, 0.5 * WorldSizeYZ);
  auto worldLogic = new vecgeom::LogicalVolume("World", worldSolid);
  vecgeom::VPlacedVolume *worldPlaced = worldLogic->Place();

  //
  // Calorimeter
  //
  auto calorSolid = new vecgeom::UnplacedBox(0.5 * CalorThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto calorLogic = new vecgeom::LogicalVolume("Calorimeter", calorSolid);
  vecgeom::Transformation3D origin;
  worldLogic->PlaceDaughter(calorLogic, &origin);

  //
  // Layers
  //
  auto layerSolid = new vecgeom::UnplacedBox(0.5 * LayerThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);

  //
  // Absorbers
  //
  auto gapSolid = new vecgeom::UnplacedBox(0.5 * GapThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto gapLogic = new vecgeom::LogicalVolume("Gap", gapSolid);
  vecgeom::Transformation3D gapPlacement(-0.5 * LayerThickness + 0.5 * GapThickness, 0, 0);

  auto absorberSolid = new vecgeom::UnplacedBox(0.5 * AbsorberThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto absorberLogic = new vecgeom::LogicalVolume("Absorber", absorberSolid);
  vecgeom::Transformation3D absorberPlacement(0.5 * LayerThickness - 0.5 * AbsorberThickness, 0, 0);

  // Create a new LogicalVolume per layer, we need unique IDs for scoring.
  double xCenter = -0.5 * CalorThickness + 0.5 * LayerThickness;
  for (int i = 0; i < NbOfLayers; i++) {
    auto layerLogic = new vecgeom::LogicalVolume("Layer", layerSolid);
    vecgeom::Transformation3D placement(xCenter, 0, 0);
    calorLogic->PlaceDaughter(layerLogic, &placement);

    layerLogic->PlaceDaughter(gapLogic, &gapPlacement);
    layerLogic->PlaceDaughter(absorberLogic, &absorberPlacement);

    xCenter += LayerThickness;
  }

  vecgeom::GeoManager::Instance().SetWorldAndClose(worldPlaced);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::SetWorldMaterial(const G4String &material)
{
  // search the material by its name
  G4Material *pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) fDefaultMaterial = pttoMaterial;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::SetNbOfLayers(G4int ival)
{
  // set the number of Layers
  //
  if (ival < 1) {
    G4cout << "\n --->warning from SetfNbOfLayers: " << ival << " must be at least 1. Command refused" << G4endl;
    return;
  }
  fNbOfLayers = ival;
  // G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor - 1)) {
    G4cout << "\n ---> warning from SetfNbOfAbsor: " << ival << " must be at least 1 and and most " << kMaxAbsor - 1
           << ". Command refused" << G4endl;
    return;
  }
  fNbOfAbsor = ival;
  // G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::SetAbsorMaterial(G4int ival, const G4String &material)
{
  // search the material by its name
  //
  if (ival > fNbOfAbsor || ival <= 0) {
    G4cout << "\n --->warning from SetAbsorMaterial: absor number " << ival << " out of range. Command refused"
           << G4endl;
    return;
  }

  G4Material *pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) fAbsorMaterial[ival] = pttoMaterial;
  // G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::SetAbsorThickness(G4int ival, G4double val)
{
  // change Absorber thickness
  //
  if (ival > fNbOfAbsor || ival <= 0) {
    G4cout << "\n --->warning from SetAbsorThickness: absor number " << ival << " out of range. Command refused"
           << G4endl;
    return;
  }
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetAbsorThickness: thickness " << val << " out of range. Command refused" << G4endl;
    return;
  }
  fAbsorThickness[ival] = val;
  // G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetfCalorSizeYZ: thickness " << val << " out of range. Command refused" << G4endl;
    return;
  }
  fCalorSizeYZ = val;
  // G4RunManager::GetRunManager()->ReinitializeGeometry();
}
