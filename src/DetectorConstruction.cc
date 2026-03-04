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
//
#include "DetectorConstruction.hh"
#include "CameraBiasingOperator.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "HexitecSD.hh"
#include "G4SDManager.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : fPixelSize(500 * um),
      fSensorThickness(5 * mm),
      fInterLayerSpacing(2.4 * cm),
      fPhantXYZ(20. * cm),
      fPhantCamDistance(20 * cm),
      fNPixel(120) {
    fSensorSizeXY = fNPixel * fPixelSize;
    fCameraXY = fSensorSizeXY;
    fCameraZ = (3 * fSensorThickness) + (2 * fInterLayerSpacing);

    // Compute and print correction factor here — master thread, always visible
    // G4double halfDiagonal = std::sqrt(2.0) * (fCameraXY / 2);
    // G4double halfAngle = std::atan(halfDiagonal / fPhantCamDistance);
    // G4double fCameraSolidAngle = 2.0 * CLHEP::pi * (1.0 - std::cos(halfAngle));
    // G4double correctionFactor = fCameraSolidAngle / (4.0 * CLHEP::pi);

    // G4cout << "\n=== Biasing Parameters ===" << G4endl;
    // G4cout << "Camera solid angle:        " << fCameraSolidAngle / sr << " sr" << G4endl;
    // G4cout << "Correction factor (Ω/4π): " << correctionFactor << G4endl;
    // G4cout << "Variance reduction:        " << 1.0 / correctionFactor << "x" << G4endl;
    // G4cout << "==========================\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
    // Materials:
    G4Material *worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material *phantomMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    auto *elCd = new G4Element("Cadmium", "Cd", 48., 112.414 * g / mole);
    auto *elZn = new G4Element("Zinc", "Zn", 30., 65.38 * g / mole);
    auto *elTe = new G4Element("Tellurium", "Te", 52., 127.6 * g / mole);
    auto *CZT = new G4Material("CZT", 5.8 * g / cm3, 3);
    CZT->AddElementByNumberOfAtoms(elCd, 1);
    CZT->AddElementByNumberOfAtoms(elZn, 1);
    CZT->AddElementByNumberOfAtoms(elTe, 1);
    G4Material *sensorMaterial = G4Material::GetMaterial("CZT");

    // Constants:
    G4double worldSizeXY = 200. * cm;
    G4double worldSizeZ = 300. * cm;


    // -----------------------------------
    // -- World
    // -----------------------------------
    auto *solidWorld = new G4Box("World", worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2);
    auto *logicWorld = new G4LogicalVolume(solidWorld, worldMaterial, "World");
    auto *physWorld = new G4PVPlacement(nullptr, G4ThreeVector(),
                                        logicWorld, "World", nullptr, false, 0);

    // -----------------------------------
    // -- Water phantom (biased volume)
    // -----------------------------------
    auto *solidTest = new G4Box("test.solid", fPhantXYZ / 2, fPhantXYZ / 2, fPhantXYZ / 2);
    auto *logicTest = new G4LogicalVolume(solidTest, phantomMaterial, "test.logical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -fPhantCamDistance),
                      logicTest, "test.phys", logicWorld, false, 0);

    // -----------------------------------
    // -- Compton camera
    // -----------------------------------
    // Main camera envelope
    auto *cameraS = new G4Box("Camera", fCameraXY / 2, fCameraXY / 2, fCameraZ / 2);
    auto *cameraLV = new G4LogicalVolume(cameraS, worldMaterial, "Camera");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fCameraZ / 2),
                      cameraLV, "Camera", logicWorld, false, 0);

    // Detector layers
    auto *layerS = new G4Box("Layer", fSensorSizeXY / 2, fSensorSizeXY / 2, fSensorThickness / 2);
    auto *layerLV = new G4LogicalVolume(layerS, worldMaterial, "Layer");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, -fInterLayerSpacing - fSensorThickness),
                      layerLV, "Layer1", cameraLV, true, 0);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                      layerLV, "Layer2", cameraLV, true, 1);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, fInterLayerSpacing + fSensorThickness),
                      layerLV, "Layer3", cameraLV, true, 2);

    // Pixel and row replicas
    auto *solidPixel = new G4Box("Pixel", fPixelSize / 2, fPixelSize / 2, fSensorThickness / 2);
    auto *logicPixel = new G4LogicalVolume(solidPixel, sensorMaterial, "Pixel");

    auto *solidRow = new G4Box("Row", fSensorSizeXY / 2, fPixelSize / 2, fSensorThickness / 2);
    auto *logicRow = new G4LogicalVolume(solidRow, worldMaterial, "Row");

    new G4PVReplica("Pixel", logicPixel, logicRow, kXAxis, fNPixel, fPixelSize);

    auto *sensorS = new G4Box("Sensor", fSensorSizeXY / 2, fSensorSizeXY / 2, fSensorThickness / 2);
    auto *sensorLV = new G4LogicalVolume(sensorS, worldMaterial, "Sensor");
    new G4PVReplica("Row", logicRow, sensorLV, kYAxis, fNPixel, fPixelSize);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), sensorLV, "Sensor", layerLV, false, 0);


    auto *phantVisAtt = new G4VisAttributes(G4Colour(0.2, 0.6, 0.9));
    phantVisAtt->SetVisibility(true);
    logicTest->SetVisAttributes(phantVisAtt);

    auto worldVisAttr = new G4VisAttributes();
    worldVisAttr->SetVisibility(false);
    logicWorld->SetVisAttributes(worldVisAttr);

    sensorLV->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
    logicRow->SetVisAttributes(new G4VisAttributes(false));
    logicPixel->SetVisAttributes(new G4VisAttributes(false));

    return physWorld;
}

void DetectorConstruction::ConstructSDandField() {
    // Compton camera sensitive detector
    auto *sensorSD = new HexitecSD("SensorSD", "SensorHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(sensorSD);
    SetSensitiveDetector("Pixel", sensorSD);

    // Biasing: force secondary gammas from the phantom toward the camera
    G4ThreeVector cameraCenter(0, 0, 0);
    G4LogicalVolume *logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume("test.logical");
    auto *biasingOp = new CameraBiasingOperator(cameraCenter, fCameraXY/2);
    biasingOp->AttachTo(logicTest);
}
