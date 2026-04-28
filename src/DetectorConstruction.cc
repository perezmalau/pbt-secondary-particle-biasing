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
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
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
      fSensorThickness(2 * mm),
      fInterLayerSpacing(2.4 * cm),
      fPhantXYZ(20. * cm),
      fPhantCamDistance(20 * cm),
      fNPixel(120) {
    fSensorSizeXY = fNPixel * fPixelSize;
    fCameraXY = fSensorSizeXY;
    fCameraZ = (3 * fSensorThickness) + (2 * fInterLayerSpacing);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
    // Materials:
    G4NistManager *nist = G4NistManager::Instance();
    G4Material *worldMaterial = nist->FindOrBuildMaterial("G4_AIR");
    // G4Material *phantomMaterial = nist->FindOrBuildMaterial("G4_WATER");
    G4Material *skullMaterial = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
    G4Material *skinMaterial = nist->FindOrBuildMaterial("G4_SKIN_ICRP");
    G4Material* tissueEquivalentMaterial = new G4Material("TissueEquivalent", 1.0*g/cm3, 4);

    tissueEquivalentMaterial->AddElement(nist->FindOrBuildElement("H"), 0.101);
    tissueEquivalentMaterial->AddElement(nist->FindOrBuildElement("C"), 0.111);
    tissueEquivalentMaterial->AddElement(nist->FindOrBuildElement("N"), 0.026);
    tissueEquivalentMaterial->AddElement(nist->FindOrBuildElement("O"), 0.762);

    // Set parts per million of boron
    auto ppm_low = 25e-6;
    auto ppm_high = 62.5e-6;

    auto* B10 = new G4Isotope("Boron10", 5, 10, 10. * g / mole);
    auto* elB10 = new G4Element("Boron10", "B10", 1);
    elB10->AddIsotope(B10, 1.0); // Abundance fraction = 100%

    auto* boronMaterial = new G4Material("boronMaterial", 2.37 * g / cm3, 1);
    boronMaterial->AddElement(elB10, 1.0); // 100% Boron-10

    G4Material* boronatedTissue = new G4Material("boronatedTissue", 1.0*g/cm3, 2);
    boronatedTissue->AddMaterial(tissueEquivalentMaterial, 1.0 - ppm_low);
    boronatedTissue->AddMaterial(boronMaterial, ppm_low);

    G4Material* boronatedTumour = new G4Material("boronatedTumour", 1.0*g/cm3, 2);
    boronatedTumour->AddMaterial(tissueEquivalentMaterial, 1.0 - ppm_high);
    boronatedTumour->AddMaterial(boronMaterial, ppm_high);

    // CZT
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
    // auto *solidTest = new G4Box("test.solid", fPhantXYZ / 2, fPhantXYZ / 2, fPhantXYZ / 2);
    // auto *logicTest = new G4LogicalVolume(solidTest, phantomMaterial, "test.logical");
    // new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -fPhantCamDistance),
    //                   logicTest, "test.phys", logicWorld, false, 0);

    // -----------------------------------
    // -- Snyder phantom (biased volume)
    // -----------------------------------
    G4double ax_brain = 6. *cm;
    G4double by_brain = 9. *cm;
    G4double cz_brain = 6.5 *cm;
    G4double ax_skull = 6.8 *cm;
    G4double by_skull = 9.8 *cm;
    G4double cz_skull = 8.3 *cm;
    G4double ax_head = 7.0 * cm;
    G4double by_head = 10.0 * cm;
    G4double cz_head = 8.50 * cm;
    G4double tumour_radius = 1 *cm;
    G4double tumour_posX = -3.5 *cm;

    // Head
    auto* headS = new G4Ellipsoid("Head", ax_head, by_head, cz_head);
    auto *headLV = new G4LogicalVolume(headS, skinMaterial, "Head");
    new G4PVPlacement(nullptr,G4ThreeVector(0, 0, -fPhantCamDistance),
        headLV,"Head", logicWorld, false, 0);

    // Skull
    auto* craniumOut = new G4Ellipsoid("CraniumOut", ax_skull, by_skull, cz_skull);
    auto* craniumIn = new G4Ellipsoid("CraniumIn", ax_brain, by_brain, cz_brain);
    auto* skullS = new G4SubtractionSolid("Skull", craniumOut, craniumIn, nullptr,
                                         G4ThreeVector(0.0, 0.0, 1. *cm));
    auto *skullLV = new G4LogicalVolume(skullS, skullMaterial, "Skull");
    new G4PVPlacement(nullptr,G4ThreeVector(0., 0., 0.),
                      skullLV, "Skull",  headLV, false, 0);

	// Brain
    auto* brainS = new G4Ellipsoid("Brain", ax_brain, by_brain, cz_brain);
    auto *brainLV = new G4LogicalVolume(brainS, boronatedTissue, "Brain");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 1 *cm),
        brainLV, "Brain", headLV, false, 0);

    // Tumour
    auto *targetS = new G4Sphere("Target", 0, tumour_radius,
        0, 2*M_PI, 0, M_PI);
    auto *targetLV = new G4LogicalVolume(targetS, boronatedTumour, "Target");
    new G4PVPlacement(nullptr, G4ThreeVector(tumour_posX, 0 *cm, 0. *cm),
        targetLV,"Target", brainLV, false, 0);

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
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
        sensorLV, "Sensor", layerLV, false, 0);


    // auto *phantVisAtt = new G4VisAttributes(G4Colour(0.2, 0.6, 0.9));
    // phantVisAtt->SetVisibility(true);
    // logicTest->SetVisAttributes(phantVisAtt);

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
    //G4LogicalVolume *logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume("test.logical");
    G4LogicalVolume *logicHead = G4LogicalVolumeStore::GetInstance()->GetVolume("Head");
    G4LogicalVolume *logicSkull = G4LogicalVolumeStore::GetInstance()->GetVolume("Skull");
    G4LogicalVolume *logicBrain = G4LogicalVolumeStore::GetInstance()->GetVolume("Brain");
    G4LogicalVolume *logicTarget = G4LogicalVolumeStore::GetInstance()->GetVolume("Target");
    auto *biasingOp = new CameraBiasingOperator(cameraCenter, fCameraXY/2);
    // biasingOp->AttachTo(logicTest);
    biasingOp->AttachTo(logicHead);
    biasingOp->AttachTo(logicSkull);
    biasingOp->AttachTo(logicBrain);
    biasingOp->AttachTo(logicTarget);
}
