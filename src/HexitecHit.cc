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
/// \file B4/B4c/src/HexitecHit.cc
/// \brief Implementation of the HexitecHit class

#include "HexitecHit.hh"

#include "G4UnitsTable.hh"


G4ThreadLocal G4Allocator<HexitecHit> *HexitecHitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexitecHit::HexitecHit() : G4VHit(),
                           fCopyNo(-1),
                           fTrackID(0),
                           fParticleName(""),
                           fProcessName(""),
                           fInteractionPos(0),
                           fEdiff(0.),
                           fLine(-1),
                           fCol(-1),
                           fEdep(0.) {
    // fTrueAngle(0.),
    // fTrueEnergy(0.){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexitecHit::~HexitecHit() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexitecHit::HexitecHit(const HexitecHit &right)
    : G4VHit() {
    fCopyNo = right.fCopyNo;
    fTrackID = right.fTrackID;
    fParticleName = right.fParticleName;
    fProcessName = right.fProcessName;
    fInteractionPos = right.fInteractionPos;
    fEdiff = right.fEdiff;
    fLine = right.fLine;
    fCol = right.fCol;
    fEdep = right.fEdep;
    // fTrueAngle = right.fTrueAngle;
    // fTrueEnergy = right.fTrueEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const HexitecHit &HexitecHit::operator=(const HexitecHit &right) {
    fCopyNo = right.fCopyNo;
    fTrackID = right.fTrackID;
    fParticleName = right.fParticleName;
    fInteractionPos = right.fInteractionPos;
    fEdiff = right.fEdiff;
    fLine = right.fLine;
    fCol = right.fCol;
    fEdep = right.fEdep;
    // fTrueAngle = right.fTrueAngle;
    // fTrueEnergy = right.fTrueEnergy;

    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int HexitecHit::operator==(const HexitecHit &right) const {
    return (this == &right) ? 1 : 0;
}
