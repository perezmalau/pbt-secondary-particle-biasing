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
/// \file src/HexitecSD.cc

#include "HexitecSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexitecSD::HexitecSD(const G4String &name, const G4String &collName)
    : G4VSensitiveDetector(name)
{
    collectionName.insert(collName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HexitecSD::Initialize(G4HCofThisEvent *hce)
{
    // Create hits collection
    fHitsCollection = new HexitecHitsCollection(SensitiveDetectorName, collectionName[0]);
    // Add this collection in hce
    auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HexitecSD::ProcessHits(G4Step *step, G4TouchableHistory*)
{
    G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    if (proc!="CoupledTransportation" and proc!="Rayl") {
        G4int parentID = step->GetTrack()->GetParentID();
        G4int trackID = step->GetTrack()->GetTrackID();
        G4String particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
        G4double edep = step->GetTotalEnergyDeposit();
        G4double epre = step->GetPreStepPoint()->GetKineticEnergy();
        G4double epost = step->GetPostStepPoint()->GetKineticEnergy();
        G4double ediff = epre - epost;
        G4int lineNumber = step->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);
        G4int columnNumber = step->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0);
        G4int copyNo = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber(3);
        G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
        G4double init_ene = step->GetTrack()->GetVertexKineticEnergy();
        G4double weight = step->GetTrack()->GetWeight();
        auto *hit = new HexitecHit();
        hit->Add(copyNo, parentID, trackID, edep, lineNumber, columnNumber, pos, init_ene, ediff, particle, weight);
        fHitsCollection->insert(hit);
    }
    return true;
}
