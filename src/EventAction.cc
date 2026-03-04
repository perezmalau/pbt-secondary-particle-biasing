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
/// \file src/EventAction.cc
/// \brief Implementation of the EventAction class, where hits are processed and output ntuples are created

#include "EventAction.hh"
#include "HexitecHit.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include <iomanip>
#include <RunAction.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexitecHitsCollection* EventAction::GetHitsCollection(G4int hcID, const G4Event *event)
{
    auto *hitsCollection = static_cast<HexitecHitsCollection *>(event->GetHCofThisEvent()->GetHC(hcID));

    if (! hitsCollection) {
        G4cerr << "Cannot access hitsCollection ID " << hcID << G4endl;
        exit(1);
    }

    return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    // Print per event (modulo n)
    //
//    auto eventID = event->GetEventID();
//    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
//    if ((printModulo > 0) && (eventID % printModulo == 0)) {
//        G4cout << "--> Starting event: " << eventID << "\n" << G4endl;
//    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event)
{
    G4int evtID = event->GetEventID();
    // Get hits collections IDs (only once)
    if (fSensorHCID == -1) {
        fSensorHCID = G4SDManager::GetSDMpointer()->GetCollectionID("SensorHitsCollection");
    }
    // Get hits collections
    auto sensorHC = GetHitsCollection(fSensorHCID, event);
    G4int nSensor = sensorHC->entries();

    bool copiesHit[3] = {false, false, false};

    for (G4int i = 0; i < nSensor; i++) {
        auto hit = (*sensorHC)[i];
        G4int copyNo = hit->GetCopyNo();
        // Mark the copy number as hit
        if (copyNo >= 0 && copyNo < 3) {
            copiesHit[copyNo] = true;
            // If all 3 copies have been hit, break early
            if (copiesHit[0] && copiesHit[1] && copiesHit[2]) {
                break;
            }
        }
    }

    if (copiesHit[0] && copiesHit[1] && copiesHit[2]) {
        std::map<std::pair<int, int>, G4double> pixelEnergyMap1;
        std::map<std::pair<int, int>, G4double> pixelEnergyMap2;
        std::map<std::pair<int, int>, G4double> pixelEnergyMap3;
        for (G4int i = 0; i < nSensor; i++) {
            auto hit = (*sensorHC)[i];
            G4int parentID = hit->GetParentID();
            G4String particle = hit->GetParticleName();
            G4int copyNo = hit->GetCopyNo();
            G4ThreeVector pos = hit->GetInteractionPos();
            G4String proc = hit->GetProcess();
            auto analysisManager = G4AnalysisManager::Instance();
            // Sensor 1 hits
            if (copyNo == 0) {
                if (parentID == 1 and particle=="gamma") {
                    analysisManager->FillNtupleIColumn(0, 0, evtID);
                    analysisManager->FillNtupleDColumn(0, 1, pos.x()/10);
                    analysisManager->FillNtupleDColumn(0, 2, pos.y()/10);
                    analysisManager->FillNtupleDColumn(0, 3, pos.z()/10);
                    analysisManager->FillNtupleSColumn(0, 4, proc);
                    analysisManager->AddNtupleRow(0);
                    //G4cout << "Gamma created in event " << evtID << " hit Sensor1 and underwent " << proc << G4endl;
                }
                G4int line = hit->GetLineNumber();
                G4int col = hit->GetColNumber();
                auto pixelCoords = std::make_pair(line, col);
                pixelEnergyMap1[pixelCoords] += hit->GetEdep();
            }
            // Sensor 2 hits
            if (copyNo == 1) {
                if (parentID == 1 and particle=="gamma") {
                    analysisManager->FillNtupleIColumn(1, 0, evtID);
                    analysisManager->FillNtupleDColumn(1, 1, pos.x()/10);
                    analysisManager->FillNtupleDColumn(1, 2, pos.y()/10);
                    analysisManager->FillNtupleDColumn(1, 3, pos.z()/10);
                    analysisManager->FillNtupleSColumn(1, 4, proc);
                    analysisManager->AddNtupleRow(1);
                    //G4cout << "Gamma created in event " << evtID << " hit Sensor2 and underwent " << proc << G4endl;
                }
                G4int line = hit->GetLineNumber();
                G4int col = hit->GetColNumber();
                auto pixelCoords = std::make_pair(line, col);
                pixelEnergyMap2[pixelCoords] += hit->GetEdep();
            }
            // Sensor 3 hits
            if (copyNo == 2) {
                if (parentID == 1 and particle=="gamma") {
                    analysisManager->FillNtupleIColumn(2, 0, evtID);
                    analysisManager->FillNtupleDColumn(2, 1, pos.x()/10);
                    analysisManager->FillNtupleDColumn(2, 2, pos.y()/10);
                    analysisManager->FillNtupleDColumn(2, 3, pos.z()/10);
                    analysisManager->FillNtupleSColumn(2, 4, proc);
                    analysisManager->AddNtupleRow(2);
                    //G4cout << "Gamma created in event " << evtID << " hit Sensor3 and underwent " << proc << G4endl;
                }
                G4int line = hit->GetLineNumber();
                G4int col = hit->GetColNumber();
                auto pixelCoords = std::make_pair(line, col);
                pixelEnergyMap3[pixelCoords] += hit->GetEdep();
            }
        }
        // Fill measured information
        auto analysisManager = G4AnalysisManager::Instance();
        // Fill Ntuples for pixel map 1, 2 and 3
        for (const auto &entry: pixelEnergyMap1) {
            G4double etot = entry.second;
            G4int line = entry.first.first;
            G4int col = entry.first.second;
            if (etot > 0.003) {
                analysisManager->FillNtupleIColumn(3, 0, evtID);
                analysisManager->FillNtupleIColumn(3, 1, line);
                analysisManager->FillNtupleIColumn(3, 2, col);
                analysisManager->FillNtupleDColumn(3, 3, etot);
                analysisManager->AddNtupleRow(3);
            }
        }
        for (const auto &entry: pixelEnergyMap2) {
            G4double etot = entry.second;
            G4int line = entry.first.first;
            G4int col = entry.first.second;
            if (etot > 0.003) {
                analysisManager->FillNtupleIColumn(4, 0, evtID);
                analysisManager->FillNtupleIColumn(4, 1, line);
                analysisManager->FillNtupleIColumn(4, 2, col);
                analysisManager->FillNtupleDColumn(4, 3, etot);
                analysisManager->AddNtupleRow(4);
            }
        }
        for (const auto &entry: pixelEnergyMap3) {
            G4double etot = entry.second;
            G4int line = entry.first.first;
            G4int col = entry.first.second;
            if (etot > 0.003) {
                analysisManager->FillNtupleIColumn(5, 0, evtID);
                analysisManager->FillNtupleIColumn(5, 1, line);
                analysisManager->FillNtupleIColumn(5, 2, col);
                analysisManager->FillNtupleDColumn(5, 3, etot);
                analysisManager->AddNtupleRow(5);
            }
        }
    }
}

