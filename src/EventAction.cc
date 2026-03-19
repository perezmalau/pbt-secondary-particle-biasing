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
    G4int nHits = sensorHC->entries();

    // --- Group hits by trackID, preserving all hits per stage ---
    std::map<G4int, std::array<std::vector<HexitecHit*>, 3>> hitsByTrack;
    for (G4int i = 0; i < nHits; i++)
    {
        HexitecHit* hit = (*sensorHC)[i];
        G4int copyNo = hit->GetCopyNo();
        G4String particle = hit->GetParticleName();
        if (copyNo >= 0 && copyNo <= 2)
            hitsByTrack[hit->GetTrackID()][copyNo].push_back(hit);
    }

    // Define pixel maps - the realistic detector information
    // std::map<std::pair<int,int>, G4double> pixelMaps[3];

    // --- Find coincident gammas: must have at least one hit in each stage ---
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    for (auto& [trackID, hitsPerStage] : hitsByTrack)
    {
        // Condition for triple coincidence:
        // bool isCoincident = !hitsPerStage[0].empty()
        //                  && !hitsPerStage[1].empty()
        //                  && !hitsPerStage[2].empty();
        // if (!isCoincident) continue;

        // Condition for any kind of coincidence involving 2 or more stages:
        G4int stagesHit = (G4int)!hitsPerStage[0].empty()
                + (G4int)!hitsPerStage[1].empty()
                + (G4int)!hitsPerStage[2].empty();
        if (stagesHit < 2) continue;

        // --- Record every individual hit across all three stages ---
        for (G4int stage = 0; stage < 3; stage++)
        {
            for (HexitecHit* hit : hitsPerStage[stage])
            {
                auto pos = hit->GetInteractionPos();
                auto ediff = hit->GetEnergyDiff();
                auto process = hit->GetProcessName();
                auto particle = hit->GetParticleName();
                if (particle == "gamma") {
                    analysisManager->FillNtupleIColumn(stage, 0, evtID);
                    analysisManager->FillNtupleIColumn(stage, 1, trackID);
                    analysisManager->FillNtupleDColumn(stage, 2, pos.x()/CLHEP::mm);
                    analysisManager->FillNtupleDColumn(stage, 3, pos.y()/CLHEP::mm);
                    analysisManager->FillNtupleDColumn(stage, 4, pos.z()/CLHEP::mm);
                    analysisManager->FillNtupleDColumn(stage, 5, ediff);
                    analysisManager->FillNtupleSColumn(stage, 6, process);
                    analysisManager->AddNtupleRow(stage);
                }
            }
        }
    }


    // Loop through all hits in the event, record the relevant information
    // for (G4int i = 0; i < nHits; i++) {
    //     auto hit      = (*sensorHC)[i];
    //     G4int copyNo  = hit->GetCopyNo();
    //     // Accumulate pixel energy
    //     G4int line = hit->GetLineNumber();
    //     G4int col  = hit->GetColNumber();
    //     pixelMaps[copyNo][{line, col}] += hit->GetEdep();
    // }

    // // Fill realistic detector information based on event coincidences and a detector threshold at 3 keV
    // constexpr G4double kEnergyThreshold = 0.003; // MeV
    //
    // // Realistic coincidence: all 3 pixel maps are non-empty
    // bool realisticCoincidence = !pixelMaps[0].empty() &&
    //                             !pixelMaps[1].empty() &&
    //                             !pixelMaps[2].empty();
    //
    // // Fill the interaction ntuples (ntuples 3, 4, 5 per sensor)
    // if (realisticCoincidence) {
    //     for (G4int sensorIdx = 0; sensorIdx < 3; sensorIdx++) {
    //         G4int ntupleIdx = 3 + sensorIdx;
    //         for (const auto &[coords, etot] : pixelMaps[sensorIdx]) {
    //             if (etot <= kEnergyThreshold) continue;
    //             analysisManager->FillNtupleIColumn(ntupleIdx, 0, evtID);
    //             analysisManager->FillNtupleIColumn(ntupleIdx, 1, coords.first);  // line
    //             analysisManager->FillNtupleIColumn(ntupleIdx, 2, coords.second); // col
    //             analysisManager->FillNtupleDColumn(ntupleIdx, 3, etot);
    //             analysisManager->AddNtupleRow(ntupleIdx);
    //         }
    //     }
    // }
}