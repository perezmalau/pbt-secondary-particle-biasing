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

    // Define pixel maps - the realistic detector information
    std::map<std::pair<int,int>, G4double> pixelMaps[3];

    // Define a gamma record with the true gamma information (set of sensor hit and true pos and process info)
    struct GammaRecord {
        // The actual set of hits in the sensors
        std::set<G4int> sensorsHit;
        // The position, process and weight are recorded as true information
        std::map<G4int, std::vector<std::tuple<G4ThreeVector, G4double, G4double, G4double >>> trueInfo;
    };

    std::map<G4int, GammaRecord> gammaRecords; // key: trackID of the gamma

    // Loop through all hits in the event, record the relevant information
    for (G4int i = 0; i < nSensor; i++) {
        auto hit      = (*sensorHC)[i];
        G4int copyNo  = hit->GetCopyNo();
        // Accumulate pixel energy
        G4int line = hit->GetLineNumber();
        G4int col  = hit->GetColNumber();
        pixelMaps[copyNo][{line, col}] += hit->GetEdep();

        // Only fill true info for primary gammas
        G4int    parentID = hit->GetParentID();
        G4String particle = hit->GetParticleName();
        if (particle == "gamma" and parentID == 1) {
            G4int gammaTrackID = hit->GetTrackID();
            auto &record = gammaRecords[gammaTrackID];
            record.sensorsHit.insert(copyNo);
            record.trueInfo[copyNo].push_back({ hit->GetInteractionPos(), hit->GetInitialEnergy(), hit->GetEnergyDiff(), hit->GetWeight()});
        }
    }

    auto analysisManager = G4AnalysisManager::Instance();

    // Fill true information based on trackID coincidences
    for (auto &[gammaTrackID, record] : gammaRecords) {
        // Require hits in all 3 sensors
        if (record.sensorsHit.size() < 3) continue;
        // Fill "true" interaction ntuples (ntuples 0, 1, 2 per sensor)
        for (G4int sensorIdx = 0; sensorIdx < 3; sensorIdx++) {
            auto it = record.trueInfo.find(sensorIdx);
            if (it == record.trueInfo.end()) continue;
            for (const auto &[pos, ene, ediff, weight] : record.trueInfo[sensorIdx]) {
                analysisManager->FillNtupleIColumn(sensorIdx, 0, evtID);
                analysisManager->FillNtupleIColumn(sensorIdx, 1, gammaTrackID);
                analysisManager->FillNtupleDColumn(sensorIdx, 2, pos.x()/CLHEP::cm);
                analysisManager->FillNtupleDColumn(sensorIdx, 3, pos.y()/CLHEP::cm);
                analysisManager->FillNtupleDColumn(sensorIdx, 4, pos.z()/CLHEP::cm);
                analysisManager->FillNtupleDColumn(sensorIdx, 5, ene);
                analysisManager->FillNtupleDColumn(sensorIdx, 6, ediff);
                //analysisManager->FillNtupleDColumn(sensorIdx, 7, weight);
                analysisManager->AddNtupleRow(sensorIdx);
            }
        }
    }

    // Fill realistic detector information based on event coincidences and a detector threshold at 3 keV
    constexpr G4double kEnergyThreshold = 0.003; // MeV

    // Realistic coincidence: all 3 pixel maps are non-empty
    bool realisticCoincidence = !pixelMaps[0].empty() &&
                                !pixelMaps[1].empty() &&
                                !pixelMaps[2].empty();

    // Fill the interaction ntuples (ntuples 3, 4, 5 per sensor)
    if (realisticCoincidence) {
        for (G4int sensorIdx = 0; sensorIdx < 3; sensorIdx++) {
            G4int ntupleIdx = 3 + sensorIdx;
            for (const auto &[coords, etot] : pixelMaps[sensorIdx]) {
                if (etot <= kEnergyThreshold) continue;
                analysisManager->FillNtupleIColumn(ntupleIdx, 0, evtID);
                analysisManager->FillNtupleIColumn(ntupleIdx, 1, coords.first);  // line
                analysisManager->FillNtupleIColumn(ntupleIdx, 2, coords.second); // col
                analysisManager->FillNtupleDColumn(ntupleIdx, 3, etot);
                analysisManager->AddNtupleRow(ntupleIdx);
            }
        }
    }
}

