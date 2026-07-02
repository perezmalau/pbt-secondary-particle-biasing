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
/// \file src/RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4Timer.hh"
#include "G4AnalysisManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction()
{
    timer = new G4Timer();
  	// Create analysis manager
  	// The choice of analysis technology is done via selection of a namespace
  	// in Analysis.hh
  	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetDefaultFileType("root");
  	analysisManager->SetFileName("HexComptCam");
  	G4cout << "Using " << analysisManager->GetType() << G4endl;
  	analysisManager->SetVerboseLevel(1);
  	analysisManager->SetNtupleMerging(true); //only with root


  	// create Ntuple for ideal gamma info in scatterer1, id = 0
	analysisManager->CreateNtuple("G4S1TrueGammaInfo", "Sensor1 ideal hits");
	analysisManager->CreateNtupleIColumn("EventID");
	analysisManager->CreateNtupleIColumn("TrackID");
	analysisManager->CreateNtupleDColumn("Xpos");
	analysisManager->CreateNtupleDColumn("Ypos");
    analysisManager->CreateNtupleDColumn("Zpos");
	analysisManager->CreateNtupleDColumn("Etot");
	analysisManager->CreateNtupleSColumn("Process");
	// analysisManager->CreateNtupleDColumn("TrueAngle");
	// analysisManager->CreateNtupleDColumn("TrueEnergy");
	analysisManager->FinishNtuple();

  	// create Ntuple for ideal gamma info in scatterer2, id = 1
	analysisManager->CreateNtuple("G4S2TrueGammaInfo", "Sensor2 ideal hits");
	analysisManager->CreateNtupleIColumn("EventID");
	analysisManager->CreateNtupleIColumn("TrackID");
	analysisManager->CreateNtupleDColumn("Xpos");
	analysisManager->CreateNtupleDColumn("Ypos");
	analysisManager->CreateNtupleDColumn("Zpos");
	analysisManager->CreateNtupleDColumn("Etot");
	analysisManager->CreateNtupleSColumn("Process");
	// analysisManager->CreateNtupleDColumn("TrueAngle");
	// analysisManager->CreateNtupleDColumn("TrueEnergy");
	analysisManager->FinishNtuple();

    // create Ntuple for ideal gamma info in scatterer3, id = 2
    analysisManager->CreateNtuple("G4S3TrueGammaInfo", "Sensor3 ideal hits");
	analysisManager->CreateNtupleIColumn("EventID");
	analysisManager->CreateNtupleIColumn("TrackID");
	analysisManager->CreateNtupleDColumn("Xpos");
	analysisManager->CreateNtupleDColumn("Ypos");
	analysisManager->CreateNtupleDColumn("Zpos");
	analysisManager->CreateNtupleDColumn("Etot");
	analysisManager->CreateNtupleSColumn("Process");
	// analysisManager->CreateNtupleDColumn("TrueAngle");
	// analysisManager->CreateNtupleDColumn("TrueEnergy");
    analysisManager->FinishNtuple();

    // create Ntuple for measured sensor1 hits, id = 3
    analysisManager->CreateNtuple("G4Sensor1Hits", "Sensor1 hits");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("Line");
    analysisManager->CreateNtupleIColumn("Col");
    analysisManager->CreateNtupleDColumn("Edep");
    analysisManager->FinishNtuple();

    // create Ntuple for measured sensor2 hits, id = 4
    analysisManager->CreateNtuple("G4Sensor2Hits", "Sensor2 hits");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("Line");
    analysisManager->CreateNtupleIColumn("Col");
    analysisManager->CreateNtupleDColumn("Edep");
    analysisManager->FinishNtuple();

    // create Ntuple for measured sensor3 hits, id = 5
    analysisManager->CreateNtuple("G4Sensor3Hits", "Sensor3 hits");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("Line");
    analysisManager->CreateNtupleIColumn("Col");
    analysisManager->CreateNtupleDColumn("Edep");
    analysisManager->FinishNtuple();

	analysisManager->SetH1Activation(true);
	analysisManager->SetH3Activation(true);
	analysisManager->CreateH3("Dose3D", "3D Dose Distribution",
	70, -7.*cm, 7.*cm,
	100, -10.*cm, 10.*cm,
	85, -8.5*cm, 8.5*cm);

	analysisManager->CreateH3("Gamma478", "3D Distribution of 478 keV Gammas",
		70, -7.*cm, 7.*cm,
		100, -10.*cm, 10.*cm,
		85, -8.5*cm, 8.5*cm);

	analysisManager->CreateH3("Gamma2200", "3D Distribution of 2.2 MeV Gammas",
		70, -7.*cm, 7.*cm,
		100, -10.*cm, 10.*cm,
		85, -8.5*cm, 8.5*cm);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
  delete timer;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	if (isMaster) {
        timer->Start();
    }
    G4AnalysisManager *analysisManager =  G4AnalysisManager::Instance();
	analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run * /*aRun*/)
{
  	auto analysisManager = G4AnalysisManager::Instance();
  	if (isMaster) {
		timer->Stop();
	  	G4cout << "Simulation time: " << timer->GetRealElapsed() << " seconds." << G4endl;
	}

  	analysisManager->Write();
  	analysisManager->CloseFile();
}
