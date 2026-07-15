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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include "G4SystemOfUnits.hh"


SteppingAction::SteppingAction()
        : G4UserSteppingAction()
{}


SteppingAction::~SteppingAction()
{}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4String vol = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
    auto analysisManager = G4AnalysisManager::Instance();
    if (vol == "test.phys") {
        // Get total energy deposit
        G4double edep = step->GetTotalEnergyDeposit()/MeV;
        if (edep > 0.) {
            G4ThreeVector pos = 0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());
            analysisManager->FillH1(0, pos.x(), edep);
        }

        // create a list of all secondaries created in this step
        const std::vector<const G4Track *> *secondaries = step->GetSecondaryInCurrentStep();
        // loop through and send them to the analysis manager
        for (auto track: *secondaries) {
            G4double ekin = track->GetKineticEnergy();
            G4String particle = track->GetParticleDefinition()->GetParticleName();
            if (particle == "gamma") {
                G4double posX = track->GetPosition().getX();
                G4double posY = track->GetPosition().getY();
                if (ekin > 2 and ekin < 7) {
                    analysisManager->FillH2(0, posX, posY, 1);
                }
                if (ekin > 2 and ekin < 2.4) {
                    analysisManager->FillH2(1, posX, posY, 1);
                }
                if (ekin > 4.2 and ekin < 4.6) {
                    analysisManager->FillH2(2, posX, posY, 1);
                }
                if (ekin > 5 and ekin < 5.4) {
                    analysisManager->FillH2(3, posX, posY, 1);
                }
                if (ekin > 5.9 and ekin < 6.3) {
                    analysisManager->FillH2(4, posX, posY, 1);
                }
            }
        }
    }
}