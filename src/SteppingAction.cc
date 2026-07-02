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
    G4double phantCamDistance = 20*cm;
    G4String vol = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
    auto analysisManager = G4AnalysisManager::Instance();
    if (vol == "Head" || vol == "Skull" || vol == "Brain" || vol == "Target") {
        // Get total energy deposit
        G4double edep = step->GetTotalEnergyDeposit()/MeV;
        if (edep > 0.) {
            G4ThreeVector pos = 0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());
            analysisManager->FillH3(0, pos.x(), pos.y(), pos.z() + phantCamDistance, edep);
        }

        // create a list of all secondaries created in this step
        const std::vector<const G4Track *> *secondaries = step->GetSecondaryInCurrentStep();
        // loop through and send them to the analysis manager
        for (auto track : *secondaries) {
            if (track->GetParticleDefinition()->GetParticleName() != "gamma") continue;

            G4double ekin = track->GetKineticEnergy();
            G4ThreeVector pos = track->GetPosition();
            G4double x = pos.x();
            G4double y = pos.y();
            G4double z = pos.z() + phantCamDistance;

            if (std::abs(ekin - 478.*keV) < 5.*keV) {
                analysisManager->FillH3(1, x, y, z, 1);
            } else if (std::abs(ekin - 2223.*keV) < 5.*keV) {
                analysisManager->FillH3(2, x, y, z, 1);
            }
        }
    }
}