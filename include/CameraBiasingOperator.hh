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

// CameraBiasingOperator.hh
#ifndef CAMERABIASINGOPERATOR_HH
#define CAMERABIASINGOPERATOR_HH
#include "G4VBiasingOperator.hh"
#include "G4ThreeVector.hh"
class CameraDirectionBiasingOperation;

class CameraBiasingOperator : public G4VBiasingOperator
{
public:
    CameraBiasingOperator(G4ThreeVector  cameraCenter, G4double cameraHalfXY);

    void StartRun() override;

    G4VBiasingOperation* ProposeOccurenceBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override
    { return nullptr; }

    G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override
    { return nullptr; }

    G4VBiasingOperation* ProposeFinalStateBiasingOperation(
        const G4Track* track,
        const G4BiasingProcessInterface* callingProcess) override;

private:
    CameraDirectionBiasingOperation* fBiasingOp;
    G4ThreeVector fCameraCenter;
    G4double      fCameraHalfXY;
};
#endif