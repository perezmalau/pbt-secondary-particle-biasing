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

// CameraDirectionBiasingOperation.hh
#ifndef CAMERADIRECTIONBIASINGOPERATION_HH
#define CAMERADIRECTIONBIASINGOPERATION_HH
#include "G4VBiasingOperation.hh"
#include "G4ThreeVector.hh"

class CameraDirectionBiasingOperation : public G4VBiasingOperation
{
public:
    CameraDirectionBiasingOperation(const G4String& name,
                                     G4ThreeVector  cameraCenter,
                                     G4double cameraHalfXY);

    // We only use final state biasing
    const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw(
        const G4BiasingProcessInterface*, G4ForceCondition&) override
    { return nullptr; }

    G4VParticleChange* ApplyFinalStateBiasing(
    const G4BiasingProcessInterface* callingProcess,
    const G4Track*                   track,
    const G4Step*                    step,
          G4bool&                    forceBiasedFinalState) override;

    G4double DistanceToApplyOperation(const G4Track*, G4double,
        G4ForceCondition*) override { return DBL_MAX; }

    G4VParticleChange* GenerateBiasingFinalState(
        const G4Track*, const G4Step*) override { return nullptr; }

private:
    static G4ThreeVector SampleDirectionInCone(const G4ThreeVector& axis,
                                         G4double halfAngle);

    G4ThreeVector fCameraCenter;
    G4double fCameraHalfXY;
};
#endif // CAMERADIRECTIONBIASINGOPERATION_HH
