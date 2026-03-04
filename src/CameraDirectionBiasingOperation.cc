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
// CameraDirectionBiasingOperation.cc
#include "CameraDirectionBiasingOperation.hh"
#include "G4ParticleChange.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <G4BiasingProcessInterface.hh>
#include <utility>

CameraDirectionBiasingOperation::CameraDirectionBiasingOperation(
    const G4String& name, G4ThreeVector  cameraCenter, G4double cameraHalfXY)
  : G4VBiasingOperation(name),
    fCameraCenter(std::move(cameraCenter)),
    fCameraHalfXY(cameraHalfXY)
{}

G4VParticleChange*
CameraDirectionBiasingOperation::ApplyFinalStateBiasing(
    const G4BiasingProcessInterface* callingProcess,
    const G4Track*                   track,
    const G4Step*                    step,
          G4bool&                    forceBiasedFinalState)
{
    // Let the underlying process generate its natural final state
    G4VParticleChange* processFinalState =
        callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);

    G4int nSecondaries = processFinalState->GetNumberOfSecondaries();
    for (G4int i = 0; i < nSecondaries; ++i)
    {
        G4Track* secondary = processFinalState->GetSecondary(i);

        if (secondary->GetDefinition() != G4Gamma::Definition()) continue;

        // Direction from production point toward camera center
        const G4ThreeVector& prodPos  = track->GetPosition();
        // Four corners of the camera front face
        G4double xc = fCameraCenter.x();
        G4double yc = fCameraCenter.y();
        G4double h  = fCameraHalfXY;
        G4double z  = 0; // Camera is positioned such that the front face is at z=0
        std::array<G4ThreeVector, 4> corners = {{
            G4ThreeVector(xc + h, yc + h, z),
            G4ThreeVector(xc + h, yc - h, z),
            G4ThreeVector(xc - h, yc + h, z),
            G4ThreeVector(xc - h, yc - h, z)
        }};

        // Direction from production point to camera centre
        G4ThreeVector toCenter = (G4ThreeVector(xc, yc, z) - prodPos).unit();

        // Half-angle = largest angle between toCenter and any corner direction
        G4double halfAngle = 0.0;
        for (const auto& corner : corners) {
            G4ThreeVector toCorner = (corner - prodPos).unit();
            G4double angle = std::acos(std::min(1.0, toCenter.dot(toCorner)));
            halfAngle = std::max(halfAngle, angle);
        }

        G4double cameraSolidAngle = 2.0 * CLHEP::pi * (1.0 - std::cos(halfAngle));

        G4ThreeVector biasedDir   = SampleDirectionInCone(toCenter, halfAngle);
        G4double weightCorrection = cameraSolidAngle / (4.0 * CLHEP::pi);

        secondary->SetMomentumDirection(biasedDir);
        secondary->SetWeight(secondary->GetWeight() * weightCorrection);
    }

    return processFinalState;
}

G4ThreeVector
CameraDirectionBiasingOperation::SampleDirectionInCone(
    const G4ThreeVector& axis, G4double halfAngle)
{
    // Sample uniformly in solid angle within cone
    // cosθ uniform in [cos(halfAngle), 1]
    G4double cosMin   = std::cos(halfAngle);
    G4double cosTheta = cosMin + (1.0 - cosMin) * G4UniformRand();
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    G4double phi      = 2.0 * CLHEP::pi * G4UniformRand();

    // Build in local frame around z-axis, then rotate to axis
    G4ThreeVector localDir(sinTheta * std::cos(phi),
                           sinTheta * std::sin(phi),
                           cosTheta);

    // Rotate local z → axis
    G4ThreeVector zAxis(0., 0., 1.);
    G4ThreeVector rotAxis = zAxis.cross(axis);
    G4double sinA = rotAxis.mag();
    G4double cosA = zAxis.dot(axis);

    if (sinA < 1e-10) // axis is already ±z
        return (cosA > 0) ? localDir : G4ThreeVector(-localDir.x(),
                                                      -localDir.y(),
                                                      -localDir.z());

    rotAxis /= sinA;
    localDir.rotate(std::atan2(sinA, cosA), rotAxis);
    return localDir;
}