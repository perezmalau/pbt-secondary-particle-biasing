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
/// \file HexitecHit.hh
/// \brief Definition of the HexitecHit class

#ifndef HexitecHit_h
#define HexitecHit_h 1

#include <utility>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class HexitecHit : public G4VHit
{
public:
    HexitecHit();
    HexitecHit(const HexitecHit &right);
    virtual ~HexitecHit();

    // operators
    const HexitecHit &operator=(const HexitecHit &);
    G4int operator==(const HexitecHit &) const;

    inline void *operator new(size_t);
    inline void  operator delete(void *aHit);

    // methods to handle data
    void Add(G4int copyNo, G4int trackID, const G4ThreeVector& intpos, G4double ediff,
        G4String particle, G4String process);

    // get methods
    G4int GetCopyNo() const { return fCopyNo; }
    G4int GetTrackID() const { return fTrackID; }
    G4String GetParticleName() const { return fParticleName; }
    G4String GetProcessName() const { return fProcessName; }
    G4ThreeVector GetInteractionPos() const { return fInteractionPos; }
    G4double GetEnergyDiff() const {return fEdiff; }

private:
    G4int           fCopyNo;
    G4int           fTrackID;
    G4String        fParticleName;
    G4String        fProcessName;
    G4ThreeVector   fInteractionPos;
    G4double        fEdiff;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using HexitecHitsCollection = G4THitsCollection<HexitecHit>;

extern G4ThreadLocal G4Allocator<HexitecHit>* HexitecHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void *HexitecHit::operator new(size_t)
{
    if (!HexitecHitAllocator)
        HexitecHitAllocator = new G4Allocator<HexitecHit>;
    void *hit;
    hit = (void *) HexitecHitAllocator->MallocSingle();
    return hit;
}

inline void HexitecHit::operator delete(void *hit)
{
    if (!HexitecHitAllocator)
        HexitecHitAllocator = new G4Allocator<HexitecHit>;
    HexitecHitAllocator->FreeSingle((HexitecHit *) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void HexitecHit::Add(G4int copyNo, G4int trackID, const G4ThreeVector& intpos, G4double ediff,
    G4String particle, G4String process)
{
    fCopyNo = copyNo;
    fTrackID = trackID;
    fParticleName = std::move(particle);
    fProcessName = std::move(process);
    fInteractionPos = intpos;
    fEdiff = ediff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
