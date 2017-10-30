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
// $Id: B1EventAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef FTMEventAction_h
#define FTMEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class FTMRunAction;

/// Event action class
///

class FTMEventAction : public G4UserEventAction
{
  public:
    FTMEventAction(FTMRunAction* runAction);
    virtual ~FTMEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    //in stepping action?
    void AddEdep(G4double edep) { fEdep += edep; };
    void AddPrimaries(G4int primaries) { fNprimaries += primaries; };
    void AddSecondaryParticlesGen(G4int secondaryPGen) { fNsecondaryParticlesGen += secondaryPGen; };
    void AddSecondaryEle(G4int secondaryEle) { fNsecondaryEle += secondaryEle; };
    void AddSecondaryEleBottom(G4int secondaryEleBottom) { fNsecondaryEleBottom += secondaryEleBottom; };
    void AddPrimariesTop(G4int primariesTop) { fNprimariesTop += primariesTop; };

    //void AddAbs(G4double de, G4double dl);


  private:
    FTMRunAction* fRunAction;
    G4double  fEdep;
    G4int  fNprimaries;
    G4int  fNsecondaryParticlesGen;
    G4int  fNsecondaryEle;
    G4int  fNsecondaryEleBottom;
    G4int  fNprimariesTop;

};

// Functions Definition
/*inline void FTMEventAction::AddAbs(G4double de, G4double dl) {
    fEnergyAbs += de;
    fTrackLAbs += dl;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
