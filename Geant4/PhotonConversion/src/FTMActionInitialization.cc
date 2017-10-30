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
// $Id: B1ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B1ActionInitialization.cc
/// \brief Implementation of the B1ActionInitialization class

#include "FTMActionInitialization.hh"
#include "FTMPrimaryGeneratorAction.hh"
#include "FTMRunAction.hh"
#include "FTMEventAction.hh"
#include "FTMSteppingAction.hh"
#include "SteppingVerbose.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMActionInitialization::FTMActionInitialization(FTMDetectorConstruction* detector)
 : G4VUserActionInitialization(),
   fDetector(detector)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMActionInitialization::~FTMActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMActionInitialization::BuildForMaster() const
{
  // Histo manager
  FTMRunAction* runAction = new FTMRunAction();
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMActionInitialization::Build() const
{    
  SetUserAction(new FTMPrimaryGeneratorAction(fDetector));

  FTMRunAction* runAction = new FTMRunAction();
  SetUserAction(runAction);
  
  FTMEventAction* eventAction = new FTMEventAction(runAction);
  SetUserAction(eventAction);
  
  SetUserAction(new FTMSteppingAction(fDetector, eventAction));
  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose* FTMActionInitialization::InitializeSteppingVerbose() const
{
    return new SteppingVerbose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
