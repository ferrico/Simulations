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
/// \file FTMDetectorMessenger.cc
/// \brief Implementation of the FTMDetectorMessenger class
//
//
// $Id: FTMDetectorMessenger.cc 98242 2016-07-04 16:57:39Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "FTMDetectorMessenger.hh"

#include "FTMDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMDetectorMessenger::FTMDetectorMessenger( FTMDetectorConstruction* Det)
: G4UImessenger(),
 fDetector(Det),
 fFTMDir(0),
 fDetDir(0),
 fr4ThicknessCmd(0)
{ 
  fFTMDir = new G4UIdirectory("/FTM/");
  fFTMDir->SetGuidance("UI commands of this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/FTM/detector/",broadcast);
  fDetDir->SetGuidance("detector control");
    
  fr4ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/FTM/detector/fr4Thickness",this);
  fr4ThicknessCmd->SetGuidance("Set Thickness of the FR4");
  fr4ThicknessCmd->SetParameterName("Size",false);
  fr4ThicknessCmd->SetRange("Size>=0.");
  fr4ThicknessCmd->SetUnitCategory("Length");
  fr4ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMDetectorMessenger::~FTMDetectorMessenger()
{
  delete fr4ThicknessCmd;
  delete fDetDir;
  delete fFTMDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  
  if( command == fr4ThicknessCmd )
   { fDetector->SetFr4Thickness(fr4ThicknessCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
