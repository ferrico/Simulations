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
// $Id: B1EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "FTMEventAction.hh"
#include "FTMRunAction.hh"
#include "FTMAnalysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMEventAction::FTMEventAction(FTMRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fNprimaries(0.),
  fNsecondaryParticlesGen(0.),
  fNsecondaryEle(0.),
  fNsecondaryEleBottom(0.),
  fNprimariesTop(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMEventAction::~FTMEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMEventAction::BeginOfEventAction(const G4Event* event)
{
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if (eventID%printModulo == 0)
    //G4cout << "--- Event ID: " << eventID << G4endl;
    
  // initialisation per event
  fEdep = 0.;
  fNprimaries = 0.;
  fNsecondaryParticlesGen = 0.;
  fNsecondaryEle = 0.;
  fNsecondaryEleBottom = 0.;
  fNprimariesTop = 0.;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMEventAction::EndOfEventAction(const G4Event* event)
{
  auto eventID = event->GetEventID();
  // accumulate statistics in run action
  fRunAction->AddEvents(fEdep);
    
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
    
  // fill histograms
  analysisManager->FillH1(0, 1); //nEvents
  analysisManager->FillH1(1, fNprimaries);
  analysisManager->FillH1(2, fNsecondaryParticlesGen);
  analysisManager->FillH1(3, fNsecondaryEle);
  analysisManager->FillH1(6, fEdep);
  analysisManager->FillH1(9, fNsecondaryEleBottom);
  analysisManager->FillH1(12, fNprimariesTop);
    
  // Print per event (modulo n)
  //
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
      //G4cout << "---> End of event: " << eventID << G4endl;
      /*G4cout
      << "   Absorber: total energy: " << std::setw(7)
      << G4BestUnit(fEnergyAbs,"Energy")
      << "       total track length: " << std::setw(7)
      << G4BestUnit(fTrackLAbs,"Length")
      << G4endl
      << "        Gap: total energy: " << std::setw(7)
      << G4BestUnit(fEnergyGap,"Energy")
      << "       total track length: " << std::setw(7)
      << G4BestUnit(fTrackLGap,"Length")
      << G4endl;*/
  }
  /*G4cout << "--- Event ID: " << eventID << " # primaries " << fNprimaries << G4endl;
  G4cout << "--- Event ID: " << eventID << " # secondary electrons generated from the initial gamma " << fNsecondaryParticlesGen << G4endl;
  G4cout << "--- Event ID: " << eventID << " # secondary electrons outcoming from FR4 " << fNsecondaryEle << G4endl;
  G4cout << "--- Event ID: " << eventID << " tot. energy in FR4 " << G4BestUnit(fEdep,"Energy") << G4endl;*/
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
