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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "FTMSteppingAction.hh"
#include "FTMEventAction.hh"
#include "FTMDetectorConstruction.hh"
#include "FTMAnalysis.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMSteppingAction::FTMSteppingAction(FTMDetectorConstruction* detector, FTMEventAction* eventAction)
: G4UserSteppingAction(),
  fDetector(detector),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMSteppingAction::~FTMSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMSteppingAction::UserSteppingAction(const G4Step* step)
{
  auto stepID = step->GetTrack()->GetCurrentStepNumber();
  //G4cout << "--- Step ID: " << stepID << G4endl;
    
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
    
  // get volume of the current step
  //G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;
  ////////////////////////////////////////////
  G4int trackID = step->GetTrack()->GetTrackID() ;			// track ID
  G4int parentID = step->GetTrack()->GetParentID() ;			// parent ID
  G4int particlePDG = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  G4double particleCharge= step->GetTrack()->GetParticleDefinition()->GetPDGCharge();
  G4double particleEk = step->GetTrack()->GetKineticEnergy();
  G4double edepStep = step->GetTotalEnergyDeposit();
  if (parentID == 0 && trackID == 1){
        EleTrkID = 0;
        edepStepEle = 0;
  }
    
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();
    
  ///////// outside FR4 (Top Layer) ///////////////////////
  if( volume == fDetector->GetTop()){
      // Primaries in Top//
    if(parentID == 0 && postStepPoint->GetStepStatus() == fGeomBoundary){
        //G4cout << "--- Step ID: " << stepID << " Primary gamma on Top " << G4endl;
        if (stepID <3 ){ // tentativo di evitare di ricontare i gamma che vengono deviati all'indietro da compton
        fEventAction->AddPrimariesTop(1);
        analysisManager->FillH1(13, 1);
        }
    }
  }
    ///////// outside FR4 (Bottom Layer) ///////////////////////
    
  ///////// inside FR4 ///////////////////////
  //if( volume != fDetector->GetFr4()) return;
  if( volume == fDetector->GetFr4()){
    //G4cout << "--- Step ID: " << stepID << " in FR4 (scoring volume) " << G4endl;
    /*G4cout << "--- Step ID: " << stepID << " track ID " << trackID << G4endl;
    G4cout << "--- Step ID: " << stepID << " parent ID " << parentID << G4endl;
    G4cout << "--- Step ID: " << stepID << " particle PDG " << particlePDG << G4endl;
    G4cout << "--- Step ID: " << stepID << " particle Charge " << particleCharge << G4endl;
    G4cout << "--- Step ID: " << stepID << " particle Ek " << G4BestUnit(particleEk,"Energy") << G4endl;
    G4cout << "--- Step ID: " << stepID << " E deposited " << G4BestUnit(edepStep,"Energy") << G4endl;*/
  

    // Use the GetStepStatus() method of G4StepPoint to get the status of the
    // current step (contained in post­step point) or the previous step
    // (contained in pre­step point):
    if(preStepPoint->GetStepStatus() == fGeomBoundary)
      //G4cout << "Step starts on geometry boundary" << G4endl;
    if(postStepPoint->GetStepStatus() == fGeomBoundary)
      //G4cout << "Step ends on geometry boundary" << G4endl;
    
    // Primitive gamma on FR4
    if(parentID == 0 && preStepPoint->GetStepStatus() == fGeomBoundary){
      //G4cout << "--- Step ID: " << stepID << " Primary gamma on FR4 " << G4endl;
      fEventAction->AddPrimaries(1);
      analysisManager->FillH1(4, 1);
    }
      
    // generated secondary electrons
    //G4int genSecondaries = step->GetSecondary()->size();
    G4int genSecondaries = step->GetSecondaryInCurrentStep()->size();
    //G4cout << "--- Step ID: " << stepID << " Secondary electron generated " << genSecondaries << G4endl;
    if(parentID == 0 ){
      const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
      fEventAction->AddSecondaryParticlesGen(genSecondaries);
      for (size_t lp=0; lp<(*secondary).size(); lp++) {
          if ((*secondary)[lp]->GetDefinition()->GetPDGEncoding() == 11){
              analysisManager->FillH1(8, (*secondary)[lp]->GetKineticEnergy());// gen ele kin. energy
              EleGenEk = (*secondary)[lp]->GetKineticEnergy();
              //G4cout << "--- Step ID: " << stepID << " Secondary electron generated with Ek = " << G4BestUnit((*secondary)[lp]->GetKineticEnergy(),"Energy")<< G4endl;
          } // if secondary generated is an electron
      }// loop on secondaries gen
    } // if primary
    
    // Collect energy dep. from secondary electrons in FR4
    if(parentID != 0 && particlePDG == 11 ) {
        if (EleTrkID != trackID) {
            edepStepEle = 0;}
        EleTrkID = trackID;
        edepStepEle += edepStep;
        //G4cout << "--- Step ID: " << stepID << " Secondary electron Edep in FR4 " << G4BestUnit(edepStepEle,"Energy") << G4endl;
        
        const std::vector<const G4Track*>* secondary_e = step->GetSecondaryInCurrentStep();
        for (size_t lpe=0; lpe<(*secondary_e).size(); lpe++) {
            //if ((*secondary_e)[lpe]->GetDefinition()->GetPDGEncoding() == 11){
                EleGenEk_e = (*secondary_e)[lpe]->GetKineticEnergy();
                //G4cout << "---! Step ID: " << stepID << " Secondary particle generated by e with Ek = " << G4BestUnit((*secondary_e)[lpe]->GetKineticEnergy(),"Energy")<< G4endl;
                edepStepEle += EleGenEk_e;
            //} // if secondary generated is an electron
        }// loop on secondaries gen
        
        //if (postStepPoint->GetStepStatus() == fGeomBoundary){
        //G4cout << "--- Step ID: " << stepID << " Secondary electron Edep in FR4 " << G4BestUnit(step->GetTotalEnergyDeposit(),"Energy") << G4endl;
        //G4cout << "--- Step ID: " << stepID << " Secondary electron DeltaE in FR4 " << G4BestUnit(step->GetDeltaEnergy(),"Energy") << G4endl;
        //G4cout << "--- Step ID: " << stepID << " Secondary electron Cumulative Edep in FR4 " << G4BestUnit(edepStepEle,"Energy") << G4endl;
        //G4cout << "--- Step ID: " << stepID << " Secondary electron Ek gen " << G4BestUnit(EleGenEk,"Energy") << G4endl;
        //G4cout << "--- Step ID: " << stepID << " Secondary electron Ek gen from Ek + Edep " << G4BestUnit(particleEk + edepStepEle,"Energy") << G4endl;
        //}
    }
      
    // Secondaries electrons exiting from FR4//
    if(parentID != 0 && particlePDG == 11 && postStepPoint->GetStepStatus() == fGeomBoundary) {
      //G4cout << "--- Step ID: " << stepID << " Secondary electron of Ek " << G4BestUnit(particleEk,"Energy") << G4endl;
      fEventAction->AddSecondaryEle(1);// n. ele per event
      analysisManager->FillH1(5, 1);// counting ele
      analysisManager->FillH1(7, particleEk);// ele kin. energy
    }
  
    /// collect energy deposited in the scoring volume
    fEventAction->AddEdep(edepStep);
      
  } // if inside FR4
  ///////// inside FR4 ///////////////////////
    
  ///////// outside FR4 (Bottom Layer) ///////////////////////
  if( volume == fDetector->GetBottom()){
      
    // Primitive gamma on FR4
    if(parentID == 0 && preStepPoint->GetStepStatus() == fGeomBoundary){
        //G4cout << "--- Step ID: " << stepID << " Primary gamma in Bottom of Ek " << G4BestUnit(particleEk,"Energy") << G4endl;
        analysisManager->FillH1(14, 1);// counting init gamma in bottom
        analysisManager->FillH1(15, particleEk);// ele kin. energy
    }
      
    // Secondaries electrons exiting from FR4//
    if(parentID != 0 && particlePDG == 11 && preStepPoint->GetStepStatus() == fGeomBoundary) {
        //if (step->GetStepLength()==0.) G4cout << "--- Step ID: " << stepID << " SOSPETTO " << G4endl;
        if(step->GetStepLength()>0.){// fill only with final exiting ele
        //G4cout << "--- Step ID: " << stepID << " Secondary electron in Bot. of Ek " << G4BestUnit(particleEk,"Energy") << G4endl;
        fEventAction->AddSecondaryEleBottom(1);// n. ele per event
        //G4cout << "--- Step ID: " << stepID << " Secondary electron Ek gen from Ek + Edep " << G4BestUnit(particleEk + edepStepEle,"Energy") << G4endl;
        //if((particleEk + edepStepEle)>340*keV && (particleEk + edepStepEle)<400*keV)G4cout << "--- Step ID: " << stepID << " !!! Secondary electron Ek gen from Ek + Edep " << G4BestUnit(particleEk + edepStepEle,"Energy") << G4endl;
        analysisManager->FillH1(10, 1);// counting ele
        analysisManager->FillH1(11, particleEk);// ele kin. energy
        analysisManager->FillH1(16, particleEk + edepStepEle);// ele kin. energy gen (only for e- exiting)
        }
    }
  }
  ///////// outside FR4 (Bottom Layer) ///////////////////////
  
    
  //////////////////////////////////////
  // collect energy deposited in this step
  //G4double stepl = 0.;
  //if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
  //stepl = step->GetStepLength();
  ///////////////////////////////////////////////////////

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

