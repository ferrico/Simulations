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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "FTMPrimaryGeneratorAction.hh"
#include "FTMDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMPrimaryGeneratorAction::FTMPrimaryGeneratorAction(FTMDetectorConstruction* DC)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fDetector(DC)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(511.*keV);
    
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Detector volume
  //G4double detectorSizeXY = fDetector->GetDetectorSizeXY();
  G4double detectorSizeZ = fDetector->GetDetectorSizeZ();
  
  //G4double size = 0.8; // for the 80% variation of the X Y position
  //G4double x0 = size * detectorSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = size * detectorSizeXY * (G4UniformRand()-0.5);
  //G4double size = 1;
  //G4double x0 = 0.0;
  //G4double y0 = 0.0;
  //G4double topThickness = 3.*mm;
  //G4double z0 = -detectorSizeZ/2;// + topThickness ;
    
  G4cout
  << "------------------------------------------------------------" << G4endl
  << "------------------------------------------------------------" << G4endl
  << "------------------------------------------------------------" << G4endl
  //<< " z0 of the gun "
  //<< G4BestUnit(z0, "Length")
  << " detectorSizeZ "
  << G4BestUnit(detectorSizeZ, "Length")
  << G4endl;
    
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*mm,0.*mm,-detectorSizeZ*0.5));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMPrimaryGeneratorAction::~FTMPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

