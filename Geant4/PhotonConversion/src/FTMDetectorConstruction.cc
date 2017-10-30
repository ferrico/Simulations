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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "FTMDetectorConstruction.hh"
#include "FTMDetectorMessenger.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GenericMessenger.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMDetectorConstruction::FTMDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true),
  fScoringVolume(0),
  fNofLayers(-1),
  fDetectorMessenger(0),
  fFR4PV(0),
  fTopPV(0),
  fBottomPV(0)
{
    // default parameter values of the detector
    fFr4Thickness = 10.*mm;

    fDetectorSizeXY = 10.*cm;
    fDetectorThickness  = 10.*cm;
    fWorldSizeXY  = 10.*cm;
    fWorldSizeZ   = 10.*cm;
    
    // materials
    DefineMaterials();
    
    // create commands for interactive definition of the detector
    fDetectorMessenger = new FTMDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FTMDetectorConstruction::~FTMDetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void FTMDetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // variables
  G4int numel(0), natoms(0) ;
  G4double density(0.), fractionMass(0.)  ;//temperature(0.), pressure(0.)
  
  // define Elements
  G4Element* elH  = nist->FindOrBuildElement(1);
  G4Element* elC  = nist->FindOrBuildElement(6);
  G4Element* elSi = nist->FindOrBuildElement(14);
  G4Element* elO  = nist->FindOrBuildElement(8);
    
  // define Materials
  // Vacuum //
  G4Material* empty = nist->FindOrBuildMaterial("G4_Galactic");
  vacuumMaterial = empty ;
    
  //Air// ??? do I need to re-define it? why?
  G4Material *Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR") ;
  airMaterial = Air;
    
  //FR4//
  //Epoxy (for FR4 )
  //from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm //???
  density = 1.2*g/cm3;
  G4Material* Epoxy = new G4Material("Epoxy" , density, numel=2);
  Epoxy->AddElement(elH, natoms=2);
  Epoxy->AddElement(elC, natoms=2);
  //SiO2 (Quarz)
  G4Material* SiO2 =  new G4Material("SiO2",density= 2.200*g/cm3, numel=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);
  //FR4 (Glass + Epoxy)
  density = 1.86*g/cm3;
  G4Material* FR4 = new G4Material("FR4"  , density, numel=2);
  FR4->AddMaterial(Epoxy, fractionMass=0.472);
  FR4->AddMaterial(SiO2, fractionMass=0.528);
  fr4Material = FR4;
    
  //Kapton//
  G4Material *Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
  kaptonMaterial = Kapton;
  //G4_GLASS_LEAD
  G4Material *GlassLead = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_LEAD");
  glassLeadMaterial = GlassLead;
  //G4_GLASS_PLATE
  G4Material *GlassPlate = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
  glassPlateMaterial = GlassPlate;
    
  //Cu//
  G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu") ;
  cuMaterial = Cu;
  
  // Ar:CO2 (70:30) @ STP conditions
  //gases at STP conditions
  G4Material* Argon = nist->FindOrBuildMaterial("G4_Ar");
  G4Material* CarbonDioxide = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4double mixtureDensity = (Argon->GetDensity() * 70/100.0 + CarbonDioxide->GetDensity() * 30/100.0) ;
  G4Material *ArCO2 = new G4Material("Ar/CO2",mixtureDensity,2) ;
  ArCO2->AddMaterial(Argon, 0.7) ;
  ArCO2->AddMaterial(CarbonDioxide, 0.3) ;
  
  // Choice of the gas
  gasMaterial = ArCO2 ;
  //  fGasMat = fAirMat ;

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* FTMDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double topThickness = 0.5*mm;
  G4double bottomThickness = 0.5*mm;
    
  fDetectorThickness = topThickness + fFr4Thickness + bottomThickness;
  fWorldSizeXY = 2.5 * fDetectorSizeXY;
  fWorldSizeZ  = 2.5 * fDetectorThickness;
    
  // Check materials
    
  if ( ! vacuumMaterial || ! fr4Material || ! cuMaterial || ! gasMaterial || ! airMaterial ) {
        G4ExceptionDescription msg;
        msg << "Cannot retrieve materials already defined.";
        G4Exception("B4DetectorConstruction::DefineVolumes()",
                    "MyCode0001", FatalException, msg);
  }
    
  //
  // World
  //
  auto worldS
  = new G4Box("World",           // its name
                fWorldSizeXY/2, fWorldSizeXY/2, fWorldSizeZ/2); // its size
    
  auto worldLV
  = new G4LogicalVolume(
                          worldS,           // its solid
                          vacuumMaterial,   // its material
                          "World");         // its name
    
  auto worldPV
  = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        worldLV,          // its logical volume
                        "World",          // its name
                        0,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
  //
  // Detector
  //
  auto detectorS
  = new G4Box("Detector",     // its name
                fDetectorSizeXY/2, fDetectorSizeXY/2, fDetectorThickness/2); // its size
    
  auto detectorLV
  = new G4LogicalVolume(
                          detectorS,        // its solid
                          vacuumMaterial,   // its material
                          "Detector");      // its name
    
  new G4PVPlacement(
                      0,                // no rotation
                      G4ThreeVector(),  // at (0,0,0)
                      detectorLV,       // its logical volume
                      "Detector",       // its name
                      worldLV,          // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps

  //
  // Air/Vacuum (Top)
  //
  fTopPV=0;
    
  auto topS
  = new G4Box("Top",           // its name
                fDetectorSizeXY/2, fDetectorSizeXY/2, topThickness/2); //its size
    
  auto topLV
  = new G4LogicalVolume(
                          topS,           // its solid
                          //airMaterial,    // its material
                          vacuumMaterial,
                          "TopLV");         // its name
    
  fTopPV = new G4PVPlacement(
                      0,                // no rotation
                      G4ThreeVector(0., 0., -fDetectorThickness/2 + topThickness/2), // its position
                      topLV,            // its logical volume
                      "Top",            // its name
                      detectorLV,       // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps
  
    
  //
  // FR4
  //
  fFR4PV=0;
    
  auto fr4S
  = new G4Box("FR4",           // its name
                fDetectorSizeXY/2, fDetectorSizeXY/2, fFr4Thickness/2); //its size
    
  auto fr4LV
  = new G4LogicalVolume(
                          fr4S,           // its solid
                          //fr4Material,    // its material
                          //kaptonMaterial,
                          //glassLeadMaterial,
                          //glassPlateMaterial,
                          cuMaterial,
                          "FR4LV");         // its name
    
  fFR4PV = new G4PVPlacement(
                      0,                // no rotation
                      G4ThreeVector(0., 0., -fDetectorThickness/2 + topThickness + fFr4Thickness/2), // its position
                      fr4LV,       // its logical volume
                      "FR4",            // its name
                      detectorLV,       // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps
    
  //
  // Air/Vacuum (Bottom)
  //
    fBottomPV=0;
    auto bottomS
    = new G4Box("Bottom",            // its name
                fDetectorSizeXY/2, fDetectorSizeXY/2, bottomThickness/2); // its size
    
    auto bottomLV
    = new G4LogicalVolume(
                          bottomS,        // its solid
                          //airMaterial, // its material
                          vacuumMaterial,
                          "BottomLV");    // its name
    
   fBottomPV = new G4PVPlacement(
                      0,                // no rotation
                      G4ThreeVector(0., 0., -fDetectorThickness/2 + topThickness + fFr4Thickness + bottomThickness/2), // its position
                      bottomLV,            // its logical volume
                      "Bottom",            // its name
                      detectorLV,       // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps
    
  //
  // Set gasLV as scoring volume
  //
  fScoringVolume = fr4LV;
    
    //
    // print parameters
    //
    G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "------------------------------------------------------------" << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> The detector is one layer of: [ "
    << topThickness/mm << "mm of " << vacuumMaterial->GetName()
    << " + "
    << fFr4Thickness/mm << "mm of " << fr4Material->GetName()
    //<< fFr4Thickness/mm << "mm of " << kaptonMaterial->GetName()
    //<< fFr4Thickness/mm << "mm of " << glassLeadMaterial->GetName()
    //<< fFr4Thickness/mm << "mm of " << glassPlateMaterial->GetName()
    //<< " + "
    //<< cuThickness/um << "um of " << cuMaterial->GetName()
    << " + "
    << bottomThickness/mm << "mm of " << vacuumMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl
    << "------------------------------------------------------------" << G4endl
    << "------------------------------------------------------------" << G4endl;

    
    //
    // Visualization attributes
    //
    worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    //detector
    auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    detectorLV->SetVisAttributes(simpleBoxVisAtt);

    // Air G4Color::White() + opacity = 0.5
    // it was air, now I have vacuum, but I'll keep it white
    auto airBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0, 0.5));
    airBoxVisAtt->SetVisibility(true);
    topLV->SetVisAttributes(airBoxVisAtt);
    bottomLV->SetVisAttributes(airBoxVisAtt);
    
    // FR4 G4Color::Red()
    //auto fr4BoxVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5, 0.5));
    auto fr4BoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.,0., 0.5));
    fr4BoxVisAtt->SetVisibility(true);
    fr4LV->SetVisAttributes(fr4BoxVisAtt);
    
    // Cu G4Color::Yellow()
    //auto cuBoxVisAtt= new G4VisAttributes(G4Color(1.0,1.0,0.0, 0.5));
    //cuBoxVisAtt->SetVisibility(true);
    //cuLV->SetVisAttributes(cuBoxVisAtt);

    //
    // Always return the physical World
    //
    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* FTMDetectorConstruction::Construct()
{  
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
    
  // Define all materials and set global variables
  // DefineMaterials() ;
  // Define volumes
  return DefineVolumes();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FTMDetectorConstruction::SetFr4Thickness(G4double val)
{
    
    fFr4Thickness = val;
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    //G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

