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
// $Id: B1DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef FTMDetectorConstruction_h
#define FTMDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

class FTMDetectorMessenger;

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction class to define materials and geometry.
class FTMDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    FTMDetectorConstruction(); // Constructor
    virtual ~FTMDetectorConstruction(); // Destructor
    
  public:
    
    void SetFr4Thickness(G4double val);

    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    
    G4double GetFr4Thickness()         {return fFr4Thickness;};
    G4double GetWorldSizeXY()          {return fWorldSizeXY;};
    G4double GetWorldSizeZ()           {return fWorldSizeZ;};
    G4double GetDetectorSizeXY()       {return fDetectorSizeXY;};
    G4double GetDetectorSizeZ()        {return fDetectorThickness;};
    
    const G4VPhysicalVolume* GetFr4()   {return fFR4PV;};
    const G4VPhysicalVolume* GetTop()   {return fTopPV;};
    const G4VPhysicalVolume* GetBottom()   {return fBottomPV;};

    
    G4LogicalVolume*   fScoringVolume;
    G4VPhysicalVolume* fFR4PV;
    G4VPhysicalVolume* fTopPV;
    G4VPhysicalVolume* fBottomPV;
    
  private:
    
    void DefineMaterials() ;
    G4VPhysicalVolume* DefineVolumes();
    
    void DefineCommands();
    
    FTMDetectorMessenger* fDetectorMessenger;  //pointer to the Messenger

    G4Material*        vacuumMaterial;
    G4Material*        airMaterial;
    G4Material*        fr4Material;
    G4Material*        cuMaterial;
    G4Material*        kaptonMaterial;
    G4Material*        glassLeadMaterial;
    G4Material*        glassPlateMaterial;
    G4Material*        gasMaterial;
    G4bool             fCheckOverlaps; // option to activate checking of volumes overlaps
    G4int              fNofLayers;     // number of layers
    G4double           fFr4Thickness;
    G4double           fWorldSizeXY;
    G4double           fWorldSizeZ;
    G4double           fDetectorSizeXY;
    G4double           fDetectorThickness;
    

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

