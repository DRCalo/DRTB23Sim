//**************************************************
// \file DRTB23SimDetectorConstruction.hh
// \brief: Definition of DRTB23SimDetectorConstruction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multplie times
//
#ifndef DRTB23SimDetectorConstruction_h
#define DRTB23SimDetectorConstruction_h 1

//Includers from Geant4
//
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

//Include geometrical parameters
#include "DRTB23SimGeoPar.hh"
//Forward declaration
//
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class DRTB23SimDetectorConstruction : public G4VUserDetectorConstruction {
  
    public:
        //Constructor
        //
        DRTB23SimDetectorConstruction(const G4bool VerRot);
        //De-constructor
        //
        virtual ~DRTB23SimDetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct();

        G4LogicalVolume* constructscinfiber(double tolerance,
                                            G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* ScinMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CherMaterial);
    
        G4LogicalVolume* constructcherfiber(double tolerance, 
                                            G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* CherMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CladCherMaterial);
        //Getters
    	//
	const G4VPhysicalVolume* GetLeakCntPV() const;
    	const G4VPhysicalVolume* GetWorldPV() const;

        //Other methods
	//
	G4int GetTowerID( const G4int& cpno ) const;
        G4int GetSiPMID(const G4int& cpno ) const; 
	G4int GetSiPMTower(const G4int& town ) const;
       
        //
	//  Build contour in x-y plane of a module as 
	//  an hexcell shape
	//
	std::vector<G4TwoVector> calcmod(double radius, int nrow, int ncol); 

    private:
        
        //Mandatory method for Geant4
        //
        G4VPhysicalVolume* DefineVolumes();
/*
	void DefineCommands();
	G4GenericMessenger* fMessenger;
	G4double fAngleX;
	G4double fAngleY;
*/				//Members
				//
        G4bool  fCheckOverlaps; // option for checking volumes overlaps
				
				G4VPhysicalVolume* fLeakCntPV; //PV: lekage counter
				G4VPhysicalVolume* fWorldPV;   //PV: wourld volume

        G4bool fVertRot;  
};

inline G4int DRTB23SimDetectorConstruction::GetTowerID( const G4int& cpno ) const {
// remap as for 2021 hardware numbering from front face
//    const G4int idmap[9]={1,2,3,4,0,5,6,7,8};
// test:remap as for output of old simulation
//    const G4int idmap[9]={3,2,1,5,0,4,8,7,6};
//    return idmap[cpno-1];
    return cpno;		
}

inline G4int DRTB23SimDetectorConstruction::GetSiPMID( const G4int& cpno ) const {
// kept for compatibility with old simulation. Dummy for now
    return cpno;		
}

inline G4int DRTB23SimDetectorConstruction::GetSiPMTower( const G4int& town ) const {
    G4int SiPMTower=-1;
    for(int i=0;i<NoModulesSiPM;i++){
      if(town==SiPMMod[i])SiPMTower=i;
    }
    return SiPMTower;		
}

inline const G4VPhysicalVolume* DRTB23SimDetectorConstruction::GetLeakCntPV() const {
    return fLeakCntPV;
}

inline const G4VPhysicalVolume* DRTB23SimDetectorConstruction::GetWorldPV() const {
    return fWorldPV;
}

#endif

//**************************************************
