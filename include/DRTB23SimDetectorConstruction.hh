//**************************************************
// \file DRTB23SimDetectorConstruction.hh
// \brief: Definition of 
//         DRTB23SimDetectorConstruction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
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

//Includers from project files
//
#include "DRTB23SimGeoPar.hh"
#include "DRTB23SimGeoMessenger.hh"

//Forward declaration
//
class G4VPhysicalVolume;

class DRTB23SimDetectorConstruction : public G4VUserDetectorConstruction {
  
    public:
        //Constructor
        //
        DRTB23SimDetectorConstruction();

        //De-constructor
        //
        virtual ~DRTB23SimDetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct();

        G4LogicalVolume* constructscinfiber(G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* ScinMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CherMaterial);
    
        G4LogicalVolume* constructcherfiber(G4double tuberadius,
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
    	G4VPhysicalVolume* GetWorldPV() const;
    	G4VPhysicalVolume* GetPSPV() const;
    	G4VPhysicalVolume* GetPSScinPV() const;
    	G4LogicalVolume* GetSAbsLV() const {return fSfiber_Abs_LV;};
    	G4LogicalVolume* GetSCoreLV() const {return fSfiber_Core_LV;};
    	G4LogicalVolume* GetSCladLV() const {return fSfiber_Clad_LV;};
    	G4LogicalVolume* GetCAbsLV() const {return fCfiber_Abs_LV;};
    	G4LogicalVolume* GetCCoreLV() const {return fCfiber_Core_LV;};
    	G4LogicalVolume* GetCCladLV() const {return fCfiber_Clad_LV;};
        G4double GetXshift() const {return fXshift;};
        G4double GetYshift() const {return fYshift;};
        G4double GetOrzrot() const {return fOrzrot;};
        G4double GetVerrot() const {return fVerrot;};

        //Setters
        //
        void SetXshift(const G4double& val) {fXshift=val;};
        void SetYshift(const G4double& val) {fYshift=val;};
        void SetOrzrot(const G4double& val) {fOrzrot=val;};
        void SetVerrot(const G4double& val) {fVerrot=val;};

        //Other methods
	//
	G4int GetTowerID( const G4int& cpno ) const;
        G4int GetSiPMID(const G4int& cpno ) const; 
	G4int GetSiPMTower(const G4int& town ) const;
       
	//Build contour in x-y plane of a module as 
	//an hexcell shape
	//
	std::vector<G4TwoVector> calcmod(double radius, int nrow, int ncol); 

    private:
        
        //Mandatory method for Geant4
        //
        G4VPhysicalVolume* DefineVolumes();
	
        //Members
	//
        G4bool  fCheckOverlaps; // option for checking volumes overlaps
				
	G4VPhysicalVolume* fWorldPV; //PV: world volume
        G4VPhysicalVolume* fPSPV; //PV: preshower volume
        G4VPhysicalVolume* fPSScinPV; //PV: preshower scintillator volume
        G4LogicalVolume*   fSfiber_Abs_LV; //LV: Absorber of S fiber
        G4LogicalVolume*   fSfiber_Core_LV; //LV: Core of S fiber
        G4LogicalVolume*   fSfiber_Clad_LV; //LV: Cladding of S fiber
        G4LogicalVolume*   fCfiber_Abs_LV; //LV: Absorber of C fiber
        G4LogicalVolume*   fCfiber_Core_LV; //LV: Core of C fiber
        G4LogicalVolume*   fCfiber_Clad_LV; //LV: Cladding of C fiber

        //Pointer to messenger for UI
        //
        DRTB23SimGeoMessenger* fGeoMessenger;

        //Parameters selectable via UI
        //
        G4double fXshift{0.}, fYshift{0.}, fVerrot{0.}, fOrzrot{0.};
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

inline G4VPhysicalVolume* DRTB23SimDetectorConstruction::GetWorldPV() const {
    return fWorldPV;
}

inline G4VPhysicalVolume* DRTB23SimDetectorConstruction::GetPSPV() const {
    return fPSPV;
}

inline G4VPhysicalVolume* DRTB23SimDetectorConstruction::GetPSScinPV() const {
    return fPSScinPV;
}

#endif

//**************************************************
