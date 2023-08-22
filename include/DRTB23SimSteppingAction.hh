//**************************************************
// \file DRTB23SimSteppingAction.hh
// \brief: Definition of DRTB23SimSteppingAction.hh
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef DRTB23SimSteppingAction_h
#define DRTB23SimSteppingAction_h 1

//Includers from Geant4
//
#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4LogicalVolume.hh"

//Forward declarations from Geant4
//
class G4OpBoundaryProcess;

//Forward declarations from project files
//
class DRTB23SimDetectorConstruction;
class DRTB23SimEventAction;

//Includers from project files
//
#include "DRTB23SimSignalHelper.hh"

class DRTB23SimSteppingAction : public G4UserSteppingAction {
    
    public:
        //Constructor
        //
        DRTB23SimSteppingAction(DRTB23SimEventAction* eventAction,
				const DRTB23SimDetectorConstruction* detConstruction );
        //De-constructor
        //
        virtual ~DRTB23SimSteppingAction();
        
        //User impementation of SteppingAction
        //
        virtual void UserSteppingAction( const G4Step* step );

        //Retrieve auxialiry info from Step
        //
        void AuxSteppingAction( const G4Step* step );

        //Fast signal simulation (no optical photon propagation)
        //
        void FastSteppingAction( const G4Step* step ); 

    private:

        DRTB23SimEventAction*  fEventAction;  

        G4OpBoundaryProcess* fOpProcess;
                
	//Pointers
	//
        const DRTB23SimDetectorConstruction* fDetConstruction;
	G4VPhysicalVolume* fWorldPV; //PV: world volume
        G4VPhysicalVolume* fPSPV; //PV: preshower volume
        G4VPhysicalVolume* fPSScinPV; //PV: preshower scintillator volume
        G4LogicalVolume*   fSfiber_Abs_LV; //LV: Absorber of S fiber
        G4LogicalVolume*   fSfiber_Core_LV; //LV: Core of S fiber
        G4LogicalVolume*   fSfiber_Clad_LV; //LV: Cladding of S fiber
        G4LogicalVolume*   fCfiber_Abs_LV; //LV: Absorber of C fiber
        G4LogicalVolume*   fCfiber_Core_LV; //LV: Core of C fiber
        G4LogicalVolume*   fCfiber_Clad_LV; //LV: Cladding of C fiber
				
        //Pointer to only existing implementation (singleton)
    	//of DRTB23SimTowerHelper
    	//
        DRTB23SimSignalHelper* fSignalHelper;

};

#endif

//**************************************************
