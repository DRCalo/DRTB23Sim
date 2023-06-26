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
				const DRTB23SimDetectorConstruction* detConstruction,
                                const G4bool FullOptic );
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
        //fFullOptic == false
        //
        void FastSteppingAction( const G4Step* step ); 

        //Slow signal simulation (optical photon propagation)
        //fFullOptic == true
        //
        void SlowSteppingAction( const G4Step* step );
    
    private:

        DRTB23SimEventAction*  fEventAction;  

        G4OpBoundaryProcess* fOpProcess;
                
	//Pointer to DRTB23SimDetectorConstruction
	//
        const DRTB23SimDetectorConstruction* fDetConstruction;
				
	G4bool fFullOptic;

        //Pointer to only existing implementation (singleton)
    	//of DRTB23SimTowerHelper
    	//
        DRTB23SimSignalHelper* fSignalHelper;

};

#endif

//**************************************************
