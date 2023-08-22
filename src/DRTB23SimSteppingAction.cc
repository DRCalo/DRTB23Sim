//**************************************************
// \file DRTB23SimSteppingAction.cc
// \brief: Implementation of 
//         DRTB23SimSteppingAction.cc
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimSteppingAction.hh"
#include "DRTB23SimEventAction.hh"
#include "DRTB23SimDetectorConstruction.hh"

//Includers from Geant4
//
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"

//Define constructor
//
DRTB23SimSteppingAction::DRTB23SimSteppingAction( DRTB23SimEventAction* eventAction )
    : G4UserSteppingAction(),
    fEventAction(eventAction),
    fWorldPV(nullptr),
    fPSPV(nullptr),
    fPSScinPV(nullptr),
    fSfiber_Abs_LV(nullptr),
    fSfiber_Core_LV(nullptr),
    fSfiber_Clad_LV(nullptr),
    fCfiber_Abs_LV(nullptr),
    fCfiber_Core_LV(nullptr),
    fCfiber_Clad_LV(nullptr){
	
        fSignalHelper = DRTB23SimSignalHelper::Instance(); 
	
}

//Define de-constructor
//
DRTB23SimSteppingAction::~DRTB23SimSteppingAction() {}

//Define UserSteppingAction() method
//
void DRTB23SimSteppingAction::UserSteppingAction( const G4Step* step ) {
    
    //Save auxiliary information
    //(comment out AuxSteppingAction for fast execution)
    AuxSteppingAction( step );

    //Save signals
    FastSteppingAction( step );
}

//Define AuxSteppingAction() method
//
void DRTB23SimSteppingAction::AuxSteppingAction( const G4Step* step ) {

    //--------------------------------------------------
    //Store auxiliary information from event steps
    //--------------------------------------------------

    //Get volumes of interest
    //
    if(!fWorldPV){
        const DRTB23SimDetectorConstruction* detector = 
            static_cast<const DRTB23SimDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fWorldPV = detector->GetWorldPV();
        fPSPV = detector->GetPSPV();
        fPSScinPV = detector->GetPSScinPV();
        fSfiber_Abs_LV = detector->GetSAbsLV();
        fSfiber_Core_LV = detector->GetSCoreLV();
        fSfiber_Clad_LV = detector->GetSCladLV();
        fCfiber_Abs_LV = detector->GetCAbsLV();
        fCfiber_Core_LV = detector->GetCCoreLV();
        fCfiber_Clad_LV = detector->GetCCladLV();
    }

    // Get step info
    //
    G4VPhysicalVolume* volume 
        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4LogicalVolume* Lvolume = volume->GetLogicalVolume();
    G4double edep = step->GetTotalEnergyDeposit();
    
    // Collect out of world leakage
    //
    if (!step->GetTrack()->GetNextVolume()) {
        fEventAction->AddEscapedEnergy(step->GetTrack()->GetKineticEnergy());
    }

    if ( Lvolume== fSfiber_Clad_LV || Lvolume== fSfiber_Core_LV || Lvolume== fSfiber_Abs_LV  ||
         Lvolume== fCfiber_Clad_LV || Lvolume== fCfiber_Core_LV || Lvolume== fCfiber_Abs_LV ) {
        fEventAction->AddVecTowerE(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3), edep );
        if (Lvolume==fSfiber_Core_LV) fEventAction->AddScin(edep);
        else if (Lvolume==fCfiber_Core_LV) fEventAction->AddCher(edep);
        else {}
    }

    //Collect energy in preshower (whole and scintillator)
    if ( volume == fPSPV || volume == fPSScinPV ){
        fEventAction->AddPSEnergy( edep );
    }
    if ( volume == fPSScinPV ){
        fEventAction->AddPSSciEnergy( edep );
    }
 
}

//Define FastSteppingAction() method
//
void DRTB23SimSteppingAction::FastSteppingAction( const G4Step* step ) { 
		
    //-----------------------------------------------------
    //Store signals from Scintillation and Cherenkov fibers
    //-----------------------------------------------------

    //Get volumes of interest, i.e. core of fibers
    //
    if(!fWorldPV){
        const DRTB23SimDetectorConstruction* detector = 
            static_cast<const DRTB23SimDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fSfiber_Core_LV = detector->GetSCoreLV();
        fCfiber_Core_LV = detector->GetCCoreLV();
    }

    // Get step info
    //
    G4LogicalVolume* Lvolume 
        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4double edep = step->GetTotalEnergyDeposit();
    G4double steplength = step->GetStepLength();

    //Return if the step is not in core of fibers
    //
    if( Lvolume != fSfiber_Core_LV && Lvolume != fCfiber_Core_LV) return;

    else if (Lvolume == fSfiber_Core_LV){ //scintillating fiber
 
        G4int signalhit = 0;

        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) {
            step->GetTrack()->SetTrackStatus( fStopAndKill ); 
            return;
        }

	if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. 
             || edep == 0. ) { return; }

        //G4VPhysicalVolume* modvolume 
        //    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(3);
        //G4cout << " grandmother name " << modvolume->GetName() << " number " << modvolume->GetCopyNo() << G4endl;
        //G4cout << " grandmother nunber " << step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3) << G4endl;

        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);

        G4int TowerID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);
	signalhit = fSignalHelper->SmearSSignal( fSignalHelper->ApplyBirks( edep, steplength ) );
        // Attenuate Signal
        signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);
	fEventAction->AddVecSPMT( TowerID, signalhit ); 
	if(TowerID == 0){ // in sipm-readout tower
            G4int SiPMID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
	    fEventAction->AddVectorScin( signalhit, SiPMID ); 
        }
    } // end of scintillating fiber

    else if (Lvolume == fCfiber_Core_LV ) { //Cherenkov fiber

	if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ){
	
	    G4OpBoundaryProcessStatus theStatus = Undefined;

	    G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

	    if (OpManager) {
    	        G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
		G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

		for ( G4int i=0; i<MAXofPostStepLoops; i++) {
		    G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
		    fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
		    if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break; }
		}
	    }

	    switch ( theStatus ){
						
	        case TotalInternalReflection: {
                    G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
		    G4int c_signal = fSignalHelper->SmearCSignal( );
                    // Attenuate Signal
                    c_signal = fSignalHelper->AttenuateCSignal(c_signal, distance_to_sipm);
		    G4int TowerID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);		
		    fEventAction->AddVecCPMT( TowerID, c_signal );

		    if(TowerID == 0){ // in sipm-readout tower
		        G4int SiPMID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
			fEventAction->AddVectorCher(SiPMID, c_signal);
	            }
		    step->GetTrack()->SetTrackStatus( fStopAndKill );
		}
		default: 
                    /*step->GetTrack()->SetTrackStatus( fStopAndKill )*/;
	    } //end of swich cases
        } //end of optical photon
        else return;
    } //end of Cherenkov fiber

    else return;

}

//**************************************************
