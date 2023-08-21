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
DRTB23SimSteppingAction::DRTB23SimSteppingAction( DRTB23SimEventAction* eventAction,
						  const DRTB23SimDetectorConstruction* detConstruction )
    : G4UserSteppingAction(),
    fEventAction(eventAction),
    fDetConstruction(detConstruction),
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
        fEventAction->AddVecTowerE(fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3)), edep );
    }

    //Collect energy in preshower (whole and scintillator)
    if ( volume == fPSPV || volume == fPSScinPV ){
        fEventAction->AddPSEnergy( edep );
    }
    if ( volume == fPSScinPV ){
        fEventAction->AddPSSciEnergy( edep );
    }
    
    //Collect energy in calo (N.B. not exact method as in 2023 the platform was included)
    if ( volume != fWorldPV && volume != fPSPV && volume != fPSScinPV ) { fEventAction->Addenergy(edep); }
 
}

//Define FastSteppingAction() method
//
void DRTB23SimSteppingAction::FastSteppingAction( const G4Step* step ) { 
		
    // Get step info
    //
    G4VPhysicalVolume* volume 
        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4double edep = step->GetTotalEnergyDeposit();
    G4double steplength = step->GetStepLength();
    
    //--------------------------------------------------
    //Store information from Scintillation and Cherenkov
    //signals
    //--------------------------------------------------
   
    std::string Fiber;
    std::string S_fiber = "S_fiber";
    std::string C_fiber = "C_fiber";
    Fiber = volume->GetName(); 
    G4int TowerID;
    G4int SiPMID = 900;
    G4int SiPMTower;
    G4int signalhit = 0;

    if ( strstr( Fiber.c_str(), S_fiber.c_str() ) ) { //scintillating fiber/tube

        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) {
            step->GetTrack()->SetTrackStatus( fStopAndKill ); 
	}

	if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return; } //not ionizing particle
//    G4VPhysicalVolume* modvolume 
//        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(3);
//	 std::cout << " grandmother name " << modvolume->GetName() << " number " << modvolume->GetCopyNo() << std::endl;
//        std::cout << " grandmother nunber " << step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3) << std::endl;

    G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);

	TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
	SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
	fEventAction->AddScin(edep);
	signalhit = fSignalHelper->SmearSSignal( fSignalHelper->ApplyBirks( edep, steplength ) );
    // Attenuate Signal
    signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);
//	if ( TowerID != 0 ) { fEventAction->AddVecSPMT( TowerID, signalhit ); }
	fEventAction->AddVecSPMT( TowerID, signalhit ); 
	if(SiPMTower > -1){ 
            SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
	    fEventAction->AddVectorScin( signalhit, SiPMTower*NoFibersTower+SiPMID ); 
        }
    }

    if ( strstr( Fiber.c_str(), C_fiber.c_str() ) ) { //Cherenkov fiber/tube

        fEventAction->AddCher(edep);

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
		    TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));		
	            SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
		    fEventAction->AddVecCPMT( TowerID, c_signal );
//		    if ( TowerID != 0 ) { fEventAction->AddVecCPMT( TowerID, c_signal ); }

		    if(SiPMTower > -1){ 
		        SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
			fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, c_signal);
	            }
		    step->GetTrack()->SetTrackStatus( fStopAndKill );
		}
		default:
		    step->GetTrack()->SetTrackStatus( fStopAndKill );
	    } //end of swich cases

        } //end of optical photon

    } //end of Cherenkov fiber
   
}

//**************************************************
