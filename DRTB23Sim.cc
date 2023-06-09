//**************************************************
// \file DRTB23Sim.cc
// \brief: main() of DRTB23Sim project
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
// @lopezzot
// \start date: 9 June 2023
//**************************************************

// Includers from project files
//
#include "DREMTubesDetectorConstruction.hh"
#include "DREMTubesActionInitialization.hh"
#include "DREMTubesPhysicsList.hh"

// Includers from Geant4
//
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// G4err output for usage error
//
namespace PrintUsageError {
    void UsageError() {
    G4cerr << "->DREMTubes usage: " << G4endl;
    G4cerr << "DREMTubes [-m macro ] [-u UIsession] [-t nThreads] [-pl PhysicsList]" 
        << G4endl;
    G4cerr << "          [-opt FullOptic]" << G4endl;
    }
}

// main() function
//
int main(int argc, char** argv) {
    
    // Error in argument numbers
    //
    if ( argc > 11 ) {
        PrintUsageError::UsageError();
        return 1;
    }
  
    // Convert arguments in G4string and G4int
    //
    G4String macro;
    G4String session;
    G4String custom_pl = "FTFP_BERT"; //default physics list
    G4bool FullOptic = false;
    #ifdef G4MULTITHREADED
    G4int nThreads = 0;
    #endif

    for ( G4int i=1; i<argc; i=i+2 ) {
        if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
        else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
        else if ( G4String(argv[i]) == "-pl") custom_pl = argv[i+1];
        else if ( G4String(argv[i]) == "-opt") FullOptic =  
                                            G4UIcommand::ConvertToBool(argv[i+1]);
        #ifdef G4MULTITHREADED
        else if ( G4String(argv[i]) == "-t" ) {
            nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
        }
        #endif
        else {
            PrintUsageError::UsageError();
        return 1;
        }
    } 

    //Print if FullOptic option is on
    //
    if (FullOptic){ G4cout<<"DREMTubes-> Run with full optical description"<<G4endl; } 
  
    // Detect interactive mode (if no macro provided) and define UI session
    //
    G4UIExecutive* ui = nullptr;
    if ( ! macro.size() ) { //if macro card is none
        ui = new G4UIExecutive(argc, argv, session);
    }

    // Construct the default run manager
    //
    #ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    if ( nThreads > 0 ) { 
        runManager->SetNumberOfThreads(nThreads);
    }  
    #else
    G4RunManager * runManager = new G4RunManager;
    #endif

    // Set mandatory initialization classes
    //
    auto DetConstruction = new DREMTubesDetectorConstruction();
    runManager->SetUserInitialization(DetConstruction);

    runManager->SetUserInitialization(new DREMTubesPhysicsList(custom_pl, FullOptic ));
  
    auto actionInitialization = new DREMTubesActionInitialization( DetConstruction, FullOptic );
    runManager->SetUserInitialization(actionInitialization);
  
    // Initialize visualization
    //
    auto visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    auto UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    //
    if ( macro.size() ) { //macro card mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command+macro);
    }
    else  { //start UI session  
    UImanager->ApplyCommand("/control/execute DRTB23Sim_init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute DRTB23Sim_gui.mac");
    }
    ui->SessionStart();
    delete ui;
    }

    // Program termination (user actions deleted by run manager)
    //
    delete visManager;
    delete runManager;

}

//**************************************************
