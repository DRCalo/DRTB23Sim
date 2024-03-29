//**************************************************
// \file DRTB23SimRunAction.cc 
// \brief: Implementation of 
//         DRTB23SimRunAction class 
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimRunAction.hh"
#include "DRTB23SimEventAction.hh"

//Includers from Geant4
//
#include "g4root.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//Includers from C++
//
#include <string>

//Define constructor
//
DRTB23SimRunAction::DRTB23SimRunAction( DRTB23SimEventAction* eventAction )
    : G4UserRunAction(),
      fEventAction( eventAction ){ 
  
    //print event number per each event (default, can be overwritten with macro)
    //
    G4RunManager::GetRunManager()->SetPrintProgress(1);     

    //Instantiate analysis manager
    //
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel( 1 );
    analysisManager->SetNtupleMerging( 1 );

    //Using ROOT as analysisManager type, print it
    //
    G4cout << "DRTB23Sim-> Using " << analysisManager->GetType() << G4endl;

    //Define ntuple structure
    //
    analysisManager->CreateNtuple("DRTB23Simout", "simoutput");
    analysisManager->CreateNtupleDColumn("EnergyScin");                     //0
    analysisManager->CreateNtupleDColumn("EnergyCher");                     //1
    analysisManager->CreateNtupleDColumn("NofCherDet");                     //2
    analysisManager->CreateNtupleDColumn("NofScinDet");                     //3
    analysisManager->CreateNtupleDColumn("EnergyTot");                      //4
    analysisManager->CreateNtupleDColumn("PrimaryParticleEnergy");          //5
    analysisManager->CreateNtupleIColumn("PrimaryPDGID");                   //6
    analysisManager->CreateNtupleDColumn("EscapedEnergy");                  //7
    analysisManager->CreateNtupleDColumn("PSEnergy");                       //8
    analysisManager->CreateNtupleDColumn("PSSciEnergy");                    //9
    analysisManager->CreateNtupleDColumn("PrimaryX");                       //10
    analysisManager->CreateNtupleDColumn("PrimaryY");                       //11
    analysisManager->CreateNtupleDColumn("NofSiPMScinDet");                 //12
    analysisManager->CreateNtupleDColumn("NofSiPMCherDet");                 //13
    analysisManager->CreateNtupleDColumn("VectorSignals", fEventAction->GetVectorSignals());
    analysisManager->CreateNtupleDColumn("VectorSignalsCher", fEventAction->GetVectorSignalsCher());
    analysisManager->CreateNtupleDColumn("VecTowerE", fEventAction->GetVecTowerE());
    analysisManager->CreateNtupleDColumn("VecSPMT", fEventAction->GetVecSPMT());
    analysisManager->CreateNtupleDColumn ("VecCPMT", fEventAction->GetVecCPMT());
    analysisManager->FinishNtuple();
      
}

//Define de-constructor
//
DRTB23SimRunAction::~DRTB23SimRunAction(){
   
    //Delete only instance of G4AnalysisManager
    //
    delete G4AnalysisManager::Instance();  

}

//Define BeginOfRunAction() and EndOfRunAction() methods
//
void DRTB23SimRunAction::BeginOfRunAction( const G4Run* Run )  { 
    
    //Save random seeds (optional)
    //
    //G4RunManager::GetRunManager()->SetRandomNumberStore( true );
    
    //Open output file, one per Run
    //
    auto analysisManager = G4AnalysisManager::Instance();
    std::string runnumber = std::to_string( Run->GetRunID() );
    G4String outputfile = "DRTB23Sim_Run"+runnumber;
    analysisManager->OpenFile( outputfile );

}

void DRTB23SimRunAction::EndOfRunAction( const G4Run* ) {
  
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->Write();
    analysisManager->CloseFile();

}

//**************************************************
