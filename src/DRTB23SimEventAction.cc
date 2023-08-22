//**************************************************
// \file DRTB23SimEventAction.cc
// \brief: Implementation of DRTB23SimEventAction 
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimEventAction.hh"
#include "DRTB23SimRunAction.hh"
#include "DRTB23SimDetectorConstruction.hh"
//Includers from Geant4
//
#include "g4root.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//Includers from C++
//
#include <iomanip>
#include <vector>
#include <numeric>

//Define constructor
//
DRTB23SimEventAction::DRTB23SimEventAction()
    : G4UserEventAction(),
    EnergyScin(0.),
    EnergyCher(0.),
    NofCherDet(0),
    NofScinDet(0),
    EnergyTot(0.),
    PrimaryPDGID(0),
    PrimaryParticleEnergy(0.),
    PrimaryX(0.),
    PrimaryY(0.),
    EscapedEnergy(0.),
    PSEnergy(0.),
    PSSciEnergy(0.),
    VectorSignals(0.),
    VectorSignalsCher(0.),
    VecSPMT(0.),
    VecCPMT(0.),
    VecTowerE(0.) {
}

//Define de-constructor
//
DRTB23SimEventAction::~DRTB23SimEventAction() {}

//Define BeginOfEventAction() and EndOfEventAction() methods
//
void DRTB23SimEventAction::BeginOfEventAction(const G4Event*) {  
    
    //Initialize data memebers at begin of each event
    //
    EnergyScin = 0.;
    EnergyCher = 0.;
    NofCherDet = 0;
    NofScinDet = 0;
    EnergyTot = 0.;
    EscapedEnergy = 0.;
    PSEnergy = 0.;
    PSSciEnergy = 0.;
    //Fields involving the primary particle
    //should NOT be cleared here,
    //they are updated by the GeneralParticleSource

    VectorSignals.clear();
    VectorSignalsCher.clear();
    VecSPMT.clear();
    VecCPMT.clear();
    VecTowerE.clear();

    VectorSignals.assign(NoFibersTower*NoModulesSiPM, 0.);
    VectorSignalsCher.assign(NoFibersTower*NoModulesSiPM, 0.);
    VecSPMT.assign(NoModulesActive, 0.);
    VecCPMT.assign(NoModulesActive, 0.);
    VecTowerE.assign(NoModulesActive, 0.);

}

void DRTB23SimEventAction::EndOfEventAction(const G4Event* ) {
 
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    //Add all p.e. in Scin and Cher fibers before calibration
    //
    for (auto& n : VectorSignals) NofScinDet += n;
    for (auto& n : VecSPMT) NofScinDet += n;
    for (auto& n : VectorSignalsCher) NofCherDet += n;
    for (auto& n : VecCPMT) NofCherDet += n;

    //Add to EnergyTot the energies in towers
    //
    EnergyTot = std::accumulate(VecTowerE.begin(),VecTowerE.end(),0.);

    //Fill ntuple event by event
    //entries with vectors are automatically filled
    //
    analysisManager->FillNtupleDColumn(0, EnergyScin);
    analysisManager->FillNtupleDColumn(1, EnergyCher);
    analysisManager->FillNtupleDColumn(2, NofCherDet);
    analysisManager->FillNtupleDColumn(3, NofScinDet);
    analysisManager->FillNtupleDColumn(4, EnergyTot);
    analysisManager->FillNtupleDColumn(5, PrimaryParticleEnergy);
    analysisManager->FillNtupleIColumn(6, PrimaryPDGID);
    analysisManager->FillNtupleDColumn(7, EscapedEnergy);
    analysisManager->FillNtupleDColumn(8, PSEnergy);
    analysisManager->FillNtupleDColumn(9, PSSciEnergy);
    analysisManager->FillNtupleDColumn(10,PrimaryX);
    analysisManager->FillNtupleDColumn(11,PrimaryY);
    analysisManager->AddNtupleRow();
    //Vector entries in ntuple are automatically filled

}

//**************************************************
