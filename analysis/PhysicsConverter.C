//**************************************************
// \file PhysicsConverter.C
// \brief: converter from merged trees to Event obj
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
//          Edoardo Proserpio (Uni Insubria)
// \start date: 20 August 2021
//**************************************************
//
////usage: root -l .x PhysicsConverter.C++
//
//
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <array>
#include <stdint.h>
#include <string>
#include <fstream>
#include "PhysicsEvent.h"
#include <string>
#include <cstring>
//
//  Script which takes the output of TB23SIM 
//  and converts it into an ntuple in the same format
//  as the data
//
//  The name of the output file is
//   
//    physics_sps2023_runxxxxx.root 
//
//    where xxxxx is a dummy run number given in calling sequence
//
ClassImp(EventOut)

void PhysicsConverter(const string run){

  //Open merge ntuples
  //
  string infile = "DRTB23Sim_Run0.root";
  std::cout<<"Using file: "<<infile<<std::endl;
  char cinfile[infile.size() + 1];
  strcpy(cinfile, infile.c_str());
  string outfile = "physics_sps2023_run"+run+".root";
  char coutfile[outfile.size() + 1];
  strcpy(coutfile, outfile.c_str());
  auto simfile = new TFile(cinfile, "READ");
  auto *simtree = (TTree*)simfile->Get( "DRTB23Simout" );
  //Create new tree and Event object
  //
  auto Outfile = new TFile(coutfile,"RECREATE");
  auto ftree = new TTree("Ftree","Ftree");
  ftree->SetDirectory(Outfile);
  auto ev = new Event();
  auto evout = new EventOut();
  ftree->Branch("Events",evout);
  //Create calibration objects
  //
  SiPMCalibration sipmCalibration;
  PMTCalibration pmtCalibration;

  //Check entries in trees
  //
  std::cout<<"Entries in PMT "<<simtree->GetEntries()<<std::endl;

  //Allocate branch pointers
        int pdg; simtree->SetBranchAddress( "PrimaryPDGID", &pdg );
        double venergy; simtree->SetBranchAddress( "PrimaryParticleEnergy", &venergy );
        double lenergy; simtree->SetBranchAddress( "EscapedEnergy", &lenergy );
        double edep; simtree->SetBranchAddress( "EnergyTot", &edep );
        double Stot; simtree->SetBranchAddress( "NofScinDet", &Stot );
        double Ctot; simtree->SetBranchAddress( "NofCherDet", &Ctot );
        double PSdep; simtree->SetBranchAddress( "PSEnergy", &PSdep );
//        double PSScidep; simtree->SetBranchAddress( "PShower", &PSScidep );
        double beamX; simtree->SetBranchAddress( "PrimaryX", &beamX );
        double beamY; simtree->SetBranchAddress( "PrimaryY", &beamY );
        vector<double>* TowerE = NULL; 
        simtree->SetBranchAddress( "VecTowerE", &TowerE );
        vector<double>* SPMT = NULL; 
        simtree->SetBranchAddress( "VecSPMT", &SPMT );
        vector<double>* CPMT = NULL; 
        simtree->SetBranchAddress( "VecCPMT", &CPMT );
        vector<double>* SSiPM = NULL; 
        simtree->SetBranchAddress( "VectorSignals", &SSiPM );
        vector<double>* CSiPM = NULL; 
        simtree->SetBranchAddress( "VectorSignalsCher", &CSiPM );
        
  for( unsigned int i=0; i<simtree->GetEntries(); i++){
    simtree->GetEntry(i);
    evout->EventID = 0;

    //Fill ev data members
    //
    ev->SPMT1 = SPMT->at(3);
    ev->SPMT2 = SPMT->at(2);
    ev->SPMT3 = SPMT->at(1);
    ev->SPMT4 = SPMT->at(5);
    ev->SPMT5 = SPMT->at(4);
    ev->SPMT6 = SPMT->at(8);
    ev->SPMT7 = SPMT->at(7);
    ev->SPMT8 = SPMT->at(6);
    ev->CPMT1 = CPMT->at(3);
    ev->CPMT2 = CPMT->at(2);
    ev->CPMT3 = CPMT->at(1);
    ev->CPMT4 = CPMT->at(5);
    ev->CPMT5 = CPMT->at(4);
    ev->CPMT6 = CPMT->at(8);
    ev->CPMT7 = CPMT->at(7);
    ev->CPMT8 = CPMT->at(6);
//
//    evout->PShower = PSScidep*69.8151951+210;
    evout->PShower = 0.;
    evout->PShowerall = PSdep;
    evout->MCounter = 0.;
    evout->C1 = 90.;
    evout->C2 = 30.;
    evout->XDWC1=beamX;
    evout->YDWC1=beamY;
    evout->XDWC2=beamX;
    evout->YDWC2=beamY;
    evout->EnergyTot=edep;
    evout->EscapedEnergy=lenergy;
    for( unsigned int j=0; j<9; j++){
      evout->VecTowerE[j]=TowerE->at(j);
    } 
    for( unsigned int j=0; j<160; j++){
      int column=j/10;
      int col1=15-column;
      int row=j%10;
      int row1=9-row;
      int j1=16*row1+col1;
      if(j1>160){
	 std::cout << j1 << " " << row1 << " " << " " << row << " " << col1 << " " << column << endl;
      }
      ev->SiPM_sci[j1]=SSiPM->at(j);
      ev->SiPM_che[j1]=CSiPM->at(j);
    }
    //Calibrate SiPMs and PMTs
    //
    ev->calibrate(sipmCalibration, evout);
    ev->calibratePMT(pmtCalibration, evout);
    evout->CompSPMTene();
    evout->CompCPMTene();
    //Write event in ftree
    //
    ftree->Fill();
    //Reset totSiPMPheC and totSiPMPheS to 0
    //
    evout->totSiPMCene = 0;
    evout->totSiPMSene = 0;
  }

  //Write and close Outfile
  //
  ftree->Write();
  Outfile->Close();

}

//**************************************************
