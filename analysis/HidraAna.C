//**************************************************
// \file HidraAna.C
// \brief:  analysis skeleton for HidraSim ntuples
// \author: Giacomo Polesello (INFN Pavia) 
//          Edoardo Proserpio (Uni Insubria)
// \start date: May 3, 2022
//**************************************************
//
////usage: root -l -b -q 'HidraAna.C(energy,"filename")'
///  where energy is the energy of the beam, and filename
//   the name of the data ntuple
//
//   It produces an histogram file 
//   Hidra+"energy"+.root
//
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <array>
#include <stdint.h>
#include <string>
#include <fstream>
#include <string>
#include <cstring>
// include file with geometry of module
#include "HidraGeo.h"
void HidraAna(double energy, const string intup){
//Open ntuples
  string infile = "/home/storage/hidradata/"+intup;
  std::cout<<"Using file: "<<infile<<std::endl;
  char cinfile[infile.size() + 1];
  strcpy(cinfile, infile.c_str());
  auto simfile = new TFile(cinfile, "READ");
  auto *simtree = (TTree*)simfile->Get( "DREMTubesout" );
  std::ostringstream os;
  os << energy;
  std::string enstr = os.str();
  string outfile="hidra"+enstr+".root";
  TFile f(outfile.c_str(), "RECREATE");
//
//  build vectors with row and column position of each of
//  the MiniModules
//
  int modcol[84];
  int modrow[84];
  for(int i=0;i<120;i++){
    int row=i/NofmodulesX;
    int col=i%NofmodulesX;
    int imod=modflag[i];
    if(imod>=0){
      modcol[imod]=col;
      modrow[imod]=row;
    }
  }
// book histograms  
  double bmin=energy-0.4*sqrt(energy)*10.;
  double bmax=energy+0.4*sqrt(energy)*10.;
  auto sciene = new TH1F("sciene", "sciene",100,bmin,bmax);
  auto cerene = new TH1F("cerene", "cerene",100,bmin,bmax);
  auto totene = new TH1F("totene", "totene",100,bmin,bmax);
  auto totenec = new TH1F("totenec", "totenec",100,bmin,bmax);
  auto totdep = new TH1F("totdep", "totdep",100,bmin,bmax);
  auto leakene = new TH1F("leakene", "leakene",100,0.,0.1);
  auto chidist = new TH1F("chidist", "chidist",100,0.,1.);
  auto mapcalo  = new TH2F("mapcalo", "mapcalo",NofmodulesX,0.,NofmodulesX,NofmodulesY,0.,NofmodulesY);

  int nentries=simtree->GetEntries();
  std::cout<<"Entries "<<nentries<<std::endl;

//Allocate branch pointers
  int pdg; simtree->SetBranchAddress( "PrimaryPDGID", &pdg );
  double venergy; simtree->SetBranchAddress( "PrimaryParticleEnergy", &venergy );
  double lenergy; simtree->SetBranchAddress( "EscapedEnergy", &lenergy );
  double edep; simtree->SetBranchAddress( "EnergyTot", &edep );
  double Stot; simtree->SetBranchAddress( "NofScinDet", &Stot );
  double Ctot; simtree->SetBranchAddress( "NofCherDet", &Ctot );
  double PSdep; simtree->SetBranchAddress( "PSEnergy", &PSdep );
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
// 
  double chi=0.38;   
  double sciPheGeV=217.501;
  double cerPheGeV=54.1621;
  double elcont=1.005;
  double picont=1.028;
// Loop on events 
  for( unsigned int i=0; i<simtree->GetEntries(); i++){
    double ecalo=energy-lenergy/1000;
    simtree->GetEntry(i);
    double totsci=0.;
    double totcer=0.;
    double tottow=0.;
// Sum energy over all MiniModules
    for(unsigned int j=0; j<SPMT->size(); j++){
      totsci+=SPMT->at(j)/sciPheGeV;
      totcer+=CPMT->at(j)/cerPheGeV;
      tottow+=TowerE->at(j);
      mapcalo->Fill(modcol[j],modrow[j],TowerE->at(j)/1000/nentries);
    }
    sciene->Fill(totsci);        
    cerene->Fill(totcer);        
    totene->Fill(elcont*0.5*(totsci+totcer));   
    totenec->Fill(picont*(totsci-chi*totcer)/(1-chi));   
    totdep->Fill(tottow/1000.);   
    leakene->Fill(lenergy/1000/energy);   
    chidist->Fill((totsci-ecalo)/(totcer-ecalo));    
  }
  totenec->Fit("gaus","Q","");
  TF1 *fit1 = totenec->GetFunction("gaus");
  double peak1=fit1->GetParameter(1);
  double epeak1=fit1->GetParError(1);
  double rms1=fit1->GetParameter(2);
  double erms1=fit1->GetParError(2);
  cout << " # " << energy << " " << peak1 << " " << epeak1 << " " << rms1 << " " << erms1 << endl;
  f.Write();
  //

}

//**************************************************
