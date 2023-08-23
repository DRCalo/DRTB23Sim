//**************************************************
// \file DRTB23SimConverter.C
// \brief: converter from DRTB23Sim format to
//         2023 test-beam format
//\history: - First implementation by Giacomo 
//          Polesello during 2023 July test-beam.
//          - Adaptation by Lorenzo Pezzotti
//          on 23 August 2023.
//**************************************************
//
////usage: root .x 'DRTB23SimConverter.C("runNo")'
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
double const sq3=sqrt(3.);

struct SiPMCalibration{
    std::array<double,1> PheGeVS,PheGeVC;
    SiPMCalibration();
};

SiPMCalibration::SiPMCalibration(){
//    testdata
    PheGeVS[0] = 262;
    PheGeVC[0] = 44;
//
}

struct PMTCalibration{
    std::array<double,1> PheGeVPS,PheGeVPC;
    PMTCalibration();
};

PMTCalibration::PMTCalibration(){
//    testdata
    PheGeVPS[0] = 262;
    PheGeVPC[0] = 44;
}



class EventOut{
  public:
	EventOut(){};
	~EventOut(){};
	uint32_t EventID;

        float SPMT1, SPMT2, SPMT3, SPMT4, SPMT5, SPMT6, SPMT7, SPMT8;
    	float CPMT1, CPMT2, CPMT3, CPMT4, CPMT5, CPMT6, CPMT7, CPMT8;
        float SiPMPheC[160] = {0};
        float SiPMPheS[160] = {0};
        float VecTowerE[9] = {0};
	float EnergyTot,EscapedEnergy;
	float totSiPMCene = 0.;
	float totSiPMSene = 0.;
	int NSiPMZero= 0.;
	float SPMTenergy = 0.;
	float CPMTenergy = 0.;
	float XDWC1,XDWC2,YDWC1,YDWC2;
	int PShower, MCounter, C1, C2;
	int PShowerall;

	void CompSPMTene(){SPMTenergy = SPMT1+SPMT2+SPMT3+SPMT4+SPMT5+SPMT6+SPMT7+SPMT8;}
	void CompCPMTene(){CPMTenergy = CPMT1+CPMT2+CPMT3+CPMT4+CPMT5+CPMT6+CPMT7+CPMT8;}
        int SiPMCol(int index){ return index%16; }
        int SiPMRow(int index){ return index/16; }
        pair<double, double> SiPMSpos(int index){
           int row = index / 16;
           int column = index%16;
           double x = (column-7)*2-1.5;
           double y = 2.*sq3*(4-row)+sq3/2;
           return pair<double,double>(x,y);
        }
        pair<double, double> SiPMCpos(int index){
           int row = index / 16;
           int column = index%16;
           double x = (column-7)*2-0.5;
           double y = 2.*sq3*(4-row)+1.5*sq3;
           return pair<double,double>(x,y);
        }
};


class Event{

	public:
		//Constructor and de-constructor
		//
		Event(){};
		~Event(){};

		//Data members
		//
		int SPMT1, SPMT2, SPMT3, SPMT4, SPMT5, SPMT6, SPMT7, SPMT8;
		int CPMT1, CPMT2, CPMT3, CPMT4, CPMT5, CPMT6, CPMT7, CPMT8;
	        double	beamX, beamY;

		int SiPM_sci[160];
		int SiPM_che[160];

		void calibrate(const SiPMCalibration&, EventOut*);
		void calibratePMT(const PMTCalibration&, EventOut*);

};

void Event::calibrate(const SiPMCalibration& calibration, EventOut* evout){

	for(uint16_t i=0;i<160;i++){    
          evout->SiPMPheC[i] = SiPM_che[i]/calibration.PheGeVC[0];
	  evout->totSiPMCene +=SiPM_che[i]/calibration.PheGeVC[0];
          evout->SiPMPheS[i] = SiPM_sci[i]/calibration.PheGeVS[0];
	  evout->totSiPMSene +=SiPM_sci[i]/calibration.PheGeVS[0];
        }
	evout->NSiPMZero=0;
}

void Event::calibratePMT(const PMTCalibration& pmtcalibration, EventOut* evout){

    //PMT calibration
    //
    evout->SPMT1=SPMT1/pmtcalibration.PheGeVPS[0];
    evout->SPMT2=SPMT2/pmtcalibration.PheGeVPS[0];
    evout->SPMT3=SPMT3/pmtcalibration.PheGeVPS[0];
    evout->SPMT4=SPMT4/pmtcalibration.PheGeVPS[0];
    evout->SPMT5=SPMT5/pmtcalibration.PheGeVPS[0];
    evout->SPMT6=SPMT6/pmtcalibration.PheGeVPS[0];
    evout->SPMT7=SPMT7/pmtcalibration.PheGeVPS[0];
    evout->SPMT8=SPMT8/pmtcalibration.PheGeVPS[0];
    evout->CPMT1=CPMT1/pmtcalibration.PheGeVPC[0];
    evout->CPMT2=CPMT2/pmtcalibration.PheGeVPC[0];
    evout->CPMT3=CPMT3/pmtcalibration.PheGeVPC[0];
    evout->CPMT4=CPMT4/pmtcalibration.PheGeVPC[0];
    evout->CPMT5=CPMT5/pmtcalibration.PheGeVPC[0];
    evout->CPMT6=CPMT6/pmtcalibration.PheGeVPC[0];
    evout->CPMT7=CPMT7/pmtcalibration.PheGeVPC[0];
    evout->CPMT8=CPMT8/pmtcalibration.PheGeVPC[0];
}

ClassImp(EventOut)

void DRTB23SimConverter(const string run){

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
