//**************************************************
// \file PhysicsEvent.h
// \brief: implementation of Event class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
// 			@lopezzot
//          Edoardo Proserpio (Uni Insubria)
//          Iacopo Vivarelli (Uni Sussex)
// \start date: 20 August 2021
//**************************************************

#include <iostream>
#include <array>
#include <stdint.h>
#include <string>
#include <fstream>

#ifndef Event_H
#define Event_H

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

#endif

//**************************************************
