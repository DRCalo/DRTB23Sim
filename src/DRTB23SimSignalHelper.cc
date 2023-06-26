//**************************************************
// \file DRTB23SimSignalHelper.cc
// \brief: Implementation of DRTB23SimSignalHelper
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 1 September 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimSignalHelper.hh"

//Includers from Geant4
#include "G4Poisson.hh"

DRTB23SimSignalHelper* DRTB23SimSignalHelper::instance = 0;

//Define (private) constructor (singleton)
//
DRTB23SimSignalHelper::DRTB23SimSignalHelper(){}

//Define Instance() method
//
DRTB23SimSignalHelper* DRTB23SimSignalHelper::Instance(){
    if (instance==0){
        instance = new DRTB23SimSignalHelper;
    }
    return DRTB23SimSignalHelper::instance;
}

//Define ApplyBirks() method
//
G4double DRTB23SimSignalHelper::ApplyBirks( const G4double& de, const G4double& steplength ) {
		
    const G4double k_B = 0.126; //Birks constant
    return (de/steplength) / ( 1+k_B*(de/steplength) ) * steplength;

}

//Define SmearSSignal() method
//
G4int DRTB23SimSignalHelper::SmearSSignal( const G4double& satde ) {
		
    return G4Poisson(satde*9.5);
		
}

//Define SmearCSignal() method
//
G4int DRTB23SimSignalHelper::SmearCSignal( ){
		
    return G4Poisson(0.153);

}

//**************************************************
