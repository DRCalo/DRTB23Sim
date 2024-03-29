//**************************************************
// \file DRTB23SimActionInitialization.cc
// \brief: Implementation of DRTB23SimActionInitialization class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimActionInitialization.hh"
#include "DRTB23SimPrimaryGeneratorAction.hh"
#include "DRTB23SimRunAction.hh"
#include "DRTB23SimEventAction.hh"
#include "DRTB23SimSteppingAction.hh"

//Constructor
//
DRTB23SimActionInitialization::DRTB23SimActionInitialization()
    : G4VUserActionInitialization()
{}

//De-constructor
//
DRTB23SimActionInitialization::~DRTB23SimActionInitialization() {}

//BuildForMaster() method
//
void DRTB23SimActionInitialization::BuildForMaster() const {
    
    auto eventAction = new DRTB23SimEventAction;
    SetUserAction( new DRTB23SimRunAction( eventAction ) );

}

//Build() method
//
void DRTB23SimActionInitialization::Build() const {
  
    auto eventAction = new DRTB23SimEventAction;
    SetUserAction(new DRTB23SimRunAction( eventAction ));
    SetUserAction(new DRTB23SimPrimaryGeneratorAction(eventAction) );
    SetUserAction(eventAction);
    SetUserAction(new DRTB23SimSteppingAction(eventAction));

}  

//**************************************************
