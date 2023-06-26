//**************************************************
// \file DRTB23SimRunAction.hh 
// \brief: Definition of DRTB23SimRunAction class 
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef DRTB23SimRunAction_h
#define DRTB23SimRunAction_h 1

//Includers from Geant4
//
#include "G4UserRunAction.hh"
#include "globals.hh"

class DRTB23SimEventAction;
class G4Run;

class DRTB23SimRunAction : public G4UserRunAction {
    
    public:
        //Constructor
        //
        DRTB23SimRunAction( DRTB23SimEventAction* eventAction );
        //De-constructor
        //
        virtual ~DRTB23SimRunAction();

        //Methods
        //
        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);

    private:
        DRTB23SimEventAction* fEventAction;

};

#endif

//**************************************************
