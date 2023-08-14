//**************************************************
// \file DRTB23SimPhysicsList.hh
// \brief: Definition of DRTB23SimPhysicsList class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef DRTB23SimPhysicsList_h
#define DRTB23SimPhysicsList_h 1

//Includers from Geant4
//
#include "G4VModularPhysicsList.hh"

//Includers from project files
//
#include "DRTB23SimOpticalPhysics.hh"

class DRTB23SimPhysicsList : public G4VModularPhysicsList{
    
    public:
        //Constructor
        //
        DRTB23SimPhysicsList(G4String);

        //De-constructor
        //
        virtual ~DRTB23SimPhysicsList();
    
        DRTB23SimOpticalPhysics* OpPhysics;
    
        G4bool AbsorptionOn;
};

#endif

//**************************************************


