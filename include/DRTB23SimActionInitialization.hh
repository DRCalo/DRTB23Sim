//**************************************************
// \file DRTB23SimActionInitialization.hh
// \brief: Definition of DRTB23SimActionInitialization class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef DRTB23SimActionInitialization_h
#define DRTB23SimActionInitialization_h 1

//Includers from Geant4
//
#include "G4VUserActionInitialization.hh"
#include "G4Types.hh"

//Includers from C++
//
#include <chrono>
#include <random>

//Forward declaration
//
class DRTB23SimDetectorConstruction;

class DRTB23SimActionInitialization : public G4VUserActionInitialization {
    
    public:
        //Constructor
        //
        DRTB23SimActionInitialization(DRTB23SimDetectorConstruction* );
        virtual ~DRTB23SimActionInitialization();

        virtual void BuildForMaster() const;
        virtual void Build() const;

    private:

	DRTB23SimDetectorConstruction* fDetConstruction;

};

#endif

//**************************************************
