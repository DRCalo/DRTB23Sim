//**************************************************
// \file DRTB23SimSignalHelper.hh
// \brief: Definition of DRTB23SimSignalHelper class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 1 September 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef DRTB23SimSignalHelper_h
#define DRTB23SimSignalHelper_h

//Includers from Geant4
//
#include "globals.hh"
#include "G4Step.hh"

class DRTB23SimSignalHelper {

    private:

        static DRTB23SimSignalHelper* instance;

	//Private constructor (singleton)
        //
	DRTB23SimSignalHelper();

    public:

    	static DRTB23SimSignalHelper* Instance();

    	G4double ApplyBirks( const G4double& de, const G4double& steplength );

	    G4int SmearSSignal( const G4double& de );

    	G4int SmearCSignal( );

		G4double GetDistanceToSiPM(const G4Step* step);

		G4int AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length);

		G4int AttenuateSSignal(const G4int& signal, const G4double& distance);

		G4int AttenuateCSignal(const G4int& signal, const G4double& distance);
};

#endif

//**************************************************
