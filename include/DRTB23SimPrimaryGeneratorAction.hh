//**************************************************
// \file DRTB23SimPrimaryGeneratorAction.hh
// \brief: Definition of 
//         DRTB23SimPrimaryGeneratorAction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including headers multiple times
//
#ifndef DRTB23SimPrimaryGeneratorAction_h
#define DRTB23SimPrimaryGeneratorAction_h 1

//Includers from Geant4
//
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

//Forwards declarations from Geant4
//
class G4GeneralParticleSource;
//class G4ParticleGun; //in case user want to switch to G4ParticleGun
class G4Event;

//Forward declarations from project
//
class DRTB23SimEventAction;

class DRTB23SimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 
    public:
        //Constructor()
        //
        DRTB23SimPrimaryGeneratorAction(DRTB23SimEventAction* evtAction);    

        //De-constructor()
        //
        virtual ~DRTB23SimPrimaryGeneratorAction();

        virtual void GeneratePrimaries(G4Event* event);
  
    private:
        G4GeneralParticleSource* fGeneralParticleSource;
        //G4ParticleGun*  fParticleGun; // G4ParticleGun
        DRTB23SimEventAction* fEventAction;

};

#endif

//**************************************************
