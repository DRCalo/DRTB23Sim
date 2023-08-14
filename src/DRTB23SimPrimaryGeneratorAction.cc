//**************************************************
// \file DRTB23SimPrimaryGeneratorAction.cc
// \brief: Implementation of 
//         DRTB23SimPrimaryGeneratorAction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimPrimaryGeneratorAction.hh"
#include "DRTB23SimEventAction.hh"

//Includers from Geant4
//
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"

//Constructor
//
DRTB23SimPrimaryGeneratorAction::DRTB23SimPrimaryGeneratorAction(DRTB23SimEventAction* evtAction)
 : G4VUserPrimaryGeneratorAction(),
   fGeneralParticleSource( nullptr ),
   fEventAction(evtAction)
   /*fParticleGun( nullptr )*/ {

    //G4int nofParticles = 1;                           //for particle gun
    //fParticleGun = new G4ParticleGun(nofParticles);   //for particle gun
    fGeneralParticleSource = new G4GeneralParticleSource();

    //default G4GeneralParticleSource parameters (can be changed via UI)
    //
    G4ParticleDefinition* particleDefinition =
        G4ParticleTable::GetParticleTable()->FindParticle("e-");

    fGeneralParticleSource->SetParticleDefinition(particleDefinition);
    fGeneralParticleSource->SetParticlePosition( G4ThreeVector( 0.,0.,0. ) );

}

//De-constructor
//
DRTB23SimPrimaryGeneratorAction::~DRTB23SimPrimaryGeneratorAction() {
  
    delete fGeneralParticleSource;
    //delete fParticleGun;

}

//GeneratePrimaries() method
//
void DRTB23SimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    
    fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
    //fParticleGun->GeneratePrimaryVertex(anEvent);

    G4cout<<fGeneralParticleSource->GetParticleEnergy()<<G4endl;
    G4cout<<fGeneralParticleSource->GetParticlePosition()<<G4endl;
    G4cout<<fGeneralParticleSource->GetParticleDefinition()->GetPDGEncoding()<<G4endl;

    //Save primary particle energy, PDGID and x-y position
    //
    fEventAction->SavePrimaryEnergy(fGeneralParticleSource->GetParticleEnergy());
    fEventAction->SavePrimaryPDGID(fGeneralParticleSource->GetParticleDefinition()->GetPDGEncoding());
    fEventAction->SavePrimaryXY(fGeneralParticleSource->GetParticlePosition().x(),
                                fGeneralParticleSource->GetParticlePosition().y());

}

//**************************************************
