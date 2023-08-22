//**************************************************
// \file DRTB23SimOpticalPhysics.cc 
// \brief: Implementation of DRTB23SimOpticalPhysics class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

// Includers from Geant4
//
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"

// Includers from project files
//
#include "DRTB23SimOpticalPhysics.hh"

// Constructor
//
DRTB23SimOpticalPhysics::DRTB23SimOpticalPhysics() 
    : G4VPhysicsConstructor("Optical"),
    theWLSProcess(nullptr),
    theCerenkovProcess(nullptr),
    theScintProcess(nullptr),
    theRayleighScattering(nullptr),
    theMieHGScatteringProcess(nullptr),
    theBoundaryProcess(nullptr) {}
    
// De-constructor
//
DRTB23SimOpticalPhysics::~DRTB23SimOpticalPhysics() {}


void DRTB23SimOpticalPhysics::ConstructParticle() {
    
    G4OpticalPhoton::OpticalPhotonDefinition();

}

void DRTB23SimOpticalPhysics::ConstructProcess() {
    
    G4cout<<"DRTB23SimOpticalPhysics:: Add Optical Physics Processes"<<G4endl;
   
    // Initialize optical processes
    //
    theWLSProcess             = new G4OpWLS();
    theWLSProcess->UseTimeProfile("delta");
    theWLSProcess->UseTimeProfile("exponential"); 
    
    theScintProcess             = new G4Scintillation();
    theScintProcess->SetScintillationYieldFactor(1.);
    //theScintProcess->SetTrackSecondariesFirst(true);
    theScintProcess->SetScintillationExcitationRatio(0.0);
    theScintProcess->SetTrackSecondariesFirst(true);
    // Use Birks Correction in the Scintillation process
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    theScintProcess->AddSaturation(emSaturation);
    
    theCerenkovProcess          = new G4Cerenkov();
    theCerenkovProcess->SetMaxNumPhotonsPerStep(1000.);
    //theCerenkovProcess->SetTrackSecondariesFirst(true);
    
    theRayleighScattering     = new G4OpRayleigh();
 
    theMieHGScatteringProcess = new G4OpMieHG();
    
    theBoundaryProcess          = new G4OpBoundaryProcess();
   
    // Get optical photon process manager 
    //
    G4ProcessManager* pManager =
        G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    // if does not exist rais a fatal exception
    //
    if (!pManager) {
        G4Exception("OpticalPhysics::ConstructProcess()","",
                    FatalException,"Optical Photon without a Process Manager");
    }
    // Add processes to optical photon with optical photon process manager
    //
    //pManager->AddDiscreteProcess(theRayleighScattering);
    //pManager->AddDiscreteProcess(theMieHGScatteringProcess);
    //pManager->AddDiscreteProcess(theWLSProcess);
    pManager->AddDiscreteProcess(theBoundaryProcess);
   
    // Loop on particles and apply scintillation and cherenkov
    // processes to any physical candidate
    //
    auto theParticleIterator = GetParticleIterator();
    theParticleIterator->reset();
    while ( (*theParticleIterator)() ){
       
        // Get particle
        //
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4String particleName = particle->GetParticleName();
    
        // Get particle process manager
        //
        pManager = particle->GetProcessManager();
        // If does not exist raise a fatal exception
        //
        if (!pManager) {
            G4Exception("OpticalPhysics::ConstructProcess()","",
                        FatalException,"No Process Manager for particle");
        }
        // Add Cherenkov process to each candidate
        //
        if(theCerenkovProcess->IsApplicable(*particle)){
            pManager->AddProcess(theCerenkovProcess);
            pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
        }
        // Add Scintillation process to each candidate
        // (only for debugging)
        //
        /*if(theScintProcess->IsApplicable(*particle)){
            pManager->AddProcess(theScintProcess);
            pManager->SetProcessOrderingToLast(theScintProcess,idxAtRest);
            pManager->SetProcessOrderingToLast(theScintProcess,idxPostStep);
        }*/
    }//end while
}

//**************************************************
