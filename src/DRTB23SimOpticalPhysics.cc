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
DRTB23SimOpticalPhysics::DRTB23SimOpticalPhysics(const G4bool FullOptic, G4bool toggle) 
    : G4VPhysicsConstructor("Optical"),
      fFullOptic( FullOptic ) {
    
    // Initialize private members
    //
    theWLSProcess                = NULL;
    theScintProcess              = NULL;
    theCerenkovProcess           = NULL;
    theBoundaryProcess           = NULL;
    theAbsorptionProcess         = NULL;
    theRayleighScattering        = NULL;
    theMieHGScatteringProcess    = NULL;
    AbsorptionOn                 = toggle;

}

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
    //theWLSProcess             = new G4OpWLS();
    //theWLSProcess->UseTimeProfile("delta");
    //theWLSProcess->UseTimeProfile("exponential"); 
    
    theScintProcess             = new G4Scintillation();
    theScintProcess->SetScintillationYieldFactor(1.);
    //theScintProcess->SetTrackSecondariesFirst(true);
    //theScintProcess->SetScintillationYieldFactor(1.);
    theScintProcess->SetScintillationExcitationRatio(0.0);
    theScintProcess->SetTrackSecondariesFirst(true);
    // Use Birks Correction in the Scintillation process
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    theScintProcess->AddSaturation(emSaturation);
    
    theCerenkovProcess          = new G4Cerenkov();
    theCerenkovProcess->SetMaxNumPhotonsPerStep(1000.);
    //theCerenkovProcess->SetTrackSecondariesFirst(true);
    
    theAbsorptionProcess        = new G4OpAbsorption();
    
    //theRayleighScattering     = new G4OpRayleigh();
    
    //theMieHGScatteringProcess = new G4OpMieHG();
    
    theBoundaryProcess          = new G4OpBoundaryProcess();
   
    // Get optical photon process manager 
    //
    G4ProcessManager* pManager =
        G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    // if does not exist rais a fatal exception
    //
    if (!pManager) {
        G4cout<<"DRTB23Sim->Optical Photon without a Process Manager"<<G4endl;
        G4Exception("OpticalPhysics::ConstructProcess()","",
                    FatalException,"Optical Photon without a Process Manager");
    }
    // Add processes to optical photon with optical photon process manager
    //
    if (AbsorptionOn) pManager->AddDiscreteProcess(theAbsorptionProcess);
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
            G4cout<< "DRTB23Sim->Particle " <<
                particleName << "without a Process Manager";
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
        // Adding Scintillation process only if fFullOptic == true
        //
        if(theScintProcess->IsApplicable(*particle)){
            if (fFullOptic) {
                pManager->AddProcess(theScintProcess);
                pManager->SetProcessOrderingToLast(theScintProcess,idxAtRest);
                pManager->SetProcessOrderingToLast(theScintProcess,idxPostStep);
            }
            else {}
        } 
    }//end while
}

//**************************************************
