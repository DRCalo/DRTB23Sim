//**************************************************
// \file DRTB23SimOpticalPhysics.hh 
// \brief: Definition of DRTB23SimOpticalPhysics class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************//

// Prevent including header multiple times
//
#ifndef DRTB23SimOpticalPhysics_h
#define DRTB23SimOpticalPhysics_h 1

// Includers from Geant4
//
#include "globals.hh"
#include "G4OpWLS.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4VPhysicsConstructor.hh"

class DRTB23SimOpticalPhysics : public G4VPhysicsConstructor {
    
    public: 
        // Constructor
        //
        DRTB23SimOpticalPhysics();
        // Deconstructor
        //
        virtual ~DRTB23SimOpticalPhysics();
    
        virtual void ConstructParticle();
        virtual void ConstructProcess();
   
    private:
    
        G4OpWLS*             theWLSProcess;
        G4Cerenkov*          theCerenkovProcess;
        G4Scintillation*     theScintProcess;
        G4OpRayleigh*        theRayleighScattering;
        G4OpMieHG*           theMieHGScatteringProcess;
        G4OpBoundaryProcess* theBoundaryProcess;
 
};

#endif

//**************************************************//
