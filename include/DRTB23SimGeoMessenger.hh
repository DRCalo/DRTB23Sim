//**************************************************
// \file DRTB23SimGeoMessenger.hh
// \brief: Definition of DRTB23GeomMessenger class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 10 August 2023
//**************************************************

#ifndef DRTB23SimGeoMessenger_h
#define DRTB23SimGeoMessenger_h 1

//Includers from Geant4
//
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//Includers from project files
//
class DRTB23SimDetectorConstruction;

class DRTB23SimGeoMessenger final : public G4UImessenger {

    public:
        //Constructor and destructor
        DRTB23SimGeoMessenger(DRTB23SimDetectorConstruction *DetConstruction);
        ~DRTB23SimGeoMessenger();
        
        //Virtual methods from base class
        void SetNewValue(G4UIcommand *command, G4String newValue) override;

    private:
        //Members
        DRTB23SimDetectorConstruction *fDetConstruction;
        G4UIdirectory *fMsgrDirectory;
        G4UIcmdWithADoubleAndUnit *fXshiftcmd;
        G4UIcmdWithADoubleAndUnit *fYshiftcmd;
        G4UIcmdWithADoubleAndUnit *fOrzrotcmd;
        G4UIcmdWithADoubleAndUnit *fVerrotcmd;

};

#endif //DRTB23SimGeoMessenger_h

//**************************************************
