//**************************************************
// \file DRTB23SimSignalHelper.cc
// \brief: Implementation of DRTB23SimSignalHelper
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 1 September 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimSignalHelper.hh"

//Includers from Geant4
#include "G4Poisson.hh"
#include "G4Tubs.hh"
#include "G4NavigationHistory.hh"

//Includers from C++
//
#include <random>

DRTB23SimSignalHelper* DRTB23SimSignalHelper::instance = 0;

//Define (private) constructor (singleton)
//
DRTB23SimSignalHelper::DRTB23SimSignalHelper(){}

//Define Instance() method
//
DRTB23SimSignalHelper* DRTB23SimSignalHelper::Instance(){
    if (instance==0){
        instance = new DRTB23SimSignalHelper;
    }
    return DRTB23SimSignalHelper::instance;
}

//Define ApplyBirks() method
//
G4double DRTB23SimSignalHelper::ApplyBirks( const G4double& de, const G4double& steplength ) {
		
    return (de/steplength) / ( 1+fk_B*(de/steplength) ) * steplength;

}

//Define SmearSSignal() method
//
G4int DRTB23SimSignalHelper::SmearSSignal( const G4double& satde ) {
		
    return G4Poisson(satde*9.5);
		
}

//Define SmearCSignal() method
//
G4int DRTB23SimSignalHelper::SmearCSignal( ){
		
    return G4Poisson(0.153);

}

//Define GetDistanceToSiPM() method
//
G4double DRTB23SimSignalHelper::GetDistanceToSiPM(const G4Step* step) {

    // Get the pre-step point
    const G4StepPoint* preStepPoint = step->GetPreStepPoint();
    // Get the global position of the pre-step point
    G4ThreeVector globalPos = preStepPoint->GetPosition();
    // Get the local position of the pre-step point in the current volume's coordinate system
    G4ThreeVector localPos = preStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(globalPos);
    // G4cout << "Local Position (X,Y,Z): (" << localPos.x()/CLHEP::mm << ", " << localPos.y()/CLHEP::mm << ", " << localPos.z()/CLHEP::mm << ") mm" << G4endl;

    // Get the logical volume of the current step
    G4LogicalVolume* currentVolume = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    // Get the solid associated with the logical volume
    G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
    // Get the dimensions of the solid (size of the volume)
    G4double size = solid->GetZHalfLength();

    G4double distance_to_sipm = size - localPos.z();
    return distance_to_sipm;

}

//Define AttenuateHelper() method
G4int DRTB23SimSignalHelper::AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length) {
    double probability_of_survival = exp(-distance/attenuation_length);

    // Seed the random number generator
    std::random_device rd;
    std::default_random_engine rng(rd());

    // Define the distribution with the given probability x
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    G4int survived_photons = 0;
    for (int i=0; i<signal; i++)
    {
        // Simulate drawing between 0 and 1 with probability x of getting 1
        if (dist(rng) <= probability_of_survival) survived_photons++;
    }

    return survived_photons;

}

//Define AttenuateSSignal() method
//
G4int DRTB23SimSignalHelper::AttenuateSSignal(const G4int& signal, const G4double& distance) {

    return AttenuateHelper(signal, distance, fSAttenuationLength);    

}

//Define AttenuateCSignal() method
//
G4int DRTB23SimSignalHelper::AttenuateCSignal(const G4int& signal, const G4double& distance) {

    return AttenuateHelper(signal, distance, fCAttenuationLength);    
    
}

//**************************************************
