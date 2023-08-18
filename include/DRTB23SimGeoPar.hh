//Namespace for test-beam geometry constants
//

#ifndef DRTB23SimGeoPar_h
#define DRTB23SimGeoPar_h

//Includers from Geant4
//
#include "G4SystemOfUnits.hh"

namespace { 

    // module constants
    constexpr G4int NofmodulesX = 3; 
    constexpr G4int NofmodulesY = 3;
    constexpr G4int modflag[9]={3,2,1,5,0,4,8,7,6};
    constexpr G4int NoModulesSiPM=1;
    constexpr G4int SiPMMod[1]={0};
    constexpr G4int NofFiberscolumn = 16;
    constexpr G4int NofFibersrow = 20;
    constexpr G4int NoModulesActive=9;
    constexpr G4double moduleZ = (1000.)*mm;
    constexpr G4bool irot=false;
    constexpr G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;

    // preshower constants
    constexpr G4double preshower_pos_x = 0.0*cm;
    constexpr G4double preshower_pos_y = 0.0*cm;
    constexpr G4double preshower_pos_z = -53.0*cm;

    // Comment on naming convention for housing variables:
    // length = direction along local z
    // height = direction along local y
    // width  = direction along local x
    // thickness = from context

    // platform constants
    constexpr G4double platform_radius = 1200*mm;    // Radius guessed for now
    constexpr G4double platform_half_height = 25*mm;     // Height guessed for now
    
    // bar/feet constants
    constexpr G4double bar_half_length = 40.0*cm/2;
    constexpr G4double bar_half_width  = 4.5*cm/2;
    constexpr G4double bar_half_height = 9.0*cm/2;
    constexpr G4double bar_wall_thickness = 10.0*mm;
    constexpr G4double bar_pos_from_front = 16.0*cm;

    // housing constants
    constexpr G4double housing_half_length = 145.5*cm/2;
    constexpr G4double housing_half_width  = 18.0*cm/2;
    constexpr G4double housing_half_height = 15.0*cm/2;
    constexpr G4double side_wall_thickness = 1.5*mm;
    constexpr G4double top_wall_thickness  = 1.4*mm;
    constexpr G4double bot_wall_thickness  = 10.0*mm;

    // support constants (longitudinal support below housing)
    constexpr G4double support_half_length = 163.5*cm/2;
    constexpr G4double support_half_height = 10.0*cm/2;
    constexpr G4double support_wall_thickness = 7.0*mm;

    // plastic cover constants
    constexpr G4double plastic_cover_full_length = 20.0*mm;
    constexpr G4double cutout_thickness = 4.5*mm;

}

#endif // DRTB23SimGeoPar_h
