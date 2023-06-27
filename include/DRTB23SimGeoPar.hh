//Namespace for test-beam geometry constants
//

#include "G4SystemOfUnits.hh"

namespace { 

    // module constants
    static constexpr G4int NofmodulesX = 3; 
    static constexpr G4int NofmodulesY = 3;
    static constexpr G4int modflag[9]={3,2,1,5,0,4,8,7,6};
    static constexpr G4int NoModulesSiPM=1;
    static constexpr G4int SiPMMod[1]={0};
    static constexpr G4int NofFiberscolumn = 16;
    static constexpr G4int NofFibersrow = 20;
    static constexpr G4int NoModulesActive=9;
    static constexpr G4double moduleZ = (1000.)*mm;
    static constexpr G4bool irot=false;
    static constexpr G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;

    // preshower constants
    static constexpr G4double preshower_pos_x = 0.0*cm;
    static constexpr G4double preshower_pos_y = 0.0*cm;
    static constexpr G4double preshower_pos_z = -53.0*cm;

    // Comment on naming convention for housing variables:
    // length = direction along local z
    // height = direction along local y
    // width  = direction along local x
    // thickness = from context

    // platform constants
    static constexpr G4double platform_half_height = 25*mm;     // Height guessed for now
    
    // bar/feet constants
    static constexpr G4double bar_half_length = 40.0*cm/2;
    static constexpr G4double bar_half_width  = 4.5*cm/2;
    static constexpr G4double bar_half_height = 9.0*cm/2;
    static constexpr G4double bar_wall_thickness = 10.0*mm;
    static constexpr G4double bar_pos_from_front = 16.0*cm;

    // housing constants
    static constexpr G4double housing_half_length = 145.5*cm/2;
    static constexpr G4double housing_half_width  = 18.0*cm/2;
    static constexpr G4double housing_half_height = 15.0*cm/2;
    static constexpr G4double side_wall_thickness = 1.5*mm;
    static constexpr G4double top_wall_thickness  = 1.4*mm;
    static constexpr G4double bot_wall_thickness  = 10.0*mm;

    // support constants (longitudinal support below housing)
    static constexpr G4double support_half_length = 163.5*cm/2;
    static constexpr G4double support_half_height = 10.0*cm/2;
    static constexpr G4double support_wall_thickness = 7.0*mm;

    // plastic cover constants
    static constexpr G4double plastic_cover_full_length = 20.0*mm;
    static constexpr G4double cutout_thickness = 4.5*mm;

}
