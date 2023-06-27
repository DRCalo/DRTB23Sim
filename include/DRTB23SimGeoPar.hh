//Namespace for test-beam geometry constants
//

#include "G4SystemOfUnits.hh"

namespace { 

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
}
