#include "G4SystemOfUnits.hh"
// TB2021 
/*
    const G4int NofmodulesX = 3; 
    const G4int NofmodulesY = 3;
    const G4int modflag[9]={3,2,1,5,0,4,8,7,6};
    const G4int NoModulesSiPM=1;
    const G4int SiPMMod[1]={0};
    const G4int NofFiberscolumn = 16;
    const G4int NofFibersrow = 20;
    const G4int NoModulesActive=9;
    const G4double moduleZ = (1000.)*mm;
    const G4bool irot=false;
*/

    const G4int NofmodulesX = 24; 
    const G4int NofmodulesY = 5;
    const G4int modflag[120]={-1,-1,-1,-1,-1,-1,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-1,-1,-1,-1,-1,-1,-1,
                              -1,-1,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,-1,-1, 
                              30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53, 
                              -1,-1,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,-1,-1, 
                              -1,-1,-1,-1,-1,-1,-1,74,75,76,77,78,79,80,81,82,83,-1,-1,-1,-1,-1,-1,-1}; 
    const G4int NoModulesSiPM=10;
    const G4int SiPMMod[10]={37,38,39,40,41,42,43,44,45,46};
    const G4int NofFiberscolumn = 64;
    const G4int NofFibersrow = 16;
    const G4int NoModulesActive=84;
    const G4double moduleZ = (2500.)*mm;
    const G4bool irot=true;

    const G4int NoFibersTower=NofFiberscolumn*NofFibersrow/2;
