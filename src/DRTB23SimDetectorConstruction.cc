//**************************************************
// \file DRTB23SimDetectorConstruction.cc
// \brief: Implementation of 
//         DRTB23SimDetectorConstruction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DRTB23SimDetectorConstruction.hh"
#include "DRTB23SimGeoMessenger.hh"

//Includers from Geant4
//
#include <random>
#include <iostream>
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Sphere.hh"
#include "G4Colour.hh"
#include "G4TwoVector.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

//Messenger constructor
//
DRTB23SimGeoMessenger::DRTB23SimGeoMessenger(DRTB23SimDetectorConstruction* DetConstruction)
    : fDetConstruction(DetConstruction){
    
    //The Messenger directory and commands must be initialized
    //within constructor
    //
    fMsgrDirectory = new G4UIdirectory("/tbgeo/");
    fMsgrDirectory->SetGuidance("Set movable parameters in test-beam geometry.");

    fXshiftcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/xshift", this);
    fXshiftcmd->SetGuidance("Shift test-beam platform x direction (default unit mm)");
    fXshiftcmd->SetDefaultUnit("mm");
    fYshiftcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/yshift", this);
    fYshiftcmd->SetGuidance("Shift test-beam platform y direction (default unit mm)");
    fYshiftcmd->SetDefaultUnit("mm");
    fOrzrotcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/horizrot", this);
    fOrzrotcmd->SetGuidance("Rotate platform (default deg)");
    fOrzrotcmd->SetDefaultUnit("deg");
    fVerrotcmd = new G4UIcmdWithADoubleAndUnit("/tbgeo/vertrot", this);
    fVerrotcmd->SetGuidance("Lift up calorimeter from back side (default deg)");
    fVerrotcmd->SetDefaultUnit("deg");
}

//Messenger destructor
//
DRTB23SimGeoMessenger::~DRTB23SimGeoMessenger(){

    //The Messenger fields should be deleted
    //in destructor
    delete fMsgrDirectory;
    delete fXshiftcmd;
    delete fYshiftcmd;
    delete fOrzrotcmd;
    delete fVerrotcmd;
}

//Messenger SetNewValue virtual method from base class
//
void DRTB23SimGeoMessenger::SetNewValue(G4UIcommand* command, G4String newValue){

    if(command == fXshiftcmd){
        fDetConstruction->SetXshift(fXshiftcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: x-shifted test-beam setup by "<<fDetConstruction->GetXshift()<<" mm"<<G4endl;
    }
    else if(command == fYshiftcmd){
        fDetConstruction->SetYshift(fYshiftcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: y-shifted test-beam setup by "<<fDetConstruction->GetYshift()<<" mm"<<G4endl;
    }
    else if(command == fOrzrotcmd){
        fDetConstruction->SetOrzrot(fOrzrotcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: orz-rotated test-beam setup by "<<fDetConstruction->GetOrzrot()<<" rad"<<G4endl;
    }
    else if(command == fVerrotcmd){
        fDetConstruction->SetVerrot(fVerrotcmd->GetNewDoubleValue(newValue));
        G4cout<<"tbgeo: ver-rotated test-beam setup by "<<fDetConstruction->GetVerrot()<<" rad"<<G4endl;
    }
}

//  sqrt3 constants used in code
//  reciprocal of sqrt3 given as number that divided 
//  by sqrt 3 gives 3 
//
const G4double sq3=1.733;
const G4double sq3m1=sq3/3.;

//Constructor
//
DRTB23SimDetectorConstruction::DRTB23SimDetectorConstruction(const G4bool VertRot)
    : G4VUserDetectorConstruction(),
    fCheckOverlaps(false),
    fLeakCntPV(nullptr),
    fWorldPV(nullptr),
    fVertRot(VertRot){

    fGeoMessenger = new DRTB23SimGeoMessenger(this);
}

//De-constructor
//
DRTB23SimDetectorConstruction::~DRTB23SimDetectorConstruction() {}

//Define Construct() method
G4VPhysicalVolume* DRTB23SimDetectorConstruction::Construct() {
  
    // Define volumes
    return DefineVolumes();
}

G4VPhysicalVolume* DRTB23SimDetectorConstruction::DefineVolumes() {

    //--------------------------------------------------
    //Define Elements, Mixtures and Materials
    //--------------------------------------------------

    auto nistManager = G4NistManager::Instance();

    //Elements
    //
    G4String name, symbol;    
    G4double a, z;            // a=mass of a mole, z=mean number of protons;  
  
    a = 1.01*g/mole;
    G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a); //Hidrogen

    a = 12.01*g/mole;
    G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a); //Carbon

    a = 16.00*g/mole;
    G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a); //Oxygen

    a = 28.09*g/mole;
    G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a); //Silicon
  
    a = 18.9984*g/mole;
    G4Element* elF  = new G4Element("Fluorine",symbol="F" , z= 9., a); //Fluorine

    //a = 63.546*g/mole;
    //G4Element* elCu = new G4Element("Copper", symbol="Cu", z=29., a); //Copper
    auto elCu = nistManager->FindOrBuildElement(29, true);

    //a = 65.38*g/mole;
    //G4Element* elZn = new G4Element("Zinc", symbol="Zn", z=30., a); //Zinc
    auto elZn = nistManager->FindOrBuildElement(30, true);

    //Materials 
    //

    // Polystyrene from elements (C5H5)
    G4Material* Polystyrene = new G4Material("Polystyrene", 1.05*g/cm3, 2);
    Polystyrene->AddElement(elC, 8);
    Polystyrene->AddElement(elH, 8); 

    // PMMA material from elements (C502H8)
    // 
    auto PMMA = new G4Material("PMMA", 1.19*g/cm3, 3); 
    PMMA->AddElement(elC, 5);
    PMMA->AddElement(elO, 2);
    PMMA->AddElement(elH, 8); 
    
    // Fluorinated Polymer material from elements (C2F2)
    // material for the cladding of the Cherenkov fibers
    auto fluorinatedPolymer = new G4Material("Fluorinated_Polymer", 1.43*g/cm3, 2);
    fluorinatedPolymer->AddElement(elC,2);
    fluorinatedPolymer->AddElement(elF,2);

    // Glass material from elements (SiO2)
    //
    auto Glass = new G4Material("Glass", 2.4*g/cm3, 2);
    Glass -> AddElement(elSi, 1);
    Glass -> AddElement(elO, 2); 

    // Mixtures
    //

    // CuZn37 (Brass)
    //
    const double BrassDensity = 8.44*g/cm3;
    auto CuZn37 = new G4Material(name="Brass", BrassDensity, 2);
    CuZn37->AddElement(elCu, 0.7);
    CuZn37->AddElement(elZn, 0.3);

    // Assign material to the calorimeter volumes
    //
    G4Material* defaultMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
    //G4Material* absorberMaterial = nistManager->FindOrBuildMaterial("G4_Cu");
    G4Material* SiMaterial = nistManager->FindOrBuildMaterial("G4_Si");
    G4Material* LeadMaterial = nistManager->FindOrBuildMaterial("G4_Pb");
    G4Material* PSScinMaterial = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");
    G4Material* absorberMaterial = G4Material::GetMaterial("Brass");
    G4Material* ScinMaterial = G4Material::GetMaterial("Polystyrene");
    G4Material* CherMaterial = G4Material::GetMaterial("PMMA");
    G4Material* GlassMaterial = G4Material::GetMaterial("Glass");
    G4Material* CladCherMaterial = G4Material::GetMaterial("Fluorinated_Polymer");
    G4Material* aluminiumMaterial = nistManager->FindOrBuildMaterial("G4_Al");
    G4Material* PVCMaterial = nistManager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

    //--------------------------------------------------
    //Define Optical Properties
    //--------------------------------------------------

    // Use Energy(eV)=1.24/waevelenght(um)
    // 2.034eV is 610nm RED 
    // 2.75eV is 450nm BLUE (peak of scintillating fibers)
    // 3.09eV is 400nm VIOLET (end of visible)
    //4.1eV is 300nm UV (cherenkov peak is 310-350nm)
    //
    const G4int ENTRIES = 32;
    G4double photonEnergy[ENTRIES] =                    
        { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,   
          2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     
          2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
          2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
          2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV, 
          3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV, 
          3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
          3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV }; 
    G4double rindexScin[ENTRIES] =
        { 1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59,
          1.59, 1.59, 1.59, 1.59 };
    /*G4double absorptionScin[ENTRIES] =
        { 400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm,
          400*cm, 400*cm, 400*cm, 400*cm };*/

    G4MaterialPropertiesTable *MPTScin = new G4MaterialPropertiesTable();
    MPTScin -> AddProperty("RINDEX", 
        photonEnergy, rindexScin, ENTRIES)->SetSpline(true);
    /*MPTScin -> AddProperty("ABSLENGTH",
         photonEnergy, absorptionScin, ENTRIES)->SetSpline(true);*/

    G4double rindexCher[ENTRIES] =
        { 1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49,
          1.49, 1.49, 1.49, 1.49 };
    /*G4double absorptionCher[ENTRIES] = 
        { 890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm,
          890*cm, 890*cm, 890*cm, 890*cm };*/

    G4MaterialPropertiesTable *MPTCher = new G4MaterialPropertiesTable();
    MPTCher -> AddProperty("RINDEX",
            photonEnergy, rindexCher, ENTRIES)->SetSpline(true);
    /*MPTCher -> AddProperty("ABSLENGTH", 
            photonEnergy, absorptionCher, ENTRIES)->SetSpline(true);*/
    CherMaterial -> SetMaterialPropertiesTable(MPTCher);

    G4double rindexCherclad[ENTRIES] =
        { 1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42,
          1.42, 1.42, 1.42, 1.42 };

    G4MaterialPropertiesTable *MPTCherclad = new G4MaterialPropertiesTable();
    MPTCherclad -> AddProperty("RINDEX", 
        photonEnergy, rindexCherclad, ENTRIES)->SetSpline(true);
    CladCherMaterial -> SetMaterialPropertiesTable(MPTCherclad);

    G4double rindexglass[ENTRIES] =
        { 1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51,
          1.51, 1.51, 1.51, 1.51 };

    G4MaterialPropertiesTable *MPTglass = new G4MaterialPropertiesTable();
    MPTglass -> AddProperty("RINDEX", 
            photonEnergy, rindexglass, ENTRIES)->SetSpline(true);
    GlassMaterial -> SetMaterialPropertiesTable(MPTglass);

    G4double rindexSi[ENTRIES] =
        { 3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42,
          3.42, 3.42, 3.42, 3.42 };

    G4double absorptionSi[ENTRIES] = 
        { 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
          0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm };

    G4MaterialPropertiesTable *MPTSi = new G4MaterialPropertiesTable();
    MPTSi -> AddProperty("RINDEX", photonEnergy, rindexSi, ENTRIES)->SetSpline(true);
    MPTSi -> AddProperty("ABSLENGHT", 
        photonEnergy, absorptionSi, ENTRIES)->SetSpline(true);
    SiMaterial -> SetMaterialPropertiesTable(MPTSi); 
  
    // Scintillating proprieties of the scintillating fiber material
    // Birks constant of the polystyrene
    //
    G4double Scin_FAST[ENTRIES] = // Emission spectrum for the fast component 
        { 0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.1,
          0.2, 0.4, 0.6, 0.8,
          1., 0.8, 0.6, 0.1,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0. };
    G4double Scin_SLOW[ENTRIES] = // Emission spectrum for the slow component
        { 0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0. };

    ScinMaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

    MPTScin -> AddProperty("FASTCOMPONENT", photonEnergy, Scin_FAST, ENTRIES);
    MPTScin -> AddProperty("SLOWCOMPONENT", photonEnergy, Scin_SLOW, ENTRIES);
    MPTScin -> AddConstProperty("SCINTILLATIONYIELD", 10000./MeV); 
    // Typical is 10000./MeV (this is what makes full simulations long as hell)
    MPTScin -> AddConstProperty("RESOLUTIONSCALE", 1.0); 
    // Broad the fluctuation of photons produced
    MPTScin -> AddConstProperty("FASTTIMECONSTANT", 2.8*ns);
    MPTScin -> AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
    MPTScin -> AddConstProperty("YIELDRATIO", 1.0); 
    // I don't want a slow component, if you want it must change
    ScinMaterial -> SetMaterialPropertiesTable(MPTScin);

    //--------------------------------------------------
    //Define Volumes
    //--------------------------------------------------
    
//  tubes are packed in Y, so Y distance between two tubes is sqrt(3)/2*2.*tuberadius
    double tolerance = 0.0*mm;
    G4double tuberadius = 1.0*mm;
    G4double dtubeY=sq3*tuberadius;
    G4double dtubeX=2.*tuberadius;
    G4double moduleX = (2*NofFiberscolumn+1)*tuberadius; 
    G4double moduleY = (NofFibersrow-1.)*dtubeY+4.*tuberadius*sq3m1;
//  side of hexagon in which tube is inscribed
//    G4double side=tuberadius*2*sq3m1;

    // Geometry parameters of the world, world is a G4Box
    //
    G4double worldX = 200 * moduleX;
    G4double worldY = 200 * moduleY;
    G4double worldZ = 60 * moduleZ;

    // Geometry parameters of the fiber
    //
    //G4double fiberradius = 0.5*mm;
    G4double fiberZ = moduleZ;
    
    // Geometry parameters of the core
    //
    G4double coreradius = 0.485*mm;
    G4double coreZ = moduleZ;

    // Geometry parameters of the cladding
    //
    G4double claddingradiusmin = 0.485*mm;
    G4double claddingradiusmax = 0.50*mm;
    G4double claddingZ = moduleZ;
    
    // Geometry parameters of the tube
    //
    //G4double tubeZ = fiberZ;

    // Geometry parameters of the module equipped
    //
    G4double moduleequippedZ = moduleZ;
    G4double moduleequippedX = moduleX; 
    G4double moduleequippedY = moduleY;

    //Preshower dimensions
    //
    G4double PSX = 9.2*cm;
    G4double PSY = 9.2*cm;
    G4double PSZ = 1.*cm;

    // Building geometries
    //
    // World
    //
    G4VSolid* worldS  = new G4Box("World", worldX/2, worldY/2, worldZ/2); 
                         
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,          
                                                   defaultMaterial, 
                                                   "World");       
  
    worldLV->SetVisAttributes(G4VisAttributes::Invisible);

    fWorldPV = new G4PVPlacement( 0,                // no rotation
                                  G4ThreeVector(),  // at (0,0,0)
                                  worldLV,          // its logical
                                  "World",          // its name
                                  0,                // its mother
                                  false,            // no boolean oper 
                                  0,                // copy number
                                  fCheckOverlaps);  // check overlaps 

    /*************************************************************************************************
    * Stage Volume:                                                                                  *
    * "Virtual" Volume only used to move platform, calorimeter and preshower in unison wrt. the beam *
    **************************************************************************************************/
    G4Box*            stage_solid   = new G4Box("stage_solid", platform_radius, platform_radius, platform_radius);
    G4LogicalVolume*  stage_logical = new G4LogicalVolume(stage_solid, defaultMaterial, "stage_logical");
    stage_logical->SetVisAttributes(G4VisAttributes::Invisible);
    G4double stage_x = fXshift; //set via G4UIMessenger, default=0.
    G4double stage_y = fYshift; //set via G4UIMessenger, default=0.
    /*G4VPhysicalVolume* stage_physical =*/ new G4PVPlacement(0,                // no rotation
                                                              G4ThreeVector(stage_x,stage_y, 0),
                                                              stage_logical,    // its logical
                                                              "Stage",          // its name
                                                              worldLV,          // its mother
                                                              false,            // no boolean oper 
                                                              0,                // copy number
                                                              fCheckOverlaps);  // check overlaps 



    //Preshower
    //
    auto PSSolid = new G4Box("Preshower", PSX/2., PSY/2., PSZ/2.);

    auto PSLV = new G4LogicalVolume(PSSolid, defaultMaterial, "Preshower");

    new G4PVPlacement( 0, 
		       G4ThreeVector(preshower_pos_x, preshower_pos_y, preshower_pos_z - PSZ/2.),
		       PSLV,
		       "Preshower",
		       stage_logical,
		       false,
		       0,
		       fCheckOverlaps);	 

    auto PSLeadSolid = new G4Box("Preshower_pb", PSX/2., PSY/2., PSZ/4.);

    auto PSLeadLV = new G4LogicalVolume(PSLeadSolid, LeadMaterial, "Preshower_pb");

    new G4PVPlacement( 0, 
		       G4ThreeVector(0.,0.,-PSZ/4.),
		       PSLeadLV,
		       "Preshower_pb",
		       PSLV,
		       false,
		       0,
		       fCheckOverlaps);	 

    G4VisAttributes* PbVisAtt = new G4VisAttributes( G4Colour::Grey() );
    PbVisAtt->SetVisibility(true);
    PbVisAtt->SetForceSolid(true);
    PSLeadLV->SetVisAttributes( PbVisAtt );

    auto PSScinSolid = new G4Box("Preshower_scin", PSX/2., PSY/2., PSZ/4.);

    auto PSScinLV = new G4LogicalVolume(PSScinSolid, PSScinMaterial, "Preshower_scin");

    new G4PVPlacement( 0, 
		       G4ThreeVector(0.,0.,PSZ/4.),
                       PSScinLV,
	               "Preshower_scin",
                       PSLV,
                       false,	
                       0,
                       fCheckOverlaps);	 

    G4VisAttributes* PSScinVisAtt = new G4VisAttributes( G4Colour::Cyan() );
    PSScinVisAtt->SetVisibility(true);
    PSScinLV->SetVisAttributes( PSScinVisAtt );
 
    //Absorber to calculate leakage
    //
    G4VSolid* leakageabsorber = new G4Sphere("leakageabsorber",                        
        7.*m, 7.1*m, 0.*deg, 360.*deg, 0.*deg, 180.*deg); 
    
    G4LogicalVolume* leakageabsorberLV = new G4LogicalVolume(leakageabsorber,
                                                             defaultMaterial,  
                                                             "leakageabsorber");        
    
    leakageabsorberLV->SetVisAttributes(G4VisAttributes::Invisible);   

    fLeakCntPV = new G4PVPlacement( 0, G4ThreeVector(),
				    leakageabsorberLV,         
                                    "leakageabsorber",
                                    worldLV,               
                                    false,          
                                    0,               
                                    fCheckOverlaps);

   // Module equipped
   //
   // Basic module structure: extrusion of an hexcell shape
    G4TwoVector offA(0,0), offB(0,0);
    G4double scaleA = 1, scaleB = 1;
    auto polygon=calcmod(tuberadius, NofFibersrow, NofFiberscolumn);
    G4VSolid* moduleequippedS = new G4ExtrudedSolid("moduleequipped", 
		                                     polygon, 
						     moduleequippedZ/2, 
						     offA, scaleA, offB, scaleB);
    G4LogicalVolume* moduleequippedLV = new G4LogicalVolume(moduleequippedS,
                                                            defaultMaterial,
                                                            "moduleequipped"); 

    // Calorimeter (matrix of modules equipped)
    // 
    G4double caloX=0.;
    G4double caloY=0.;
    G4double caloZ=0.;
    if(irot) {
      caloX=moduleequippedY*NofmodulesX/2;
      caloY=moduleequippedX*NofmodulesY/2;
      caloZ=moduleequippedZ/2;      
    }
    else {
      caloX=moduleequippedX*NofmodulesX/2;
      caloY=moduleequippedY*NofmodulesY/2;
      caloZ=moduleequippedZ/2;      
    }
    G4VSolid* CalorimeterS = new G4Box("CalorimeterS",caloX,caloY,caloZ);

    G4LogicalVolume* CalorimeterLV = new G4LogicalVolume( CalorimeterS,
                                                          defaultMaterial,
                                                          "CalorimeterLV");

    // Modules equipped placement
    //
    G4int copynumbermodule = 0;
    G4double m_x, m_y;
    G4ThreeVector vec_m;
    for(int column=0; column<NofmodulesX; column++){ 
        for(int row=0; row<NofmodulesY; row++){
	    G4int ii=column+row*NofmodulesX;
            if(irot) {
              m_x = -dtubeY*NofFibersrow*column+ dtubeY*NofFibersrow*(NofmodulesX-1)/2;
              m_y = row*NofFiberscolumn*dtubeX-NofFiberscolumn*dtubeX*(NofmodulesY-1)/2;
            }
            else {
              m_x = -dtubeX*NofFiberscolumn*column+ dtubeX*NofFiberscolumn*(NofmodulesX-1)/2;
              m_y = row*NofFibersrow*dtubeY-NofFibersrow*dtubeY*(NofmodulesY-1)/2;
            }	     	    
            if(modflag[ii]>=0) {        
              copynumbermodule = modflag[ii];
//              copynumbermodule = (1+column)+(row*NofmodulesX);
//	      std::cout << " column " << column << " row " << row << " cpnm " << copynumbermodule << std::endl;
// setup for 90 deg rotation of modules
              if(irot){
                G4RotationMatrix rotmod  = G4RotationMatrix();
                rotmod.rotateZ(90.*degree);
                G4ThreeVector posm;
                posm.setX(m_x);
                posm.setY(m_y);
                posm.setZ(0.);
                G4Transform3D transfm = G4Transform3D(rotmod,posm);
                new G4PVPlacement(transfm,
                                  moduleequippedLV,     
                                  "moduleequipped",                        
                                  CalorimeterLV,                      
                                  false,                          
                                  copynumbermodule); 
	      }
	      else {
// simple positioning
                vec_m.setX(m_x);
                vec_m.setY(m_y);
                vec_m.setZ(0.);
                new G4PVPlacement(0,
                                  vec_m,              
                                  moduleequippedLV,     
                                  "moduleequipped",                        
                                  CalorimeterLV,                      
                                  false,                          
                                  copynumbermodule); 
              } // irot
	    } //ifmod
        };
    }; 
 
    // Calorimeter placement (with rotation wrt beam axis)
    //
    G4RotationMatrix rotm  = G4RotationMatrix();
    G4double xrot=2.5*deg;
    G4double yrot=2.5*deg;
    rotm.rotateX(xrot);  
    rotm.rotateY(yrot);
    G4ThreeVector position;
    G4double ycomp=-1090*mm*sin(xrot);
    G4double xcomp=1090*mm*sin(yrot);
    position.setX(xcomp);
    position.setY(ycomp);
    position.setZ(0.);
    G4Transform3D transform = G4Transform3D(rotm,position); 

    /***********************************************************
    * Volumes for the iron platform the prototype is placed on *
    ************************************************************/

    // Rotation to bring the G4Tubs into the right Orientation
    G4RotationMatrix rot_vol_rotmat  = G4RotationMatrix();
    G4double rot_vol_xrot=90*deg;
    rot_vol_rotmat.rotateX(rot_vol_xrot);

    // Horizontal rotation of the platform (including prototype)
    G4RotationMatrix platform_rotmat  = G4RotationMatrix();
    double horiz_rot = fOrzrot; //set via G4UIMessenger, default=0.
    platform_rotmat.rotateY(horiz_rot);
    G4Material* platformMaterial = nistManager->FindOrBuildMaterial("G4_Fe");
    G4Tubs* iron_platform_solid = new G4Tubs("iron_platform_solid", 0, platform_radius, platform_half_height, 0., 2.*pi);

    G4LogicalVolume* iron_platform_logical = new G4LogicalVolume(iron_platform_solid,
                                                                 platformMaterial,
                                                                 "iron_platform_logical");

    // Placement of platform later utilising dimensions of calorimeter
    
    /************************************************
    * Volumes for two bars the housing is placed on *
    *************************************************/

    G4Box* outer_bar_solid = new G4Box("outer_bar_solid", bar_half_length, bar_half_height, bar_half_width);

    double subtract_bar_width = bar_half_width - bar_wall_thickness;
    double subtract_bar_height = bar_half_height - bar_wall_thickness;
    G4Box* subtract_bar = new G4Box("subtract_bar", bar_half_length, subtract_bar_height, subtract_bar_width);

    G4SubtractionSolid* bar_solid = new G4SubtractionSolid("bar_solid", outer_bar_solid, subtract_bar);

    // Placement of bars is done with the housing. Some dimensions from the housing are needed
    
    /******************************************
    * Housing and support for the calorimeter *
    *******************************************/
    // Variables needed for both housing and support
    G4RotationMatrix* unit_rotation = new G4RotationMatrix();

    // housing of calorimter
    G4Box* housing_solid = new G4Box("housing_solid", housing_half_width, housing_half_height, housing_half_length);

    double subtract_box_half_width = housing_half_width - side_wall_thickness;
    double subtract_box_half_height = housing_half_height - (top_wall_thickness+bot_wall_thickness)/2;

    // "air box" to hollow housing
    G4Box* subtract_box_solid = new G4Box("subtract_box", subtract_box_half_width, subtract_box_half_height, housing_half_length);

    G4LogicalVolume* subtract_box_logical = new G4LogicalVolume(subtract_box_solid,
                                                         defaultMaterial,
                                                         "subtract_box_logical");

    G4ThreeVector subtract_box_pos = G4ThreeVector(0, (bot_wall_thickness-top_wall_thickness)/2, 0);
    G4Transform3D subtract_box_transfm = G4Transform3D(*unit_rotation, subtract_box_pos); 
    // Placement done later with full housing logical volume

    //G4SubtractionSolid* housing_solid = new G4SubtractionSolid("housing_solid", outer_housing_solid, subtract_box_solid, unit_rotation, subtract_box_pos);


    // approximated support structure on which housing is placed

    G4Box* outer_support_solid = new G4Box("outer_support_solid", housing_half_width, support_half_height, support_half_length);

    double subtract_support_width = housing_half_width - support_wall_thickness;
    double subtract_support_height = support_half_height - support_wall_thickness;
    G4Box* subtract_support = new G4Box("subtract_support", subtract_support_width, subtract_support_height, support_half_length);

    G4SubtractionSolid* support_solid = new G4SubtractionSolid("support_solid", outer_support_solid, subtract_support);


    // Union of housing and support to only have to place one volume (easier for vertical rotation)
    double union_shift_y = -(housing_half_height+support_half_height);
    double union_shift_z = support_half_length - housing_half_length;
    G4ThreeVector union_pos = G4ThreeVector(0, union_shift_y, union_shift_z);
    G4UnionSolid* fullbox_nobars = new G4UnionSolid("fullbox_nobars", housing_solid, support_solid, unit_rotation, union_pos);


    // Union with front and back bar
    double bar_y = -(housing_half_height+2*support_half_height+bar_half_height);
    G4ThreeVector bar_front_pos = G4ThreeVector(0, bar_y, -(housing_half_length-bar_half_width-bar_pos_from_front));
    G4UnionSolid* fullbox_frontbar = new G4UnionSolid("fullbox_frontbar", fullbox_nobars, bar_solid, unit_rotation, bar_front_pos);

    G4ThreeVector bar_back_pos  = G4ThreeVector(0, bar_y, (2*support_half_length-housing_half_length-bar_half_width-53*cm));
    G4UnionSolid* fullbox_solid = new G4UnionSolid("fullbox_solid", fullbox_frontbar, bar_solid, unit_rotation, bar_back_pos);




    G4LogicalVolume* fullbox_logical = new G4LogicalVolume( fullbox_solid,
                                                            aluminiumMaterial,
                                                            "fullbox_logical");

    // Placement of iron platform because placement of calorimter should be relative to that
    double iron_plat_Y = -(caloY + bot_wall_thickness + 2*support_half_height + 2*bar_half_height + platform_half_height);
    G4ThreeVector iron_plat_pos = G4ThreeVector(0, iron_plat_Y, 0);
    G4Transform3D iron_plat_transfm = G4Transform3D(platform_rotmat*rot_vol_rotmat, iron_plat_pos);

    /*G4VPhysicalVolume* iron_platform_physical =*/ new G4PVPlacement(iron_plat_transfm,
                                                                      iron_platform_logical,
                                                                      "IronPlatform",
                                                                      stage_logical,
                                                                      false,
                                                                      0,
                                                                      fCheckOverlaps);


    // Inverse rotation of rotating_volume for calo to be pointed towards z 
    //G4RotationMatrix fullbox_rotmat = rot_vol_rotmat.inverse();
    G4RotationMatrix fullbox_rotmat = G4RotationMatrix();
    // Vertical rotation of module
    //double vert_rot = fVertRot ? -2.5*deg : 0.0*deg; //to be removed
    double vert_rot = -fVerrot; //set via G4UIMessenger, default=0.
                                //- sign needed to rotate as
                                //done at the test-beam.
    fullbox_rotmat.rotateX(vert_rot);

    
    //double fullbox_centre_Y = housing_half_height + 2*support_half_height + 2*bar_half_height - bar_pos_from_front*tan(-vert_rot); // Centre of union volume based on first solid

    //G4ThreeVector fullbox_pos = G4ThreeVector(0, 0, fullbox_shift);
    double rotation_R = housing_half_length - caloZ - plastic_cover_full_length; //distance between housing centre and calo centre
    double fullbox_X = sin(horiz_rot)*rotation_R;

    double rotation_Z = housing_half_length; 
    double fullbox_centre_Y = platform_half_height + 2*support_half_height + 2*bar_half_height + housing_half_height - bar_pos_from_front*tan(-vert_rot);
    double fullbox_shift = (-sin(vert_rot)*rotation_Z + cos(vert_rot)*fullbox_centre_Y);
    double fullbox_Y = iron_plat_Y + fullbox_shift;

    double fullbox_Z = rotation_R*(2-cos(horiz_rot));

    G4ThreeVector fullbox_pos = G4ThreeVector(fullbox_X, fullbox_Y, fullbox_Z);

    G4Transform3D fullbox_transfm = G4Transform3D(platform_rotmat*fullbox_rotmat, fullbox_pos); 
    /*G4VPhysicalVolume* fullbox_physical =*/ new G4PVPlacement(fullbox_transfm,
                                                                fullbox_logical,
                                                                "FullBox",
                                                                stage_logical,
                                                                false,
                                                                0,
                                                                fCheckOverlaps);


    // Placement of hollowing air volume for housing
    /*G4VPhysicalVolume* subtract_box_physical =*/ new G4PVPlacement(subtract_box_transfm,
                                                                subtract_box_logical,
                                                                "SubtractBox",
                                                                fullbox_logical,
                                                                false,
                                                                0,
                                                                fCheckOverlaps);


    /*****************
     * Plastic Cover *
     *****************/
    double cover_half_length = plastic_cover_full_length/2;
    double cover_half_width  = subtract_box_half_width;
    double cover_half_height = subtract_box_half_height;
    G4Box* cover_without_cutout_solid = new G4Box("cover_without_cutout_solid", cover_half_width, cover_half_height, cover_half_length);

    double cutout_half_length = cover_half_length - cutout_thickness/2;
    double cutout_half_width =  caloX;
    double cutout_half_height = caloY;
    G4Box* cutout_solid = new G4Box("cutout_solid", cutout_half_width, cutout_half_height, cutout_half_length);

    G4ThreeVector cutout_pos = G4ThreeVector(0, cutout_half_height-cover_half_height, cutout_half_length-cover_half_length);
    G4SubtractionSolid* cover_solid = new G4SubtractionSolid("cover_solid", cover_without_cutout_solid, cutout_solid, unit_rotation, cutout_pos);

    

    G4LogicalVolume* cover_logical = new G4LogicalVolume(cover_solid,
                                                         PVCMaterial,
                                                         "cover_logical");

    
    G4ThreeVector cover_pos = G4ThreeVector(0, 0, cover_half_length-housing_half_length);
    G4Transform3D cover_transfm = G4Transform3D(*unit_rotation, cover_pos); 
    /*G4VPhysicalVolume* cover_physical =*/ new G4PVPlacement(cover_transfm,
                                                                cover_logical,
                                                                "PlasticCover",
                                                                subtract_box_logical,
                                                                false,
                                                                0,
                                                                fCheckOverlaps);



    /**************
    * Calorimeter *
    ***************/

    G4RotationMatrix calo_rotmat = G4RotationMatrix();
    double fullbox_floor = housing_half_height - bot_wall_thickness;
    double calo_shift_y = -(subtract_box_half_height - caloY);
    double calo_shift_z = -(housing_half_length - caloZ - plastic_cover_full_length);
    G4ThreeVector calo_pos = G4ThreeVector(0, calo_shift_y, calo_shift_z);

    G4Transform3D calo_transfm = G4Transform3D(calo_rotmat, calo_pos);    
    /*G4VPhysicalVolume* CalorimeterPV =*/ new G4PVPlacement(calo_transfm,
                                                         CalorimeterLV,
                                                         "Calorimeter",
                                                         subtract_box_logical,
                                                         false,
                                                         0,
                                                         fCheckOverlaps);

    // Module (same shape as moduleequipped)
    //
    G4VSolid* moduleS = new G4ExtrudedSolid("module", 
		                            polygon, 
					    moduleZ/2, 
					    offA, scaleA, offB, scaleB);
                         
    G4LogicalVolume* moduleLV = new G4LogicalVolume(moduleS,
                                                    defaultMaterial,
                                                    "module");

    G4ThreeVector pos_module;
    pos_module.setX(0.);
    pos_module.setY(0.);
    pos_module.setZ(-0.18);
                              
    /*G4VPhysicalVolume* modulePV =*/ new G4PVPlacement(0,
                                                    pos_module,
                                                    moduleLV,
                                                    "module",
                                                     moduleequippedLV,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);

    // Tubes with Scintillating fibers
    //
    G4int copynumber = 0;

    for(int column=0; column<NofFiberscolumn; column++){
        
        std::stringstream S_fiber_column;
        S_fiber_column.str("");
        S_fiber_column << column;

        for(int row=0; row<NofFibersrow; row++){
            
            std::stringstream S_fiber_row;
            S_fiber_row.str("");
            S_fiber_row << row;
            std::string S_name;
            S_name = "S_column_" + S_fiber_column.str() + "_row_" + S_fiber_row.str(); 

            G4double S_x, S_y;
            G4ThreeVector vec_S_fiber;

            if(row%2==0){
                S_x = +moduleX/2 - tuberadius - (tuberadius*2+2*tolerance)*column;
                S_y = -moduleY/2 + tuberadius + (sq3*tuberadius+2*tolerance*mm)*(row)+tuberadius*(2.*sq3m1-1.);

                vec_S_fiber.setX(S_x);
                vec_S_fiber.setY(S_y);
                vec_S_fiber.setZ(0.);

                copynumber = ((NofFibersrow/2)*column+row/2);
                auto logic_S_fiber = constructscinfiber(tolerance,
                                                        tuberadius,
                                                        fiberZ,
                                                        absorberMaterial,
                                                        coreradius,
                                                        coreZ,
                                                        ScinMaterial,
                                                        claddingradiusmin,
                                                        claddingradiusmax,
                                                        claddingZ,
                                                        CherMaterial);
                // Tubes with scintillating fiber placement
                //
                /*physi_S_fiber[column][row] =*/ new G4PVPlacement(0,
                                                               vec_S_fiber,
                                                               logic_S_fiber,
                                                               S_name,
                                                               moduleLV,
                                                               false,
                                                               copynumber); 
            }
        };
    };

    // Tubes with Cherenkov fibers
    //
  
    for(int column=0; column<NofFiberscolumn; column++){
        
        std::stringstream C_fiber_column;
        C_fiber_column.str("");
        C_fiber_column << column;
        for(int row=0; row<NofFibersrow; row++){
            
            std::stringstream C_fiber_row;
            C_fiber_row.str("");
            C_fiber_row << row;
            std::string C_name;
            C_name = "C_column_" + C_fiber_column.str() + "_row_" + C_fiber_row.str(); 

            G4double C_x, C_y;
            G4ThreeVector vec_C_fiber;

            if(row%2 != 0){
                C_x = moduleX/2 - tuberadius - tuberadius - (tuberadius*2+2*tolerance)*column;
                C_y = -moduleY/2 + tuberadius + (sq3*tuberadius+2*tolerance)*row+tuberadius*(2.*sq3m1-1.);

                vec_C_fiber.setX(C_x);
                vec_C_fiber.setY(C_y);
                vec_C_fiber.setZ(0.);

                copynumber = ((NofFibersrow/2)*column+row/2);
                        
                auto logic_C_fiber = constructcherfiber(tolerance,
                                                        tuberadius,
                                                        fiberZ,
                                                        absorberMaterial,
                                                        coreradius,
                                                        coreZ,
                                                        CherMaterial,
                                                        claddingradiusmin,
                                                        claddingradiusmax,
                                                        claddingZ,
                                                        CladCherMaterial);
                /*physi_C_fiber[column][row] =*/ new G4PVPlacement(0,
                                                         vec_C_fiber,
                                                         logic_C_fiber,
                                                         C_name,
                                                         moduleLV,
                                                         false,
                                                         copynumber);
            }      
        };
    };

    // Return physical world
    //
    return fWorldPV;

}

// Define constructscinfiber method()
//
G4LogicalVolume* DRTB23SimDetectorConstruction::constructscinfiber(double tolerance,
    G4double tuberadius, G4double fiberZ, G4Material* absorberMaterial, 
    G4double coreradius, G4double coreZ, G4Material* ScinMaterial, 
    G4double claddingradiusmin, G4double claddingradiusmax, G4double claddingZ,
    G4Material* CherMaterial){
  
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()m
    std::uniform_real_distribution<> dis(0.0, tolerance);
    double outradiussmear = dis(gen);
    tuberadius = tuberadius+outradiussmear;
  
    // Tube for scintillating fibers
    //
    G4Tubs* S_fiber = new G4Tubs("S_fiber", 0., tuberadius, fiberZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_S_fiber = new G4LogicalVolume(S_fiber,
                                                         absorberMaterial,
                                                         "S_fiber");
//    logic_S_fiber->SetVisAttributes(G4VisAttributes::Invisible);
	
    G4Tubs* Abs_S_fiber = new G4Tubs("Abs_Scin_fiber", claddingradiusmax, tuberadius, fiberZ/2,0.,2.*pi);

    G4LogicalVolume* logic_Abs_S_fiber = new G4LogicalVolume(Abs_S_fiber,
                                                             absorberMaterial,
                                                             "Abs_Scin_fiber");
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                               G4ThreeVector(0.,0.,0.),
                                               logic_Abs_S_fiber,
                                               "Abs_Scin_fiber",
                                               logic_S_fiber,
                                               false,
                                               0,
                                               fCheckOverlaps);

    G4Tubs* Core_S_fiber = new G4Tubs("Core_S_fiber", 0., 
                                      coreradius, coreZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Core_S_fiber = new G4LogicalVolume(Core_S_fiber,
                                                              ScinMaterial,
                                                              "Core_S_fiber");

    G4VisAttributes* ScincoreVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.8)); //blue
    ScincoreVisAtt->SetVisibility(true);
    ScincoreVisAtt->SetForceWireframe(true);
    ScincoreVisAtt->SetForceSolid(true);
    logic_Core_S_fiber->SetVisAttributes(ScincoreVisAtt);
//    logic_Core_S_fiber->SetVisAttributes(G4VisAttributes::Invisible);
    G4ThreeVector vec_Core_S;
    vec_Core_S.setX(0.);
    vec_Core_S.setY(0.);
    vec_Core_S.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     vec_Core_S,
                                                     logic_Core_S_fiber,
                                                     "Core_S_fiber",
                                                     logic_S_fiber,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);
 

    G4Tubs* Clad_S_fiber = new G4Tubs("Clad_S_fiber", claddingradiusmin, 
        claddingradiusmax, claddingZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Clad_S_fiber = new G4LogicalVolume(Clad_S_fiber,
                                                              CherMaterial,
                                                              "Clad_S_fiber");

    G4VisAttributes* ScincladVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    //light blue
    ScincladVisAtt->SetVisibility(true);
    ScincladVisAtt->SetForceWireframe(true);
    ScincladVisAtt->SetForceSolid(true);
    logic_Clad_S_fiber->SetVisAttributes(ScincladVisAtt);
//    logic_Clad_S_fiber->SetVisAttributes(G4VisAttributes::Invisible);

    G4ThreeVector vec_Clad_S;
    vec_Clad_S.setX(0.);
    vec_Clad_S.setY(0.);
    vec_Clad_S.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     vec_Clad_S,
                                                     logic_Clad_S_fiber,
                                                     "Clad_S_fiber",
                                                      logic_S_fiber,
                                                      false,
                                                      0,
                                                      fCheckOverlaps);


    G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(0.6,0.3,0.3)); //blue
    TubeVisAtt->SetVisibility(true);
    TubeVisAtt->SetForceWireframe(true);
    TubeVisAtt->SetForceSolid(true);
    logic_Abs_S_fiber->SetVisAttributes(TubeVisAtt);
//    logic_Abs_S_fiber->SetVisAttributes(G4VisAttributes::Invisible);
    
    return logic_S_fiber;

}

// Define constructcherfiber() method
//
G4LogicalVolume* DRTB23SimDetectorConstruction::constructcherfiber(double tolerance,
    G4double tuberadius, G4double fiberZ, G4Material* absorberMaterial,
    G4double coreradius, G4double coreZ, G4Material* CherMaterial, 
    G4double claddingradiusmin, G4double claddingradiusmax, G4double claddingZ, 
    G4Material* CladCherMaterial){ 
 
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()m
    std::uniform_real_distribution<> dis(0.0, tolerance);
    double outradiussmear = dis(gen);
    tuberadius = tuberadius+outradiussmear;

    G4Tubs* C_fiber = new G4Tubs("C_fiber", 0., tuberadius, fiberZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_C_fiber = new G4LogicalVolume(C_fiber,
                                                         absorberMaterial,
                                                         "C_fiber");

//    logic_C_fiber->SetVisAttributes(G4VisAttributes::Invisible);
    G4Tubs* Abs_C_fiber = new G4Tubs("Abs_Cher_fiber", claddingradiusmax, tuberadius, fiberZ/2,0.,2.*pi);

    G4LogicalVolume* logic_Abs_C_fiber = new G4LogicalVolume(Abs_C_fiber,
                                                             absorberMaterial,
                                                             "Abs_Cher_fiber");
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     G4ThreeVector(0.,0.,0.),
                                                     logic_Abs_C_fiber,
                                                     "Abs_Cher_fiber",
                                                     logic_C_fiber,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);

    G4Tubs* Core_C_fiber = new G4Tubs("Core_C_fiber", 0., coreradius, coreZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Core_C_fiber = new G4LogicalVolume(Core_C_fiber,
                                                              CherMaterial,
                                                              "Core_C_fiber");

    G4VisAttributes* ChercoreVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
    ChercoreVisAtt->SetVisibility(true);
    ChercoreVisAtt->SetForceWireframe(true);
    ChercoreVisAtt->SetForceSolid(true);
    logic_Core_C_fiber->SetVisAttributes(ChercoreVisAtt);
//    logic_Core_C_fiber->SetVisAttributes(G4VisAttributes::Invisible);
    G4ThreeVector vec_Core_C;
    vec_Core_C.setX(0.);
    vec_Core_C.setY(0.);
    vec_Core_C.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                    vec_Core_C,
                                                    logic_Core_C_fiber,
                                                    "Core_C_fiber",
                                                    logic_C_fiber,
                                                    false,
                                                    0,
                                                    fCheckOverlaps);

    G4Tubs* Clad_C_fiber = new G4Tubs("Clad_C_fiber", claddingradiusmin,
        claddingradiusmax, claddingZ/2, 0., 2.*pi);

    G4LogicalVolume* logic_Clad_C_fiber = new G4LogicalVolume(Clad_C_fiber,
                                                              CladCherMaterial,
                                                              "Clad_C_fiber");

    G4VisAttributes* ChercladVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    //yellow 
    ChercladVisAtt->SetVisibility(true);
    ChercladVisAtt->SetForceWireframe(true);
    ChercladVisAtt->SetForceSolid(true);
    logic_Clad_C_fiber->SetVisAttributes(ChercladVisAtt);
//    logic_Clad_C_fiber->SetVisAttributes(G4VisAttributes::Invisible);

    G4ThreeVector vec_Clad_C;
    vec_Clad_C.setX(0.);
    vec_Clad_C.setY(0.);
    vec_Clad_C.setZ(0.); 
                             
    /*G4VPhysicalVolume* =*/ new G4PVPlacement(0,
                                                     vec_Clad_C,
                                                     logic_Clad_C_fiber,
                                                     "Clad_C_fiber",
                                                     logic_C_fiber,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);

    G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(0.6,0.3,0.3)); //blue
    TubeVisAtt->SetVisibility(true);
    TubeVisAtt->SetForceWireframe(true);
    TubeVisAtt->SetForceSolid(true);
    logic_Abs_C_fiber->SetVisAttributes(TubeVisAtt);
//    logic_Abs_C_fiber->SetVisAttributes(G4VisAttributes::Invisible);

    return logic_C_fiber;

}
//
//   method to define polygonal contour of module to extrude
//
std::vector<G4TwoVector> DRTB23SimDetectorConstruction::calcmod(double radius, int nrow, int ncol) {
   G4double dxlr=radius;
   G4double dtubeY=sq3*radius;
   G4double moduleX = ncol*2.*radius+radius;
   G4double side=radius*2*sq3m1;
   G4double moduleY = (nrow-1.)*dtubeY+4.*radius*sq3m1;
   int nyp=2*nrow+1;
   int nxp=2*ncol-2;
   double yp[10000];
   double xp[10000];
   yp[0]=0.;
   xp[0]=0.;
   for(int i=1;i<nyp;i++){
     yp[i]=((i-1)%2+1)*0.5*side+yp[i-1];
     xp[i]=-((i+1)/2)%2*dxlr;
   }
   for(int i=nyp;i<(nyp+nxp);i++){
     int j=i-nyp;
     yp[i]=yp[nyp-1]+0.5*(i%2)*side;
     xp[i]=xp[nyp-1]+dxlr*(j+1);
   }
   for(int i=(nyp+nxp);i<(2*nyp+nxp);i++){
     int j=i-nyp-nxp;
     yp[i]=yp[nyp-j];
     xp[i]=xp[nyp+nxp-1]+dxlr+((j+1)/2)%2*dxlr;
   }
   for(int i=(2*nyp+nxp);i<(2*nyp+2*nxp);i++){
     int j=i-2*nyp-nxp;	
     yp[i]=0.5*side*(i%2);
     xp[i]=xp[nyp+nxp-j-1];
   }
//   xp[2*nyp+2*nxp]=0.;
//   yp[2*nyp+2*nxp]=0.;
   std::vector<G4TwoVector> polygon1(2*nyp+2*nxp);
   for(int i=0;i<(2*nyp+2*nxp);i++){
     polygon1[i].set(-xp[i]+moduleX/2.-radius,yp[i]-moduleY/2.);
   }
   return polygon1;
}
