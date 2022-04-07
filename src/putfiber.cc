
    // Optical Surface properties between the glass and the Si of the SiPM
    G4OpticalSurface* OpSurfaceGlassSi = new G4OpticalSurface("OpSurfaceGlassSi");
    OpSurfaceGlassSi -> SetType(dielectric_metal);
    OpSurfaceGlassSi -> SetModel(glisur);
    OpSurfaceGlassSi -> SetFinish(polished);
    G4double efficiencyOpSurfaceGlassSi[ENTRIES] =     //100% detection efficiency 
                                    { 1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1,
                                      1, 1, 1, 1 };

    /*G4double efficiencyOpSurfaceGlassSi[ENTRIES] =     //0% detection efficiency 
                                    { 0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0 };*/

    G4double reflectivityOpSurfaceGlassSi[ENTRIES] =  // 0% reflection
                                    { 0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0. };

    G4MaterialPropertiesTable* MPTOpSurfaceGlassSi = new G4MaterialPropertiesTable();
    MPTOpSurfaceGlassSi -> AddProperty("EFFICIENCY", 
        photonEnergy, efficiencyOpSurfaceGlassSi, ENTRIES)->SetSpline(true);
    MPTOpSurfaceGlassSi -> AddProperty("REFLECTIVITY", 
            photonEnergy, reflectivityOpSurfaceGlassSi, ENTRIES)->SetSpline(true);
    OpSurfaceGlassSi -> SetMaterialPropertiesTable(MPTOpSurfaceGlassSi);

    // SiPM
    //
    G4VSolid* SiPMS = new G4Box("SiPM", SiPMX/2, SiPMY/2, SiPMZ/2);
                         
    G4LogicalVolume* SiPMLV = new G4LogicalVolume(SiPMS, GlassMaterial,"SiPM");

    // Here I build the Si of the SiPM
    // 
    G4VSolid* SiS = new G4Box("Si", SiX/2, SiY/2, SiZ/2);
                         
    G4LogicalVolume* SiLV = new G4LogicalVolume( SiS, SiMaterial, "Si");

    // Si placement inside SiPM
    //
    G4ThreeVector vec_Si;
    vec_Si.setX(0.);
    vec_Si.setY(0.);
    vec_Si.setZ(SiPMZ/2-SiZ/2); // Si at the end of SiPM
                             
    /*G4VPhysicalVolume* SiPV =*/ new G4PVPlacement(0,
                                                vec_Si,  
                                                SiLV,
                                                "Si",
                                                SiPMLV,
                                                false,
                                                0,
                                                fCheckOverlaps);
 
    G4VisAttributes* SiVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0)); //green
    SiVisAtt->SetVisibility(true);
    SiVisAtt->SetForceWireframe(true);
    SiVisAtt->SetForceSolid(true);
    SiLV->SetVisAttributes(SiVisAtt);

    // Logical Skin Surface placement around the silicon of the SiPM
    //
    /*G4LogicalSkinSurface* OpsurfaceSi =*/ new G4LogicalSkinSurface("OpsurfaceSi", 
        SiLV, OpSurfaceGlassSi);

    // Optical Surface properties between the scintillating fibers
    // and the default material
    // I'm trying to define an optical surface completly blacked 
    // as if we absorb the light at one end of fibers
    //
    G4OpticalSurface* OpSurfacedefault = new G4OpticalSurface("OpSurfacedefault");
    OpSurfacedefault -> SetType(dielectric_dielectric);
    OpSurfacedefault -> SetModel(unified);
    OpSurfacedefault -> SetFinish(polishedbackpainted); 
    // Painted from inside the fibers, light is absorbed

    // Tubes with scintillating fibers and SiPM next to them
    //
    // Attention: I place an optical surface painted (blacked) from the moduleequippedPV 
    // to the SiPMPV, in so doing I completly avoid any cross talk between SiPMs
    //
    //G4VPhysicalVolume* physi_S_fiber[NofFibersrow][NofFiberscolumn];
    //G4VPhysicalVolume* physi_SiPM[NofFibersrow][NofFiberscolumn];  
    //G4LogicalBorderSurface* logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn];
		
    G4int copynumber = 0;

    for(int row=0; row<NofFibersrow; row++){
        
        std::stringstream S_fiber_row;
        S_fiber_row.str("");
        S_fiber_row << row;

        for(int column=0; column<NofFiberscolumn; column++){
            
            std::stringstream S_fiber_column;
            S_fiber_column.str("");
            S_fiber_column << column;
            std::string S_name;
            std::string SiPM_name;
            S_name = "S_row_" + S_fiber_row.str() + "_column_" + S_fiber_column.str(); 
            SiPM_name = "S_SiPM"; 
            //SiPM_name = "SiPMS_row" + S_fiber_row.str() + "_column_" + S_fiber_column.str();

            G4double S_x, S_y;
            G4ThreeVector vec_S_fiber;
            G4ThreeVector vec_SiPM;

            if(column%2==0){
                S_x = +moduleX/2 - tuberadius - (tuberadius*2+2*tolerance)*row;
                S_y = -moduleY/2 + tuberadius + (1.733+2*tolerance*mm)*(column);
            
                vec_S_fiber.setX(S_x);
                vec_S_fiber.setY(S_y);
                vec_S_fiber.setZ(0.);

                vec_SiPM.setX(S_x);
                vec_SiPM.setY(S_y);
                vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);
            
                copynumber = ((NofFiberscolumn/2)*row+column/2);

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
                /*physi_S_fiber[row][column] =*/ new G4PVPlacement(0,
                                                               vec_S_fiber,
                                                               logic_S_fiber,
                                                               S_name,
                                                               moduleLV,
                                                               false,
                                                               copynumber); 

                // SiPM placement
                //
                /*physi_SiPM[row][column] =*/ new G4PVPlacement(0,
                                                            vec_SiPM,
                                                            SiPMLV,
                                                            SiPM_name,
                                                            moduleequippedLV,
                                                            false,
                                                            copynumber); //same copynumber of fibers 
          
                /*logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn] =
                    new G4LogicalBorderSurface("logic_OpSurface_defaultAir",
                                               CalorimeterPV, 
                                               physi_SiPM[row][column],
                                               OpSurfacedefault);*/
            }
        };
    };

    // Tubes with Cherenkov fibers and SiPM next to them
    //
    //G4VPhysicalVolume* physi_C_fiber[NofFibersrow][NofFiberscolumn];
  
    for(int row=0; row<NofFibersrow; row++){
        
        std::stringstream C_fiber_row;
        C_fiber_row.str("");
        C_fiber_row << row;
        for(int column=0; column<NofFiberscolumn; column++){
            
            std::stringstream C_fiber_column;
            C_fiber_column.str("");
            C_fiber_column << column;
            std::string C_name;
            std::string SiPM_name;
            C_name = "C_row_" + C_fiber_row.str() + "_column_" + C_fiber_column.str(); 
            SiPM_name = "C_SiPM"; 
            //SiPM_name = "SiPMC_row" + C_fiber_row.str() + "_column_" + C_fiber_column.str();

            G4double C_x, C_y;
            G4ThreeVector vec_C_fiber;
            G4ThreeVector vec_SiPM;

            if(column%2 != 0){
                C_x = moduleX/2 - tuberadius - tuberadius - (tuberadius*2+2*tolerance)*row;
                C_y = -moduleY/2 + tuberadius + (1.733+2*tolerance*mm)*column;
         
                vec_C_fiber.setX(C_x);
                vec_C_fiber.setY(C_y);
                vec_C_fiber.setZ(0.);

                vec_SiPM.setX(C_x);
                vec_SiPM.setY(C_y);
                vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);

                copynumber = ((NofFiberscolumn/2)*row+column/2);
                        
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
                /*physi_C_fiber[row][column] =*/ new G4PVPlacement(0,
                                                         vec_C_fiber,
                                                         logic_C_fiber,
                                                         C_name,
                                                         moduleLV,
                                                         false,
                                                         copynumber);

                /*physi_SiPM[row][column] =*/ new G4PVPlacement(0,
                                                        vec_SiPM,
                                                        SiPMLV,
                                                        SiPM_name,
                                                        moduleequippedLV,
                                                        false,
                                                        copynumber); //same copynumber of fiber 

                /*logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn] =
                    new G4LogicalBorderSurface("logic_OpSurface_defaultAir",
                                               CalorimeterPV, 
                                               physi_SiPM[row][column],
                                               OpSurfacedefault);*/
            }      
        };
    };

    // Return physical world
    //
    return fWorldPV;

}
