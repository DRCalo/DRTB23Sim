##**************************************************
## \file CMakeLists.txt
## \brief: CMakeLists of DRTB23Sim project
## \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
## @lopezzot
## \start date: 9 June 2023
##**************************************************

# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(G4DRTB23Sim)

# Find Geant4 package, activating all available UI and Vis drivers by default
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and prohect include directories
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Locate project source and headers
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable and link it to the Geant4 libraries
#
add_executable(DRTB23Sim DRTB23Sim.cc ${sources} ${headers})
target_link_libraries(DRTB23Sim ${Geant4_LIBRARIES})

# Copy all scripts to the build directory 
#
set(DRTB23Sim_SCRIPTS
  DRTB23Sim_gui.mac
  DRTB23Sim_init_vis.mac
  DRTB23Sim_run.mac
  )

foreach(_script ${DRTB23Sim_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#For internal Geant4 use 
#
add_custom_target(G4DRTB23Sim DEPENDS DRTB23Sim)

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS DRTB23Sim DESTINATION bin)

##**************************************************
