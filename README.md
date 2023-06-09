# DRTB23Sim
**A Geant4 simulation of the 2023 Dual-Readout em-sized tubes prototype beam test.**

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-description">Project description</a></li>
     <li>
      <a href="#documentation-and-results">Documentation and results</a>
      <ul>
        <li><a href="#selected-presentations">Selected presentations</a></li>
      </ul>
    </li>
    <li><a href="#available-datasets-and-analyses">Available datasets and analyses</a></li>
    <li>
      <a href="#how-to">How to</a>
      <ul>
        <li><a href="#build-compile-and-execute-on-maclinux">Build, compile and execute on Mac/Linux</a></li>
        <li><a href="#build-compile-and-execute-on-lxplus">Build, compile and execute on lxplus</a></li>
      </ul>
    </li>
     </li><li><a href="#my-quick-geant4-installation">My quick Geant4 installation</a></li>
  </ol>                                           
</details>

<!--Project desription-->
## Project description
The project targets a standalone Geant4 simulation of the dual-readout electromagnetic-sized calorimeter 2023 test beam at the CERN-SPS H8 beam line.

<!--Documentation and results-->
## Documentation and results

### Selected presentations
- To be added.

<!--Available datasets and analyses-->
### Available datasets and analyses
- To be added.

<!--How to:-->
## How to

### Build, compile and execute on Mac/Linux
1. git clone the repo
   ```sh
   git clone https://github.com/DRCalo/DRTB23Sim.git
   ```
2. source Geant4 env
   ```sh
   source /relative_path_to/geant4-install/bin/geant4.sh
   ```
3. cmake build directory and make (example using geant4.10.07_p01)
   ```sh
   mkdir build && cd build
   cmake -DGeant4_DIR=/absolute_path_to/geant4.10.07_p01-install/lib/Geant4-10.7.1/ relative_path_to/DRTB23Sim/
   make
   ```
4. execute (example with DRTB23Sim_run.mac macro card, 2 thread, FTFP_BERT physics list and no optical propagation)
   ```sh
   ./DRTB23Sim -m DRTB23Sim_run.mac -t 2 -pl FTFP_BERT -opt false
   ```
Parser options
   * -m macro.mac: pass a Geant4 macro card 
   * -t integer: pass number of threads for multi-thread execution (example -t 3, default t 2)
   * -pl Physics_List: select Geant4 physics list (example -pl FTFP_BERT)
   * -opt FullOptic: boolean variable to switch on (true) the optical photon propagation in fibers (example -opt true, default false)

### Build, compile and execute on lxplus
1. git clone the repo
   ```sh
   git clone https://github.com/DRCalo/DRTB23Sim.git
   ```
2. cmake, build directory and make (example using geant4.10.07_p01, check for gcc and cmake dependencies for other versions)
   ```sh
   mkdir build && cd DREMTubes-build
   source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8.3.0/x86_64-centos7/setup.sh 
   source /cvmfs/geant4.cern.ch/geant4/10.7.p01/x86_64-centos7-gcc8-optdeb-MT/CMake-setup.sh 
   export CXX=`which g++`
   export CC=`which gcc`
   cmake3 -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.7.p01/x86_64-centos7-gcc8-optdeb-MT/lib64/Geant4-10.7.1 ../DRTB23Sim/
   make
   ```
3. execute (example with DRTB23Sim_run.mac macro card, 2 threads and FTFP_BERT physics list)
   ```sh
   ./DRTB23Sim -m DRTB23Sim_run.mac -t 2 -pl FTFP_BERT
   ```

<!--My quick Geant4 installation-->
## My quick Geant4 installation
Here is my standard Geant4 installation (example with Geant4.10.7.p01) starting from the unpacked geant4.10.07.tar.gz file under the example path "path/to".

1. create build directory alongside source files
      ```sh
   cd /path/to
   mkdir geant4.10.07-build
   cd geant4.10.07-build
   ```
2. link libraries with CMAKE (example with my favourite libraries)
   ```sh
   cmake -DCMAKE_INSTALL_PREFIX=/Users/lorenzo/myG4/geant4.10.07_p01-install \
   -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_BUILD_MULTITHREADED=ON \
   -DGEANT4_USE_GDML=ON ../geant4.10.07.p01
   ```
3. make it (using N threads, e.g. -j 4)
   ```sh
   make -jN
   make install
   ```
