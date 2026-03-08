# hibeam_g4

Geant4 program to simulate the HIBEAM detector and WASA.

## Prerequisites

* cmake
* Geant4 (v11.0 or newer) https://geant4.web.cern.ch/
* ROOT (6.14 or newer) https://root.cern/
* VMC (tested with v5.2 and newer) https://vmc-project.github.io/vgm-documentation/

## Installation

Go to the path you want to install hibeam_g4 and clone the repository:

```
cd /install/path
git clone https://github.com/HIBEAM-NNBAR/hibeam_g4.git 
```

Create a build directory. In the build directory, first run cmake then compile:

```
mkdir hibeam_g4_build
cd hibeam_g4_build
cmake -DVGM_DIR=/path/to/vgm_install/lib64/VGM-X.Y.Z ../hibeam_g4
make
```

## Running

Run the program with
```
./hibeam_g4 [OPTIONS]
```

Supported options are:

* [-h] [--help] : provide command line help.
* [-g] [--gui] : set flags for graphical interface (if Geant4 is built with QT).
* [-t] [--threads] : set number of threads to be used (default=1).
* [-i inputfile] [--input inputfile] : provide the inputfile i.e. the input file from external event generators. Be aware that this need the PrimaryGeneratorAction for your input file to be provided, included and compiled into the Geant4 code.
* [-m run.mac] [--mac run.mac] : provide Geant4 macro file for Geant4 settings (for example number of events generated).
* [-c config.par] [--config config.par] : configuration file for the code, can set different option that are implemented within the simulation code.
* Outputfile.root : mandatory ! provide name of ROOT output file.

## Config parameters

A number of parameters of the simulation can be set using the configuration file. Supported parameters are:

* Gui : set flags for graphical interface (if Geant4 is built with QT). Same as command line option "-g"
* Threads = set number of threads to be used (default=1). Same as command line option "-t"
* Gui : Legacy parameter. Has to be set to "Wasa".
* Geometry_Namefile : Path to root or gdml file containing simulation geometry. Default="geometry.root".
* CheckOverlaps : Flag to check loaded geometry for overlapping volumes. Default=0 
* Source : Define type of particle source to be used. Options are  "gps"/"GPS" - Geant4 general particle source, "mcpl"/"MCPL" - Read particles from input MCPL file, and "scattering"/"Scattering" - Generate proton-deuteron pair according to elastic scattering kinematics. Further parameters for each source have to be set using a macro file. Default="gps".
* MCPL_Inputfile : Path to mcpl file containing particles. Necessary for MCPL source. Default="none".
* Detectors : Comma-separated list of volumes which are supposed to be treated as sensitive detectors in the simulation. Default=empty.
* Sampling_Detectors : Comma-separated list of volumes which are supposed to be treated as sampling detectors in the simulation. Sampling detectors can be used to output particle position and direction at a certain point in the setup. Default=none.
* WriteHistograms : Flag to write flux histograms to the output file. Default=1.
* WriteTree : Flag to write root tree containing particle trajectories in sensitive detecors and sampling detectors to the output file. Can lead to large output files! Default=0.

An example of a configuration file can be found in the macros folder. The folder also contains example macros for setting up the different particle sources.


