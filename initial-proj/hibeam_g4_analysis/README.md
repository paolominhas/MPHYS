# hibeam_g4_analysis

Program to convert data of simulations using hibeam_g4 (https://github.com/HIBEAM-NNBAR/hibeam_g4).

## Prerequisites

* cmake
* ROOT (6.14 or newer)

## Installation

Go to the path you want to install hibeam_g4_analysis and clone the repository:

```
cd /install/path
git clone git@github.com:HIBEAM-NNBAR/hibeam_g4_analysis.git
```

Create a build directory. In the build directory, first run cmake then compile:

```
mkdir hibeam_g4_analysis_build
cd hibeam_g4_analysis_build
cmake ../hibeam_g4_analysis
make
```

## Running

Run the program with
```
./hibeam_ana --in=/path/to/input/file/from/hibeam_g4.root --out=/path/to/output/file.root
```

The program detects automatically data from which detectors is in the input file based on the root branch names. Supported detectors/branches are:
* ``Primary*``: Source/particle generator informatio.
* ``TARGET*``: Particle hits in HIBEAM target.
* ``TPC*``: Data from HIBEAM/NNBAR TPC.
* ``ProtoTPC*``: Data from TPC prototype at LU.
* ``SECE*``: Data from WASA calorimeter.
* ``CV_bar*``: Data from HIBEAM/NNBAR cosmic veto detector.
* ``Sci_bar*``: Data from HIBEAM/NNBAR HRD scintillators with WLS readout.
* ``Scintillator*``: Data from generic scintillator detector.
