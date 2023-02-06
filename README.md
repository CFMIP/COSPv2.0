# About COSP

The CFMIP Observation Simulator Package (COSP) takes the models representation of the
atmosphere and simulates the retrievals for several passive (ISCCP, MISR and MODIS)
and active (CloudSat (radar) and CALIPSO (lidar)) sensors.

An overview of COSP is provided in the [COSP1 BAMS paper](https://doi.org/10.1175/2011BAMS2856.1).

COSP Version 2 (COSP2) is a major reorganization and modernization of the previous generation of COSP. For a detailed description, see the [COSP2 GMD paper](https://doi.org/10.5194/gmd-11-77-2018).

The simulators in COSP (ISCCP, MISR, MODIS, radar/CloudSat and lidar/CALIPSO) have been developed by many institutions and agencies:
- Met Office Hadley Centre
- LLNL (Lawrence Livermore National Laboratory)
- LMD/IPSL (Laboratoire de Meteorologie Dynamique/Institut Pierre Simon Laplace)
- CSU (Colorado State University)
- UW (University of Washington)
- CU/CIRES (University of Colorado/Cooperative Institute for Research In Environmental Sciences)

# Conditions of use

The code is distributed under [BSD License](https://opensource.org/licenses/BSD-3-Clause). Each source file includes a copy of this license with details on the Owner, Year and Organisation. The license in the file `quickbeam/README` applies to all the files in the directory quickbeam.

# What's in the distribution

The repository include directories
- `src/` contains the COSP layer and the underlying satellite simulators
- `model-interface/` contains routines used by COSP to couple to the host model. Edit these before building.
- `subsample_and_optics_example/` contains an example implementation, following COSP 1, of the routines to map model-derived grid-scale physical properties to the subgrid-scale optical properties needed by COSP.
- `driver/` contains codes that run COSP on example inputs and scripts that compare the current implementation to a reference.
- `build/` contains a Makefile describing build dependencies. Users may build a COSP library and other targets from this Makefile.
- `unit_testing/` contains small programs meant to exercise some of the simulators.

Users incorporating COSP into a model will need all routines found within `src/`, appropriately-edited versions of the routines in `model-interface/`, and routines that provide the functionality of those in `subsample_and_optics_example/`.

As described in [Swales et al. (2018)](https://doi.org/10.5194/gmd-11-77-2018), COSP2 requires inputs in the forms of subcolumn-sampled optical properties. The drivers map the model physical state to these inputs using the routines in `driver/subsample_and_optics/`, which are consistent with the fixed choices made in COSP1. We anticipate that users incorporating COSP into models will develop a model-specific mapping between the model's physical state and the inputs required for COSP that is consistent with the host model.

The offline drivers read sample snapshots from the Met Office Unified Model, use the routines in `subsample_and_optics_example/` to compute COSP inputs, and record the results from COSP in netCDF files. The default driver calls COSP 2 directly and produces netCDF output. The layer mimicking the COSP 1.4 interface is tested with a separate driver. A third driver uses the CMOR1 infrastructure to write a subset of fields to individual netCDF files following conventions for publication on the Earth System Grid Federation.

# Running the offline tests

 1. Build the drivers.
    1. Edit the files in `model-interface/` if necessary. By default COSP is built using double-precision real variables and printing any error messages to the standard output.
    1. In `build/` edit `Makefile.conf` to reflect the choice of compiler, compiler flags, and library names and locations.
    1. In `build/`, `make driver` will build a COSP library, a separate library with the example mapping from model state to COSP inputs, and the `cosp2_test` executable, which is then copied to `driver/run`.

1. Running the test program
    1. Directory `driver/run` contains namelists and other files needed by the test programs. If the executables have been built they should run in this directory using these files as supplied.
    1. The behavior of COSP can be changed via the input namelists (e.g. `driver/src/cosp2_input_nl.txt`) and output (`driver/src/cosp2_output_nl.txt`) namelists. The input namelist controls the COSP setup (i.e. Number of subcolumns to be used, etc...) and
    simulator specific information (i.e. Radar simulator frequency). The output namelist controls the output diagnostics. The test program receives the input namelist as an argument from the command line, i.e.:

      `./cosp2_test cosp2_input_nl.um_global.txt`

Currently, there are 2 input namelists: `cosp2_input_nl.txt` uses a small test input file (`cosp_input.nc`), stored in the github repository; `cosp2_input_nl.um_global.txt` uses a global field (`cosp_input.um_global.nc`). The global input file and its associated known good output (`cosp2_output.um_global.gfortran.kgo.nc`) are stored externally in Google Drive, and they can be retrieved by running `download_test_data.sh` from within the `driver/` directory.

1. Regression testing (comparing to reference data)
   1. Reference data or known good output (KGO) for a small test case is provided with COSP2. The data can be found at `driver/data/outputs/UKMO/`.
   1. To compare driver output to the KGO, in `driver/`, invoke Python 3 script `compare_to_kgo.py`. This script requires the following Python modules: numpy, netCDF4, argparse, and sys. The following example shows how to call this script from the command line:

      `python compare_to_kgo.py data/outputs/UKMO/cosp2_output_um.gfortran.kgo.nc data/outputs/UKMO/cosp2_output_um.nc --atol=1.0e-20 --rtol=0.0005`

   The script accepts thresholds for absolute and relative tolerances, named `atol` and `rtol` respectively. By default the script will report all differences, i.e. `--atol=0.0 --rtol=0.0`. The example above allows for a relative tolerance of 0.05\% in the subset of absolute differences larger or equal than 1.0e-20. Previous tests indicate that these thresholds are appropriate when gfortran is used to build the driver program. If a different compiler is available, it is encouraged to build the driver using that compiler to increase the robustness of the tests. For non-gfortran compilers, the differences may be larger. This is not necessarily an issue, but further investigations may be required.
