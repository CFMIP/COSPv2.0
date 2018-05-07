%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1) About COSP

The CFMIP Observation Simulator Package (COSP) takes the models representation of the
atmosphere and simulates the retrievals for several passive (ISCCP, MISR and MODIS)
and active (CloudSat (radar) and CALIPSO (lidar)) sensors.

An overview of COSP is provided in https://doi.org/10.1175/2011BAMS2856.1

COSP Version 2 (COSP2) is a major reorganization and modernization of the previous
generation of COSP. For a detailed description, see https://doi.org/10.5194/gmd-2017-148

The simulators in COSP (ISCCP, MISR, MODIS, radar/CloudSat and lidar/CALIPSO) have been
developed by many institution and agencies:
*) Met Office Hadley Centre
*) LLNL (Lawrence Livermore National Laboratory)
*) LMD/IPSL (Laboratoire de Meteorologie Dynamique/Institut Pierre Simon Laplace)
*) CSU (Colorado State University)
*) UW (University of Washington)
*) CU/CIRES (University of Colorado/Cooperative Institute for Research In Environmental Sciences)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2) CONDITIONS OF USE

The code is distributed under BSD License (http://www.opensource.org/licenses/bsd-license.php).
Each source file includes a copy of this license with details on the Owner, Year
and Organisation. The license in the file quickbeam/README applies to all the files in
the directory quickbeam.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3) What's in the distribution

The repository include directories
  src/ contains the COSP layer and the underlying satellite simulators
  model-interface/ contains routines used by COSP to couple to the host model.
    Edit these before building.
  subsample_and_optics_example/ contains an example implementation, following
    COSP 1, of the routines to map model-derived grid-scale physical properties
    to the subgrid-scale optical properties needed by COSP
  cosp-1.4-interface/ contains a wrapper mapping the calling structure and arguments
    from COSP 1 to COSP 2, making it possible to call COSP 2 using COSP 1 arguments.
  driver/ contains codes that run COSP on example inputs and scripts that compare
    the current implementation to a reference.
  build/ contains a Makefile describing build dependencies. Users may build a COSP
    library and other targets from this Makefile.
  unit_testing/ contains small programs meant to exercise some of the simulators

Users incorporating COSP into a model will need all routines found within src/,
appropriately-edited versions of the routines in model-interface/, and routines
that provide the functionality of those in subsample_and_optics_example/.

As described in https://doi.org/10.5194/gmd-2017-148 COSP2 requires inputs in the
forms of subcolumn-sampled optical properties. The drivers map the model physical
state to these inputs using the routines in driver/subsample_and_optics/, which
are consistent with the fixed choices made in COSP1. We anticipate that users
incorporating COSP into models will develop a model-specific mapping between the
model's physical state and the inputs required for COSP that is consistent with
the host model.

The offline drivers read sample snapshots from the Met Office Unified Model, use
the routines in subsample_and_optics_example/ to compute COSP inputs, and record
the results from COSP in netCDF files. The default driver calls COSP 2 directly
and produces netCDF output. The layer mimicking the COSP 1.4 interface is tested
with a separate driver. A third driver uses the CMOR1 infrastructure to write a
subset of fields to individual netCDF files following conventions for publication
on the Earth System Grid Federation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4) Running the offline tests

  a) Build the drivers.
    *) Edit the files in model-interface/ if necessary. By default COSP is built
    using double-precision real variables and printing any error messages to the
    standard output.
    *) In build/ edit Makefile.conf to reflect the choice of compiler, compiler flags,
    and library names and locations. If you intend to build the CMOR driver this
    includes all libraries necessary for CMOR. Building the CMOR driver also
    requires a change to the definition of the DRIVER variable in Makefile.
    *) In build 'make driver' will build a COSP library, a separate library with
    the example mapping from model state to COSP inputs, and the cosp2_test
    executable, which is then copied to driver/run. 'make driver_COSP1.4' is
    analogous but builds a cosp1_test executable that uses the COSP 1.4 calling
    conventions.

  b) Running the test program
    *) Directory test/run contains namelists and other files needed by the test
    programs. If the executables have been built they should run in this
    directory using these files as supplied.
    *) The behavior of COSP can be changed via the input (driver/src/cosp2_input_nl.txt)
    and output (driver/src/cosp2_output_nl.txt) namelists. The input namelist
    controls the COSP setup (i.e. Number of subcolumns to be used, etc...) and
    simulator specific information (i.e. Radar simulator frequency). The output
    namelist controls the output diagnostics.

  c) Regression testing (comparing to reference data)
    *) Reference data for a small test case is provided with COSP2. The data can be
    found at driver/data/outputs/UKMO/. CMOR compliant reference data is also provided,
    driver/data/outputs/UKMO/cmor/ref1D/.
    *) To compare driver output to reference data. In driver/, invoke Python script
    test_cosp2imp.py. This script requires the following Python modules: os, numpy, netCDF4,
    argparse, warnings, fnmatch, sys. Examples are below.
       -) For standard netCDF output:
       python test_cosp2Imp.py data/outputs/UKMO/cosp2_output_um.ref.nc data/outputs/UKMO/cosp2_output_um.nc
       -) For CMOR compliant output:
       python test_cosp2Imp.py data/outputs/UKMO/cmor/ref1D/ data/outputs/UKMO/cmor/1D --cmor 1D
    By default the script will only report relative differences which are greater than 1e-5. This
    can be changed by the user through the optional argument "--zeroThresh".
