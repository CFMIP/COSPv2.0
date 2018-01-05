%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1) About COSP

The CFMIP Observation Simulator Package (COSP) takes the models representation of the
atmosphere and simulates the retrievals for several passive (ISCCP, MISR and MODIS) and active
(CLUDSAT(radar) and CALIPSO(lidar)) sensors.

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
  src/ containing the COSP layer and the underlying satellite simulators
  driver/ containing codes that run COSP on sample inputs and write results
  testing/ containing scripts that compare results from the local installation with reference answers

Users incorporating COSP into a model should include all routines found within src/. The Makefile within this directory
describes the dependencies.

Two offline drivers are provided. Both read the same input files; these come from the Met Office
Unified Model and contain snapshots of the model's state. One driver records every field output from
COSP in a single netCDF file; the other uses the CMOR infrastructure to write a subset of fields to individual
netCDF files following conventions for publication on the Earth System Grid Federation.

As described in https://doi.org/10.5194/gmd-2017-148 COSP2 requires inputs in the forms of
subcolumn-sampled optical properties. The drivers map the model physical state to these inputs using
the routines in driver/subsample_and_optics/, which are consistent with the fixed choices made in COSP1.

We anticipate that users incorporating COSP into models will develop a model-specific mapping between the
model's physical state and the inputs required for COSP that is consistent with the host model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4) Running the offline tests

a) Compilation
      For offline testing, a driver is provided with COSP2. This driver can be built with the
      capability of creating CMOR compliant netcdf output, or output to a single netCDF file.

      *) single netCDF output
        cd driver/src
        make cosp     ! Build COSP source code
        make          ! Build driver for COSP

      *) output using CMOR
        This driver requires the CMOR2 libraries. (www2-pcmdi.llnl.gov/cmor/download/)
        cd driver/src/cmor
        make cosp     ! Build COSP source code
        make          ! Build driver for COSP

b) Running
   i)  Set up COSP input (driver/src/cosp2_input_nl.txt) and output (driver/src/cosp2_output_nl.txt)
       namelists. The input namelist controls the COSP setup (i.e. Numebr of subcolumns to be
       used, etc...) and simulator specific information (i.e. Radar simulator frequency). The
       output namelist controls the output diagnostics.

   ii) Run test code
       ./cosp2_test
       This will run COSP and create outputs for the variables selected in the output namelist.

c) Compare to reference data.
   Provided with COSP2 is a small reference dataset, for both CMOR and non-CMOR compliant output,
   along with regression tools for comparison with your output.

   i) testing/test_cosp2Imp.py:
       python test_cosp2Imp.py ../driver/data/outputs/UKMO/cosp2_output_um.ref.nc ../driver/data/outputs/UKMO/cosp2_output_um.nc

   ii) testing/test_cosp2Imp.cmor.py: (SAME as above, but for CMOR compliant output)
       python test_cospImp.py ../driver/data/outputs/UKMO/cmor/ref1D/ ../driver/data/outputs/UKMO/cmor/1D/ 1D
