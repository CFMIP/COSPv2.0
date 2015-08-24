
General notes:
vN.N is version number (e.g. v1.0)
We will assume that the software will be installed in ~/cosp.

This README has been written by Alejandro Bodas (alejandro.bodas at metoffice.gov.uk).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1.- ABOUT THE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COSP stands for CFMIP Observation Simulator Package.
This is a code that takes information from atmospheric models as input
and simulates the signal of CloudSat and CALIPSO. It computes summary
statistics from these outputs that can be compared with similar statistics
computed from the observations. It also includes the ISCCP simulator.
There are several groups that have been actively involved in the code development:

  * Met Office Hadley Centre

  * LLNL (Lawrence Livermore National Laboratory)

  * LMD/IPSL (Laboratoire de Meteorologie Dynamique/Institut Pierre Simon Laplace)

  * CSU (Colorado State University)

  * UW (University of Washington)

The approach is to create a modular code in FORTRAN90 (some lower level 
routines are in F77) that can be plugged in different types of models, 
from general circulation models (GCMs) to cloud cloud-resolving models (CRMs).
The logical structure of the code is split into 4 sections:

  * Sub-grid cloud distribution (SCOPS:Subgrid Cloud Overlap Profile Sampler,
    which is already used by the ISCCP simulator)

  * Sub-grid precipitation.

  * Instrument simulators:
      ISCCP simulator
      Radar simulator (QuickBeam)
      Lidar and PARASOL simulators (ACTSIM)
      MISR simulator
      MODIS simulator
      RTTOV simulator

  * Summary statistics diagnostics 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2.- CONDITIONS OF USE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The code is distributed under BSD License (http://www.opensource.org/licenses/bsd-license.php).
Each source file includes a copy of this license with details on the Owner, Year 
and Organisation. The license in the file quickbeam/README applies to all the files in 
the directory quickbeam.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3.- DOWNLOADING AND UNPACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Since v1.3.1, the code lives in Google code:
http://code.google.com/p/cfmip-obs-sim/

We encourage you to download the latest stable release using a subversion client on your local computer.
You can accomplish this via the command line at a terminal prompt:

svn co http://cfmip-obs-sim.googlecode.com/svn/stable/current COSP

When COSP is stable we will update this directory, so you can update your own code using 'svn update'. 

The input and test NetCDF files have to be downloaded separately, from:
http://code.google.com/p/cfmip-obs-sim/downloads/list

Copy cosp-test-files-vN.N.N.tgz to your COSP directory and run:
tar -txzvf cosp-test-files-vN.N.N.tgz

Now you can follow the steps in Section 4.

Older versions: 1.3.0 and before
------------------------------------
These versions are not available in the Google code svn repository, only in the CFMIP website.
a) Download cosp.vN.N.tar.gz in your system. I'll assume that you have chosen ~/cosp as your
   target directory.
b) cd ~/cosp
c) Extract the source code
   tar -xvzf cosp.vN.N.tar.gz

You may want to delete the *.tar file as it is no longer needed.

This creates the directory COSP, where the source code for COSP is found.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4.- COMPILATION AND TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

4.1 - Compilation
--------------------------------------------------------------------------------------------

In order to build COSP you'll need the following libraries:
  a) HDF5 libraries (http://www.hdfgroup.org/HDF5/)
  b) NetCDF4 libraries (http://www.unidata.ucar.edu/software/netcdf/)
  c) CMOR2 libraries (http://www2-pcmdi.llnl.gov/cmor/)

Before building COSP, you will need to edit the Makefile above the following lines
 
########################################################################
#              End of modifications
########################################################################

You will need to change the environment variables so that they point to where the required 
libraries are located in your system.

You can still use CMORv1.3. Please, see section 6.7 for details.

I reccommend to compile the code first with no optimisatons and run the tests. Then, once 
you are confident that is working properly in your system you can increase the level of
optimisations.

Compilation:
make clean
make

This will create the executable cosp_test.

4.2 - Testing
--------------------------------------------------------------------------------------------
a) Setting up the namelists
The input parameters for cosp_test are provided through the  namelist COSP_INPUT,
which can be found in the file cosp_input_nl.txt. You can run two tests by changing
the input file:

 FINPUT='cosp_input_um.nc', ! NetCDF file with 1D inputs
 or
 FINPUT='cosp_input_um_2d.nc', ! NetCDF file with 2D inputs

You will also need to modify the CMOR namelist
in ~/cosp/cosp.vN.N/cmor/cosp_cmor_nl.txt. Specifically, you will need to use
   
   TABLE = 'COSP_table_1D' for 'cosp_input_um.nc'
   TABLE = 'COSP_table_2D' for 'cosp_input_um_2d.nc'
   
COSP_table_2D contains the definitions of the variables in a lon/lat grid.

b) Run cosp_test. This will create the NetCDF outputs in ~/cosp/cosp.vN.N/outputs.
c) Compare these files with the expected outputs, located in:
  ~/cosp/cosp.vN.N/outputs_test
The script cosp_test.sh can be used to check that your outputs are equal to 
the ones in outputs_test. The usage is as follows:
   ./cosp_test.sh arg1 arg2
   arg1: 1D or 2D
   arg2 must be one of these options
     0: use cmp
     1: use nccmp
     2: use nco + diff
     3: use nco to create differences

NOTE: If you get the message "command not found" when trying to run cosp_test.sh, 
change the first line (#!/usr/bin/sh) to point where sh is located in your system.

Option 1 uses nccmp, which is available from http://nccmp.sourceforge.net/. Options 2 and 3 use
the NetCDF Operators (NCO), available from http://nco.sourceforge.net/.
     
I recommend using option 1 first (nccmp). The following lines show an example of output using nccmp,
where the test found differences in two files, meantbisccp_2D.nc and tclisccp_2D.nc.

    Comparing  albisccp_2D.nc ...
    Comparing  atb532_2D.nc ...
    Comparing  boxptopisccp_2D.nc ...
    Comparing  boxtauisccp_2D.nc ...
    Comparing  cfad_dbze94_2D.nc ...
    Comparing  cfad_lidarsr532_2D.nc ...
    Comparing  clMISR_2D.nc ...
    Comparing  clcalipso2_2D.nc ...
    Comparing  clcalipso_2D.nc ...
    Comparing  clhcalipso_2D.nc ...
    Comparing  clisccp2_2D.nc ...
    Comparing  cllcalipso_2D.nc ...
    Comparing  clmcalipso_2D.nc ...
    Comparing  cltcalipso_2D.nc ...
    Comparing  cltlidarradar_2D.nc ...
    Comparing  ctpisccp_2D.nc ...
    Comparing  dbze94_2D.nc ...
    Comparing  meantbclrisccp_2D.nc ...
    Comparing  meantbisccp_2D.nc ...
    DIFFER : VALUES : RECORD : 1 : RECORD VARIABLE : meantbisccp
    Comparing  parasol_refl_2D.nc ...
    Comparing  tauisccp_2D.nc ...
    Comparing  tclisccp_2D.nc ...
    DIFFER : VALUES : RECORD : 1 : RECORD VARIABLE : tclisccp

This is telling you that there are differences at bit-level. Sometimes these differences are 
caused by small differences in floating-point arithmetics when using different compilers,
and are usually negligible. Therefore, if you find differences it is advised to 
run cosp_test.sh again, now with arg2=2, to see if those differences are significant.
cosp_test.sh will show the differences that are greater than the fourth signifcant digit.
You can adapt this by changing the variable PRECISION within cosp_test.sh.

Probably, a better option is to run cosp_test.sh with arg2=3. This computes the differences 
between the new and test files and store them in ~/cosp/cosp.vN.N/temp so that you can look 
at them with a visual analysis program (e.g IDL, Matlab, VCDAT). If you find that these 
differences are scientifically significant, please report them to me.

If nccmp and nco utilities are not available and cmp is used instead (arg2=0), 
then the output will look something like this:
    [...]
    Comparing  tauisccp_2D.nc ...
    572  60  64
    573  61  71
    575  63  61
    576  61  60
    1742  60  64
    1743  61  71
    1745  63  61
    1746  61  60
    Comparing  tclisccp_2D.nc ...
    572  60  64
    573  61  71
    575  63  61
    576  61  60
    1766  60  64
    1767  61  71
    1769  63  61
    1770  61  60
    2621  77  76
    2622  63 372
    2623  63 341
    2624  65 112
    2630 200  63
    2631   0  63
    2632   1  64
    2634 263 172
    2635  63 341
    2636  64 111
    2686 346 241
    2687 146 107
    2688 150 257
    2689  77  76
    2690  63 372
    2691  63 341
    2692  65 112
    2694 200  63
    2695   0  63
    2696   1  64
    2697  77  76
    [...]

A few lines (~10) of differences as in tauisccp_2D.nc are fine, because the time stamp 
in the 'history' attribute will differ. A long list of differences like in tclisccp_2D.nc
means that the files differ. In such a case you will have to look at them with a visual 
analysis program and see if the differences are scientifically significant.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
5.- CHANGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

5.1 - Changes in v0.2 (not released to the general public)
--------------------------------------------------------------------------------------------
1) Optimisation of initialisation of variables.
Mie tables were reloaded on each call. Variable initialisations moved into CONSTRUCT_GBX,
and GBX strucutre modified as necessary (Code provided by R. Marchand). This also includes
a correct initialisation of the temperature arrays for Ze scaling, which is only done when
not using mie table.
  
2) Minor changes to definition of variables to allow compilation in GNU Fortran 95.

3) Implemented capability of passing effective radius to the radar simulator.

4) Gaseous absorption.
Gaseous absorption is now only calculated once per gridbox. Optional arguments 
have been added to radar_simulator to allow for this.

5) Output on height levels
The standard output is on model levels. A new capability has been added to output the results on 
height levels.

6) Treatment of non-spherical particles in lidar computations
Polynomial coefficients for non spherical particles derived from a composite of
Ray-tracing theory for large particles (e.g. Noel et al., Appl. Opt., 2001)
and FDTD theory for very small particles (Yang et al., JQSRT, 2003).

7) CMOR-compliant outputs
The outputs are now written using the CMOR library, with one variable per file.

8) New outputs have been implemented:
   a) lidar_only_freq_cloud(Npoints,Nlevels): Cloud frequency of 
       occurrence as seen by CALIPSO but not CloudSat
   b) radar_lidar_tcc(Npoints): Total cloud amount as seen by radar or lidar

9) Bug fixes
Pressure was passed in hPa to lidar simulator and lmd_ipsl_stats. Now in Pa.

10) Optional outputs
A new namelist with logical flags that control which variables are written to file
has been added. This includes some flexibility in the memory allocation.
If a particular instrument is not used, the arrays that hold the output 
variables for that instrument are not allocated.

11) Change of order of levels for radar variables
In V0.1, Ze_tot and cfad_ze were swapped in cosp_io.f90, when they were written
onto disk. Now, the swap is done after calling Quickbeam for Ze_tot. This also changes 
the order in cfad_ze.

12) Deleted fr_hydro from type cosp_subgrid.
That array is not used and it would grab a lot of memory.

13) The ISCCP simulator (v3.7) has been added.

5.2 - Changes in v0.3
--------------------------------------------------------------------------------------------
1) Sanity checks of input variables
Several calls to the subroutine cosp_input_check are included in subroutine cosp (cosp.f90)
that check the limits of many of the inputs.

2) Bug fixes
k2 was defined as integer instead of real. Flip of Reff before calling quickbeam so that the
first row is at the top. Assigment of dimensions in cosp_stats.

3) ISSCP simulator outputs are flipped before they are written out. In V0.2, this was done
by the CMOR library.

4) Changes in actsim/lidar_simulator.f90:
* definition q_part changed (now in-cld water content in input of the routine, so
we remove the division by the cloud fraction). Cloud fraction arguments ls_cldf 
and cv_cldf removed from input arguments because of this.
* depolarisation removed (also from the output arguments)
* header: misleading (0 or 1) removed in description of frac_out (because it 
should be the actual output from scops, which is either 0, 1 or 2)
* parasol module added (see below)

5) Redefined bin boundaries in scattering ratio CFAD (in actsim/lmd_ipsl_stats.f90)

6) In cosp_lidar.f90:
* Instructions "Call for large-scale cloud" removed -> sgx%frac_out is used instead.
* Call to lidar_simulator changed (lsca, gbx%cca and depol removed; frac_out
changed in sgx%frac_out)

7) The labels of layered cloud fractions was wrong, now corrected.
(before: low was actually mid, mid was high, high was total, total was low)

8) New PARASOL diagnostics added. This can be switched off by setting ok_parasol=.false. in
actsim/lidar_simulator.f90. PARASOL-like reflectance at 865 nm for 5 solar zenith
angles (SZAs) added to the list of outputs.

5.3 - Changes in v0.4 (not released to the general public)
--------------------------------------------------------------------------------------------
1) Changes to Makefile.
 + -I$(ISCCP_PATH) option added to scops.o that prevented COSP
 to be build in some compilers (e.g. g95).
 + Deleted the use of CDMS_LIB and UDUNITS_LIB (see reason below)
 
2) Change of units in boxptopisccp before is passed to CMOR.
CMOR can be built without CDMS and UDUNITS. This makes the
build process of COSP much easier as UDUNITS and CDAT are not
needed any more. boxptopisccp was the only variable that was using
that capability.

3) New module that includes a MISR simulator.
This module produces a new diagnostic, clMISR. This is the Tau-CTH histogram of
the cloud fraction as calculated by the MISR Simulator. The axes have 7 and 16
bins for tau and cloud top height, respectively. This module and the output
diagnostic can be switched off/on using the logical flags Lmisr_sim and LclMISR, 
in the COSP_OUTPUT namelist.

4) ISCCP simulator changed to version 3.8.
This version of the simulator adds the capability of outputting two new 
diagnostic variables: the grid-box mean of the 10.5 micron 
infrared brightness temperature for all-sky and clear-sky conditions. 
These variables are called "meantb" and "meantbclr".
See ISCCP simulator readme file for more details.

5) Effective radius and mixing ratios treated consistently.
Type COSP_SUBGRID now includes Reff as 4-D array (Npoints,Ncolumns,Nlevels,Nhydro), 
which makes it consistent with mixing the ratios. This will allow for sub-grid 
variability to be included in the future.

6) Possibility of using precipitation mixing ratios.
There is a new logical flag in the COSP_INPUT namelist (use_precipitation_fluxes) 
that controls the use of precipitation fluxes/mixing ratios.

7) Graupel added to precipitating hydrometeor.
Graupel has been included as precipitating hydrometeor, via fl_lsgrpl in the 
netcdf input files, and stored internally in gbx%grpl_ls.
It is treated as large-scale precipitation in the subgrid overlapping distribution.

8) New script, cosp_test.sh, that compares your outputs with the expected outputs.
Please see Section 4 for details.

9) Optimisations
Subroutine pf_to_mr has been optimised slightly.

10) New logical flag, use_reff, added.
This flag controls the use of the effective radius array by the radar simulator.
If true, then that array is used, otherwise the defaults are used. Note that 
this array is always used by the lidar simulator, so it has to be properly filled 
irrespective of the value of use_reff.
        
11) Use of COSP in 'cloud resolving mode' (Ncolumns=1) using the same interface.
In previous versions, if users wanted to use COSP in this mode, they had to bypass 
the module COSP and call directly COSP_SIMULATOR. Now this can be done just by 
using Ncolumns=1. There are a few things you need to be aware when using this mode:
   a) The use of precipitation fluxes is not supported, so the input variables fl_lsrain, 
      fl_lssnow, fl_lsgrpl, fl_ccrain, and fl_ccsnow must contain mixing ratios. The 
      logical flag use_precipitation_fluxes has to be set to .false. or an error 
      message will be issued.
   b) The optical depth (emmissivity) of all clouds must be passed using dtau_s (dem_s).
      This is relevant if the ISCCP or MISR simulators are switched on.
   c) The lidar simulator is currently limited to 4 cloud hydrometeors (2 liquid, 2 ice). 
      All the cloud mixing ratios have to be condensed into 4, and passed by using the
      arrays mr_lsliq, mr_ccliq for liquid and mr_lsice,mr_ccice for ice.

12) Changes in integration the lidar equation
The accumulated optical depth for layer k contained the contribution from layers 
k to nlevels. This caused problems in the cloud identification, like false 
positives in clear areas under medium-thickness clouds. The vertical integration 
of the lidar equation has been modified to get more realistic values of 
ATB (and SR) in presence of optically thick clouds.
As a consequence of this, several changes have been done:
 a) The cloud detection based on the "first fully attenuated layer 
    encountered from above" has been removed. This removes the false 
    cloud detection mentioned above.
 b) The SR threshold value used to detect cloud has been changed from SR=3
    to SR=5 (to allow a better comparison with day-time Calipso data that are
    more noisy than night-time ones).
 c) The values of SR used to built the diagnostic "diagSR" have been changed
    again (but not the number of bins) to cover the variability
    observed and simulated (after the lidar equation modification).

13) New outputs.
The subcolumn cloud configuration produced by SCOPS (frac_out) and the lidar molecular 
backscatter (beta_mol532) have been added to the list of outputs. They are controlled 
by the logical flags Lfrac_out and Lbeta_mol532, respectively. These variables can 
be useful for debugging purposes, but are not recommended as standard output.

14) Rank of beta_mol reduced to save memory. The variable that holds beta_mol had dimensions 
(Npoints,Ncolumns,Nlevels). The second dimension is redundant as it is a gridbox-mean 
quantity. That dimension has been deleted and the statistical routines modified accordingly.

15) Bug fixes when USE_VGRID=.true. Changes in cosp_stats.f90.
  a) Use of correct dimension for ph_in(Npoints,1,Nlevels) in the 
     allocate statement.
  b) Calls to cosp_change_vertical_grid protected by if statements.
  c) Fix treatment of log_units in cosp_change_vertical_grid.

16) Added default grid when for CSAT_VGRID=.false. (and USE_VGRID=.true.).
In previous versions, the user was responsible for filling the structure 
vgrid (type (cosp_vgrid)) if the standar CloudSat grid was not used.
Now, a default vertical grid with constant step between 0 and 20 km is created.
The grid spacing is controlled by the namelist variable NLR.

17) Section 7 has been added to this README file.

18) Changes required to compile and run with NAG f95. 
  a) Inclusion of KIND statements in Quickbeam for portability.
     KR8=selected_real_kind(15,300) statements have been added 
     in several routines to allow compilation in f95 (NAG compiler), 
     instead of the  _8 format, which was not portable. This affects 
     format_input.f90, math_lib.f90, radar_simulator.f90. It would be 
     desirable to put this in an external module to avoid repetition 
     of the same statement.
  b) The 'intent' of the structures that carry output arguments have 
     been changed to inout. This seems to avoid automatic deallocation 
     that caused runtime errors.

19) The input and output namelists have been split in two different files.

20) Corrected bug in when flipping y%fq_isccp in cosp_isccp_simulator.f90.

21) Minor changes/bugs in the lidar simulator.
A bug found in the conversion of the parasol radiance to reflectance has
been corrected. The code has been cleaned, and depolarization is no 
longer available.

22) The COSP_INPUT and COSP_OUTPUT namelists are now in separate files:
cosp_input_nl.txt, and cosp_output_nl.txt.

5.4 - Changes in v1.0
--------------------------------------------------------------------------------------------
1) Corrected bug in cosp.f90. Reff for precipitation were filled in incorrectly.

2) Default effective radius (30 micron) for lidar simulations when use_reff==.false. and
Reff is set to 0.

3) nf90_close() used instead of ncclos() in cosp_io.f90.

4) S_cld = 5.0 in llnl_stats.f90, consistent with changes in lmd_ipsl_stats.f90.

5) Half_levels should be below full_levels.

6) ISCCP simulator version 4.0.

7) Added isccp_emsfc_lw as input argument to construct_cosp_gridbox.

8) cdms target added to Makefile. This links cosp_test with cdms and udunits, for 
those users that have a complete instalation of the CMOR library.

9) Correction of a minor bug in MISR simulator.

10) Deleted array beta_hydro from type cosp_sglidar as it was not used.

11) Added new capability of calling cosp_simulator iteratively,
with the number of points passed in in each iteration
controlled by the variable NPOINTS_IT in namelist COSP_INPUT. This has involved 
significant changes in cosp.f90, as well as the following changes:
- Deleted Reff and mr_hydro from cosp_subgrid.
- Deleted Nhydro as input argument of construct_cosp_subgrid.

12) Section 2 of this README file changed to reflect the changes in the license terms.

13) Minor change in prec_scops.f to make sure there will be convective precipitation 
when the number of columns is set to a small value (less than 20).

14) ISCCP total cloud sometimes produce values slightly greater than 1 (this is probably 
compiler-dependent). Those values are set to 1.0 in cosp_isccp_simulator.f90.

15) The look-up tables rlumA and rlumB in subroutine parasol have been updated with new values.

16) Changes to MIP tables: the attributes "axis: 1" have been removed to comply with CF conventions.

17) Minor bug fix in computation of pmol and pnorm of upper layer in lidar simulator.

5.5 - Changes in v1.1
--------------------------------------------------------------------------------------------
1) Deleted 'include netcdf.inc' from cosp_io.f90.

2) Bug fixes for the lidar simulator :
     - Computation of pmol and pnorm, thanks to Masaki Satoh: a factor 2 
       was missing. This affects the ATB values but not the cloud fraction. 

3)  Bug fixed in llnl/pf_to_mr.f. In the 5 lines defining  the output fields mx_rain_ls, 
mx_snow_ls, etc., the division by rho is placed beyond column 72, which 
means in this f77 file that it is ignored. These lines have been broken with a continuation
lines. Thanks to Richard Hemler for reporting this.

4) Section 6.3 added to this README on running COSP in Cloud Resolving Mode.

5) Section 6.4 added to this README on the configuration for CFMIP-2.

5.6 - Changes in v1.2
--------------------------------------------------------------------------------------------
1) HCLASS_CP set to 0 in cosp_constants.f90. Comment added to note that this is not used by COSP.

2) Added call to construct_cosp_sghydro in cosp.f90. This only affects when NCOLUMNS=1.

3) cosp.f90: The error exit on condition use_reff=.true. and Reff==0.0 has been deleted. This 
 caused trouble when calling COSP iteratively passing a small number of profiles. In this 
 case, there may be iterations with no cloud in all profiles, causing COSP to stop.

4) The timing messages from cosp_simulator.f90 have been commented out to avoid too much
information get printed to the screen. Users who need to know the timing can uncomment them.

5) lmd_ipsl_stats.f90: Warning message regarding PARASOL being valid only over ocean deleted.
This avoids too many warning messages get printed to the screen, specially for those running
in CRM mode.

6) Subroutine cosp_lidar_only_cloud (llnl_stats.f90): Variables lidar_only_freq_cloud and tcc
were not initialised. This produced the wrong results when calling COSP iteratively without 
reinitialising the output structures.

7) Values of SZA in PARASOL reflectance were indices. Now correct values in degrees (0,20,40,60,80).

8) cosp_stats.f90: in subroutine COSP_CHANGE_VERTICAL_GRID, r is now initalised to R_UNDEF. This avoids
a  Floating Point Exception (divide by zero).

9) New MODIS simulator. The simulator is controlled with the logical flag Lmodis_sim in cosp_output_nl.txt.
Fifteen new variables have been added to the list of outputs.

10) Cloud fractions are now in %. The change of units is done in cosp_simulator.f90.

11) The allocation and deallocation of mt_ttl and mt_tti has been deleted. These global arrays have been 
made of fixed shape (non-allocatable) to avoid problems when running in shared memory mode. This implies that
now there is no possibility of running Quickbeam with Mie tables (this option was not used before, anyway).

12) Possibility of running RTTOV within COSP. Please see section 6.6 below.

13) Changes to comply with CMIP5 experiments. The MIP table COSP_table_1D has been replaced by CMIP5_cf3hr. This table
contains the variable definitions of one of the latest revisions of CMOR2, plus additional variables for the simulators
that are not requested by CMIP5. The main changes between these two tables are:
   - Some of the variables have been renamed.
   - Coordinate height_mlev changed to alevel, and units defined as 1.
   - scat_ratio renamed to scatratio.
   - Last value in ISCCP_TAU changed from 50000.0 to 100.0.
   - Coordinate pressure2 renamed to plev7.
   - Coordinate sza renamed to sza5.
   - MISR_CTH_BNDS and MISR_CTH converted from km to m.
   - Coordinate cloud_top_height renamed to cth16.
   - Coordinate time renamed to time1.
   - START_DATE added to CMOR namelist (cosp_cmor_nl.txt).
   - cloud fractions changed to %.
   
14) The namelists in cfmip2 directory have been updated with the new variable names.

5.7 - Changes in v1.2.1
--------------------------------------------------------------------------------------------
1) Some of the file extensions have changed to .F90 so the they are preprocessed automatically.

2) Bug fix to avoid CALIPSO cloud fractions of zero set to undefined. Changes applied to 
cosp_simulator.f90. Thanks to J. Cole for reporting this bug.

3) Three new diagnostics added to MODIS: low-, mid-, and high-level cloud fractions.

5.8 - Changes in v1.2.2
--------------------------------------------------------------------------------------------
1) Changes to the CloudSat and CALIPSO CFADs.
When an evenly-spaced vertical grid was used (USE_VGRID=.true), 
then the CFAD values at the levels in the new grid that fell below the surface were set to 0.
This is not correct, and now they are set to R_UNDEF, to distinguish them from those bins whose 
value is actually 0. This impacted many land points, basically all with with altitude greater
than 480 m. This was particularly relevant when computing the average CFAD over a region that
contained both land and sea points. The frecuency of occurrence at lower levels would be
underestimated due to the contribution of the zeroes from the land regions.

2) Bug-fix in cosp_stats.f90. sglidar%beta_mol is copied within the Llidar_sim if statement to
avoid a mismatch in the shape of the arrays when Llidar_sim=.false..

3) Bounds added to time coordinate.

4) ISCCP_TOPHEIGHT_DIRECTION set to 2 in the COSP_INPUT namelist. This is the default value since
 V4.0 of the ISCCP simulator. Thanks to R. Heimler for reporting this. This has been updated in the 
 namelists in ./cfmip2.

5) nsizes changed from integer*4 to integer in quickbeam/dsd.f90. Thanks to J. Cole for reporting this.

6) cosp_test.sh: two additional arguments added. Test and output directories are now 
requested as arguments.

7) cosp_constants.f90: DBZE_MAX has been increased from 30.0 to 80.0. This may introduce soma changes in your
  CloudSat CFADs if you model produces reflectivities grater than 30 (very unlikely).

8) The directory 'outputs_test' has been replaced by two directories, 'outputs_test.mlev' and
 'outputs_test.slev'. The .mlev directory contains the test outputs for model level outputs, and 
 the .slev for the standard 40 levels outputs.

9) MODIS simulator: Minor bug fixes. Retrieved 
CloudTopPressure, retrievedSize and retrievedTau set
 to R_UNDEF (rather than 0.0) in subroutine modis_L2_simulator_twoTaus when cloudMask is 
.false..
 This would be consistent with the way ISCCP handles the equivalent fields, and makes it apparent in the plots 
 where MODIS is seeing clouds.
 
10) USE_VGRID=.true. use as default in COSP_INPUT namelist.

5.9 - Changes in v1.3
--------------------------------------------------------------------------------------------
1) MODIS simulator: LWP_conversion multiplied by 1000. This affects iwpmodis and lwpmodis, which are now output in the
  correct units (kg m-2).

2) Bug fix in NC_READ_INPUT_FILE: surface emmisivity was not read properly. This impacted some of the ISCCP and MODIS
  diagnostics. This affects the following variables: albisccp, boxptopisccp, boxtauisccp, clisccp, cllmodis, clmmodis, clmodis, cltisccp, meantbclrisccp, meantbisccp, pctisccp, pctmodis, tauisccp. This only affects results that were 
  created with the off-line version.

3) In NC_READ_INPUT_FILE: added error trapping in all calls to NetCDF routines.

4) sunlit was set always to 1 in cosp_test.F90. This has been fixed and it is now properly filled from the input file.

5) ISCCP simulator upgraded to version 4.1. This fixes a bug where the cloud temperature is greater than any temperature in the troposphere.

6) Bug fix in NC_WRITE_COSP_1D/2D. Trim input arguments to cmor_setup and cmor_dataset. This solves problems with CMOR2rc9/10.

6) Bux fix in MODIS. Small bug in the way the MODIS simulator included in COSP calculates cloud top pressure, in function cloud_top_pressure. Also, a slightly increase efficiency by changing how Cloud_Top_Pressure_Total_Mean is computed in modis_L3_simulator.

7) RTTOV instructions have been updated, in Section 6.6 of this Makefile.

8) New global attribute cosp_version in output files generated with CMOR2.

9) The vertical interpolation routine COSP_CHANGE_VERTICAL_GRID has been optimised for NEC SXs (thanks, Tokuta!). Please see how to switch them on in Section 6.8.

10) Tabs have been replaced with spaces to avoid issues with some compilers.

11) Most of the unused variables have been deleted.

12) rttov_multprof has been put into a module, cosp_rttov.F90.

13) Bug fix in READ_COSP_OUTPUT_NL. The logical flags for the radar&lidar joint diagnostics were not correctly set. They both were set to .true. whenever one of them was requested (and the radar and lidar simulators were on, as needed).

14) If statement added in NC_WRITE_COSP_1D to avoid calling cmor_axis for axes that are not defined in CMIP5_cf3hr. This is not very elegant but works for the CFMIP2 offline experiments.

15) New variables BRANCH_TIME, PARENT_EXPERIMENT_ID, FORCING, INSTITUTE_ID added to CMOR namelist, as requested by the CMOR2 MIP tables for CMIP5.

5.10 - Changes in v1.3.1
--------------------------------------------------------------------------------------------
1) Changes in CMOR output to produce data compliant with the CMIP5 table cf3hr. Latitude and longitude are now
   coordinate variables in each file. New 1D variable produced, toffset.
2) Changes in cosp_rttov to make it compatible with integers promoted to 64-bit. This involves some type-casting.
3) cosp.F90. Added comments that allow to set the seed the way it was done in the ISCCP simulator.
4) The effective radius is now computed by cosp_precip_mxratio, if not supplied as input. This fixes an inconsistency
   when Reff=0.0 and use_precipitation_fluxes=.true. Two new data statements have been added to cosp that define
   gamma_3 and gamma_4, as explained in Section 3 of the user's manual
5) New capability for processing multiple time steps. Please see Section 6 of the user's manual.
6) New text file that gives information on installation in a Mac. See section 6.9 of this README.txt. Thanks to Mike Bauer.
7) New optional variable (verbosity) added to NC_READ_INPUT_FILE. Not used by default. Set this to 1 in the call to 
   NC_READ_INPUT_FILE in cosp_test if you want more information on variables read.
8) Bug fix. Initialise values of axis ids that may not be used in NC_READ_INPUT_FILE.
9) Bug fix. Change of lidar CFAD scattering ratio bounds.
10) New CMOR global attributes: parent_experiment_rip, initialization_method, and physics_version. They are set
    in CMOR namelist.
11) Several standard names have been added to table COSP_table_2D.
12) COSP_table_1D has been updated to keep it parallel to CMIP5_cf3hr.
13) Changes to user's manual.
14) Minor modifications in MODIS simulator to help compilation in some compilers.
15) Extra tau bin (tau<0.3) has been included in MODIS simulator, this makes the MODIS and ISCCP histograms consistent.
    Necessary modifications have been made to the 1D and 2D MIP tables.
16) The CFMIP2 namelists in the cfmip2 directory have been updated to reflect the latest changes.

5.11 - Changes in v1.3.2
--------------------------------------------------------------------------------------------
1) Bug fix in modis_simulator.F90. Runtime error in Lahey compilers when using negative arguments to
  log10. No impact in results as these values were masked out. More details in message sent
  to the users group (2011/02/25).
2) Deleted call to flush subroutine in modis_simulator.F90.
3) Minor correction of this README, section 4.2.
4) Changes to cosp_io.F90. If surface emissivity is not present in input file then it is set to 1.0,
   and a warning message is issued.
5) N_MAX_INPUT_FILES increased to 10000 in cosp_test.F90. The value of this constant can be modified
   to suit the user requirements.
6) Bug fix in cosp_stats.F90/COSP_CHANGE_VERTICAL_GRID. dz goes out of bounds in certain configurations.
   More details in message sent to the users group (2011/03/22).

5.11 - Changes in v1.4
--------------------------------------------------------------------------------------------
1) Improvements in radar simulator:
  - New attenuation integration scheme, possibility of using a look up table, partial support of two-moment microphysics.
  - Several optimisations.
  - Optimisations in gaseous absorption.
  - Selection of microphysics via cosp_defs.h.
2) Optimised version of cosp_change_vertical_grid in cosp_stats.

3) New timing variables that estimate the performance of each simulator.

4) New CALIPSO cloud phase diagnostics.
  - A large list of new cloud diagnostics developed by Gregory Cesana (reference: Cesana G. and H. Chepfer (2013), J. Geophys. Res., doi: 10.1002/jgrd.50376.)
5) Deleted support for CMOR1.
6) Restructure of calls to CMOR functions. A new, improved append mode is used in the the off-line version (see changes in cosp_test.F90).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
6.- ADDITIONAL NOTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
6.1 - Using the NAG f95 compiler
--------------------------------------------------------------------------------------------
In addition to the changes described in item section 5.3(18) for v0.4, the NAG compiler does not support
the FLUSH subroutine (it is not f95 standard). If you want to use this compiler, then you'll need 
to comment out the call to this subroutine in icarus.f (ISCCP simulator).

6.2 - Running COSP in Cloud Resolving Mode (Ncolumns=1)
--------------------------------------------------------------------------------------------
In addition what it is said in section 5.3.11, if you want to run COSP in CRM mode, you will 
need to fill in gbx%mr_hydro with the appropriate precipitation mixing ratios in the driver
program. This can be done by changing the following lines in cosp_test.f90

from
        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice
        gbx%rain_ls = fl_lsrain
        gbx%snow_ls = fl_lssnow
        gbx%grpl_ls = fl_lsgrpl
        gbx%rain_cv = fl_ccrain
        gbx%snow_cv = fl_ccsnow

to
        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice
        gbx%mr_hydro(:,:,I_LSRAIN) = fl_lsrain
        gbx%mr_hydro(:,:,I_LSSNOW) = fl_lssnow
        gbx%mr_hydro(:,:,I_LSGRPL) = fl_lsgrpl

6.3 - Configuration of COSP for CFMIP-2
--------------------------------------------------------------------------------------------
The directory ./cfmip2 contains the namelists with the configuration for the CFMIP-2 experiments.
These files are also available on the CFMIP web site.
There are two different configurations:
a) Long time series (*long_inline.txt). This is the configuration for the 30 yr monthly and daily 
means from ISCCP and CALIPSO/PARASOL. These are global gridded data computed from model gridded inputs, 
with the simulators run inline.

b) Short time series (*short_offline.txt). This is the configuration for the 1 yr time series, both
for the curtain outputs and global gridded monthly means from curtain outputs.
Outputs from CloudSat and CALIPSO/PARASOL are requested.

6.4 - Using your own cloud generator/subgrid variability in hydrometeor water contents
--------------------------------------------------------------------------------------------
If your model requires different overlapping assumptions or the water contents, particle sizes, optical thickness and emissivities, are not horizontally homogeneous, then you will have to make modifications to the code. Below is a list 
of things that Jason Cole did to make it work in his model. Below I attach an the e-mail where he explains the steps 
he followed to do this.

1. Commented out calls to SCOPS in cosp.f90.

2. COSP_TYPES
Added new arrays to subgrid data type "COSP_SUBGRID"

      real,dimension(:,:,:),pointer :: mr_liq ! Non-precipitating cloud liquid mixing ratio
      real,dimension(:,:,:),pointer :: mr_ice ! Non-precipitating cloud ice mixing ratio
      real,dimension(:,:,:),pointer :: reff_liq ! Non-precipitating cloud liquid effective radius
      real,dimension(:,:,:),pointer :: reff_ice ! Non-precipitating cloud ice effective radius
      real,dimension(:,:,:),pointer :: dtau_s   ! Stratiform cloud optical thickness at ~0.65 microns
      real,dimension(:,:,:),pointer :: dtau_c   ! Convective cloud optical thickness at ~0.65 microns
      real,dimension(:,:,:),pointer :: dems_s   ! Stratiform cloud emissivity at ~10.5 microns
      real,dimension(:,:,:),pointer :: dems_c   ! Convective cloud emissivity at ~10.5 microns

3. CONSTRUCT_COSP_SUBGRID
Allocate the new arrays
      allocate(y%mr_liq(Npoints,Ncolumns,Nlevels))
      allocate(y%mr_ice(Npoints,Ncolumns,Nlevels))
      allocate(y%reff_liq(Npoints,Ncolumns,Nlevels))
      allocate(y%reff_ice(Npoints,Ncolumns,Nlevels))
      allocate(y%dtau_s(Npoints,Ncolumns,Nlevels))
      allocate(y%dtau_c(Npoints,Ncolumns,Nlevels))
      allocate(y%dems_s(Npoints,Ncolumns,Nlevels))
      allocate(y%dems_c(Npoints,Ncolumns,Nlevels))

4. FREE_COSP_SUBGRID
Make sure to deallocate the arrays
      deallocate(y%prec_frac, y%frac_out, y%mr_liq, y%mr_ice,
     1           y%reff_liq, y%reff_ice, y%dtau_s, y%dtau_c,
     2           y%dems_s, y%dems_c)

5. I have a "driver" code that is based heavily on "cosp_test.f90" which is provided the fields from the GCM which are assigned to the appropriate variables before calling "cosp" and then assigns the fields from cosp into variables used by the GCM.  For example,

      SUBROUTINE COSP_SIM_DRIVER(
     +                           ICCF,       ! TOTAL CLOUD FRACTION 
ISCCP/COSP (OUTPUT)
...
         call construct_cosp_misr(cfg,Npoints,misr)

        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Compute various fields from GCM input
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         print *, 'Generating necessary fields....'
         call create_fields_cosp(gbx,      ! OUTPUT
     +                           sgx,
     +                           CLW_SUB,  ! INPUT
     +                           CIC_SUB,   
     +                           REL_SUB,   
....
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Call simulator
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         print *, 'Calling simulator...'
         call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,
     1             isccp,misr,stradar,stlidar)
...

      END SUBROUTINE COSP_SIM_DRIVER

Within the subroutine "create_fields_cosp" I make sure that the duplicate fields in the structures "sgx" and "sghydro" are consistent.  
For example, sgx%reff_liq and sghydro%Reff(:,:,:,I_LSCLIQ) have the same values, etc..

6. COSP_LIDAR
Instead of passing gbx%Reff(:,:,I_LSCLIQ), etc pass the subgrid-scale values located in sghydro%Reff(:,i,:,I_LSCLIQ)

7. ICARUS
The call to "icarus" are changed slightly so that the subgrid-scale emissivities and optical thicknesses are passed in:

Original:
  call icarus(0,0,gbx%npoints,sunlit,gbx%nlevels,gbx%ncolumns, &
            pfull,phalf,qv,cc,conv,dtau_s,dtau_c, &
            gbx%isccp_top_height,gbx%isccp_top_height_direction, &
            gbx%isccp_overlap,frac_out, &
            
gbx%skt,gbx%isccp_emsfc_lw,at,dem_s,dem_c,y%fq_isccp,y%totalcldarea, &
            y%meanptop,y%meantaucld,y%meanalbedocld, &
            y%meantb,y%meantbclr,y%boxtau,y%boxptop)

New:
      call icarus(0,0,gbx%npoints,sunlit,gbx%nlevels,gbx%ncolumns,
     1 pfull,phalf,qv,cc,conv,sgx%dtau_s,sgx%dtau_c,
     2 gbx%isccp_top_height,gbx%isccp_top_height_direction,
     3 gbx%isccp_overlap,frac_out,
     4 gbx%skt,gbx%isccp_emsfc_lw,at,sgx%dems_s,sgx%dems_c,y%fq_isccp,
     5 y%totalcldarea, y%meanptop,y%meantaucld,y%meanalbedocld,
     6 y%meantb,y%meantbclr,y%boxtau,y%boxptop)

The file icarus_modified.f shows the changes needed in icarus.f.

6.5 - Incorporating RTTOV in COSP
--------------------------------------------------------------------------------------------
RTTOV is a fast radiative transfer code that provides radiances (and optionally the jacobians) for many infrared and passive microwave sensors. As it is not distributed under an open software license, it is not distributed as part of COSP, althoug it  can be included as an option.
Version 9_1 of RTTOV was released in March 2008 and is available to licensed users free of charge. To become
a licensed user of RTTOV v9, first you will need to register to access the specialist science area of the Met Office:
http://www.metoffice.gov.uk/science/creating/register.html

Then, please send a request using the RTTOV-9 Request Form:
http://www.metoffice.com/research/interproj/nwpsaf/request_forms/request_rttov_9.html

Please follow the instructions to build  RTTOV in your system. 

To build RTTOV within COSP, follow these steps:

1) Build the RTTOV library, following RTTOV instructions. Basically, you will need to execute the following commands in rttovXX/src directory 

../build/Makefile.PL
export ARCH=<architecture>
make all


2) Modify the following environment variables in your COSP Makefile to point
to your local build of RTTOV:
RTTOV_PATH     = /data/cr2/hadac/software/rttov
RTTOV_LIB_PATH = $(RTTOV_PATH)/rttov92.$(F90)/lib 
RTTOV_INC_PATH = $(RTTOV_PATH)/rttov92.$(F90)/include 
RTTOV_MOD_PATH = $(RTTOV_PATH)/rttov92.$(F90)/mod 

3) Modify the following environment variables in your COSP Makefile to point to your local build of RTTOV:
RTTOV_PATH     = /data/cr2/hadac/software/rttov
RTTOV_LIB_PATH = $(RTTOV_PATH)/rttovXX.$(F90)/lib RTTOV_INC_PATH = $(RTTOV_PATH)/rttovXX.$(F90)/include
RTTOV_MOD_PATH = $(RTTOV_PATH)/rttovXX.$(F90)/mod 

4) Depending on the RTTOV version you use, you may need to change the name of the RTTOV library in the make file. Search for llrttov9.1 and replace it with the correct version.

5) Uncomment the following line in cosp_defs.h:
#define RTTOV rttov

6) Run 'make rttov'

7) From your COSP directory, make soft links to the spectral files:
ln -s $RTTOV_PATH/rtcoef/*.dat .

We have tested COSP with RTTOV9, so it is discouraged to use any other version of RTTOV.

The logical Lrttov_sim, in cosp_output_nl.txt switches on/off RTTOV within COSP. 
Ltbrttov switches on/off the RTTOV outputs, i.e. the brighness temperatures. 
The COSP input arguments that control which instrument and channels you want to 
simulate are found in cosp_input_nl.txt, below the lines
  !----------------------------------------------------------------------------------
  !-------------- RTTOV inputs
  !----------------------------------------------------------------------------------
The possible values for the arguments Platform, Satellite, Instrument, Nchannels, Channels can
be found in the RTTOV documentation.

NOTE: the current implementation of RTTOV in COSP only computes clear-sky brightness temperatures. All sky support will be
added in the future.

6.6 - Running COSP on vector computers
--------------------------------------------------------------------------------------------
Since v1.3, COSP includes some optimisations for the vector architectures (tested on NEC SXs). These are protected by #ifdef SYS_SX directives. If you want to switch them on, please uncomment the line 

#define SYS_SX sys_sx

in cosp_defs.h.

Thanks to Tokuta Yokohata, Teruyuki Nishimura and Koji Ogochi for providing these optimisations.

6.7 - Installing COSP in Mac computers
--------------------------------------------------------------------------------------------
The file mac_info.txt contains information on the installation process in a Mac computer. This is a step-by-step 
process for a particular installation in one computer, so it may not be general. It contains information on the 
installation of all the libraries required. Thanks to Mike Bauer for providing this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
7.- OBSERVATIONAL DATASETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

7.1 - GCM Oriented Calipso Cloud Product
--------------------------------------------------------------------------------------------
LMD have prepared the Calipso/Parasol observational datasets, which
are fully consistent with COSP outputs. This work has been done in
collaboration with the PIs of the two missions (D. Winker for Calipso and
D. Tanre for Parasol). These observational datasets are available on
http://climserv.ipsl.polytechnique.fr/cfmip-atrain.html.
