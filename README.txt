%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Written by Dustin Swales (dustin.swales@noaa.gov) 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

General notes:
vN.N is version number (e.g. v2.0)
We will assume that the software will be installed in ~/cosp.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1) ABOUT THE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The CFMIP Observation Simulator Package (COSP) takes the models representation of the
atmosphere and simulates the retrievals for several passive (ISCCP, MISR and MODIS) and active
(CLUDSAT(radar) and CALIPSO(lidar)) sensors.

COSP Version 2.0 (COSPv2.0) is a major reorganization and modernization of the previous
generation. Some of the major differences between v2.0 and previous versions of COSP (e.g.
1.4.0) are...
*) COSPv2.0 now expects subcolumn inputs. In previous versions, the subclumns were drawn
   internally by COSP and provided as an output. In COSP2, we provide the resources to
   calculate subcolumn optical-inputs as were done in previous versions of COSP, however we
   encourage users to provide COSP with as much information about the host-model as possible.
*) Explicit initialization of static fields
*) Paramaterized working precision.

The simulators in COSP (ISCCP, MISR, MODIS, RADAR (cloudsat) and LIDAR (calipso) have been
developed by many institution and agencies:
*) Met Office Hadley Centre
*) LLNL (Lawrence Livermore National Laboratory)
*) LMD/IPSL (Laboratoire de Meteorologie Dynamique/Institut Pierre Simon Laplace)
*) CSU (Colorado State University)
*) UW (University of Washington)

The logical flow of COSP is as follows:
*) Subcolumn retrievals (all simulators).
*) Column retrievals (all simulators).
*) Joint-instrument products.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2) CONDITIONS OF USE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The code is distributed under BSD License (http://www.opensource.org/licenses/bsd-license.php).
Each source file includes a copy of this license with details on the Owner, Year 
and Organisation. The license in the file quickbeam/README applies to all the files in 
the directory quickbeam.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3) DOWNLOADING AND UNPACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The code is hosted by GitHub (https://github.com/CFMIP). To download the code, simply "clone"
the repository using the following command:

   git clone https://github.com/CFMIP/COSPv2.0.git

This will create a local copy, by default COSPv2.0/, of the COSP source code on your machine.

Legacy versions of COSP (v1.3.1, v1.3.2, v1.4.0 and v1.4.1) are also now archived on GitHub. These
4 versions are archived as different "commits" in the same repository. To download one of these
versions you will need to clone the COSPv1 repository and then revert to the correct commit.

For example, to download COSP v1.3.2 you would first clone the COSPv1 repository:

   git clone https://github.com/CFMIP/COSPv1.git

This will create a local copy of the repository, by default COSPv1/. Then to view the different
"commits" in that repository you can use "git log". For example:

   git log
     commit bae896c6d09af4f34493422974cfb86722f9ead5
     Author: Dustin Swales <dustin.swales@noaa.gov>
     Date:   Fri Nov 6 08:26:35 2015 -0700

          COSP version 1.4.1

     commit 572c53ab941c8019135d147fbdf01a6633487aaa
     Author: Dustin Swales <dustin.swales@noaa.gov>
     Date:   Mon Sep 14 11:53:52 2015 -0600

          COSP version 1.4.0

     commit 9e7d84d735249d95521852e0db0c8eed16b4f070
     Author: Dustin Swales <dustin.swales@noaa.gov>
     Date:   Mon Sep 14 11:15:51 2015 -0600

         COSP version 1.3.2

     commit ce4130ca7334bc30875e7cfc620d6195bf237c73
     Author: Robert.Pincus <Robert.Pincus@colorado.edu>
     Date:   Fri Mar 18 20:52:12 2011 +0000

         Copying v1.3.1 to archive of stable releases.   

This shows the commit IDs, author of the commit, date of the commit and a brief description
for the three archived releases. By default in git, when you clone a repository, you will be
at the latest commit. In this case COSPv1.4.1. To revert to an older
commit you need to "checkout" that commit (i.e. git checkout <commit>). So for example, to
revert to COSP v1.3.2, you would:

   git checkout 9e7d84d735249d95521852e0db0c8eed16b4f070

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4) COMPILATION AND TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a) Compilation

   i) Compiling COSP2 for in-line use (e.g. within GCM)
      The main COSP source code is located in src/. The Makefile will need to be
      modified to fit your systems architecture.
        cd src/
        make clean
        make install
      This will create all of the object files and modules and place them in a new directory, obj/.
      
   ii) Compiling COSP2 for off-line use.
      For offline testing, a driver is provided with COSP2. This driver can be built with the
      capability of creating CMOR compliant netcdf output, or output to a single netCDF file.
      These drivers both require the netCDF library (http://www.unidata.ucar.edu/software/netcdf/).

      *) non-CMOR compliant netCDF output
        cd driver/src
        make cosp     ! Build COSP source code
        make          ! Build driver for COSP
	
      *) CMOR2
        This driver requires the CMOR2 libraries. (www2-pcmdi.llnl.gov/cmor/download/)
        cd driver/src/cmor
        make cosp     ! Build COSP source code
        make          ! Build driver for COSP

b) Running
   i)  Set up COSP input (driver/src/cosp2_input_nl.txt) and output (driver/src/cosp2_output_nl.txt)
       namelists. The input namelist controls the COSP setup (i.e. Numebr of subcolumns to be
       used, etc...) and simulator specific information (i.e. Radar simulator frequency). The
       output namelist contains a list of logicals, one for each COSP diagnostic output field.
	  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
5) CHANGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
6) NOTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
