##########################################################################################
# Copyright (c) 2015, Regents of the University of Colorado
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are 
# permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of 
#    conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
#    of conditions and the following disclaimer in the documentation and/or other 
#    materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be 
#    used to endorse or promote products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
# THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# History
# May 2015 - D. Swales - Original version
# 
##########################################################################################
#
# The purpose of this program is to compare COSP simulator data.
#
# Coded by Dustin Swales 2014 CIRES/NOAA-ESRL-PSD
#
#                                           INPUTS
# VARIABLE NAME         TYPE            DESCRIPTION 
# dirOUT                string          location of COSP simulator output data
# dirREF                string          location of COSP reference data
# dataset               string          either 1D or 2D
# -zeroThresh           float           error threshold(optional)
#
##########################################################################################
import os,fnmatch,numpy,netCDF4,sys,argparse,warnings

# Get command line arguments
parser = argparse.ArgumentParser(description='Read user specifications')
parser.add_argument('dirOUT',metavar='Input Directory 1',type=str,nargs=1,
                    help='Directory of output data')
parser.add_argument('dirREF',metavar='Input Directory 2',type=str,nargs=1,
                    help='Directory of reference data')
parser.add_argument('dataset',metavar='Sample Dataset',type=str,nargs=1,
                    help='Reference dataset for comparison (1D or 2D)')
parser.add_argument('-zeroThresh',nargs='?',type=float,
                    default=0.00001,help='Error tolerance threshold')
args       = parser.parse_args()
dirOUT     = ''.join(args.dirOUT)
dirREF     = ''.join(args.dirREF)
dataset    = ''.join(args.dataset)

# Supress Python warnings
warnings.filterwarnings("ignore")

print "############################################################################################"
print 'Treating relative differences less than {0:.10f}% as insignificant'.format(100*args.zeroThresh)
##############################################################################################
# 1) Find output files from user supplied directory and corresponding reference files
#    for comparison. If the directories provided by the users do not contain data, the code
#    will quit
##############################################################################################
# Simple directory searching function.
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))

    return result

# Create list of the output data file names from the output directory
filesOUT  = find('*'+dataset+'*.nc',dirOUT);
nfilesOUT = len(filesOUT);

# Create list of the reference data file names from the reference directory
fileflag = [0]*nfilesOUT;       # 1 for correcponding file / -1 for no corresponding refernce file
filesREF = ['']*nfilesOUT
noref    = [0]*nfilesOUT
for ij in range(0,nfilesOUT):
    test1 = filesOUT[ij];
    test2 = find(test1[len(dirOUT)+1:len(test1)-55]+'_*'+dataset+'*.nc',dirREF);
    if (test2 == []):
        print "Variables in output diretory, but missing from reference directory"
        print test1
        fileflag[ij] = 0
        filesREF[ij] = ''
    else:
        fileflag[ij] = 1
        filesREF[ij] = ''.join(test2)#dirREF+"/"+test1[len(dirOUT)+1:len(test1)]

# HOOK 1) User provided reference data directory does not contain data
if (sum(fileflag) == 0):
    print "ERROR: User provided reference data directory is empty"
    print "EXITING!!!"
    quit()

# HOOK 2) User provided output data directory does not contain any data
if (nfilesOUT == 0):
    print "ERROR: User provided output data directory is empty"
    print "EXITING!!!"
    quit()

##############################################################################################
# 2) Loop over all output data files and compare each field within the file to the reference
#    data. Differences are determined by the numerical spacing in the reference dataset.
##############################################################################################

# Loop over each output data file.
fileError  = [0]*nfilesOUT
valError   = [0]*nfilesOUT
minError   = [0]*nfilesOUT
maxError   = [0]*nfilesOUT
varList    = ['']*nfilesOUT
for ij in range(0,nfilesOUT):

    # Only perform comparison if there is reference file to compare to
    if (fileflag[ij] == 1):

        # Open files
        fileID_out = netCDF4.Dataset(filesOUT[ij],'r')
        fileID_ref = netCDF4.Dataset(filesREF[ij],'r')
 
        # Pull out name of variable to read-in from filename
        test1       = filesOUT[ij]
        varList[ij] = test1[len(dirOUT)+1:len(test1)-55]
        if (ij == 0): print "Comparing variables"
        print varList[ij]

        # Read-in output data 
        var_out = fileID_out.variables[varList[ij]][:]

        # Read-in reference data
        var_ref = fileID_ref.variables[varList[ij]][:]

        # Compare output data to reference data.
        diff  = var_out/var_ref-1
        diff2 = var_ref-var_out   
        diff[numpy.isinf(diff)] = 1.0
        diff[numpy.isnan(diff)] = 0.0
        diff[numpy.where(diff == -1)]=0.0
        diff[numpy.where(diff2 == 0)]=0.0
        noref[ij]    = float(len(numpy.argwhere(var_ref/var_out == 0.0)))/var_ref.size
        diff         = numpy.reshape(diff,diff.size)        
        a            = numpy.argwhere(abs(diff) >= args.zeroThresh)
            
        valError[ij] = float(len(a))/diff.size
        if (len(a) > 0):
            fileError[ij] = 1
            minError[ij] = min(diff[a])
            maxError[ij] = max(diff[a])

# Print out where differences exist
if (sum(fileError) >= 1):
    for ij in range(0,nfilesOUT):
        if (ij == 0): print "Differences exist in reference and test files for:"
        if (fileError[ij] == 1):
            err  = valError[ij]*100.00
            print "  "+(varList[ij]+":").ljust(16)+(" %1.2f" % err).rjust(10),"% of values differ,"+\
                "relative range:".rjust(16)+(" %1.2e " % minError[ij]).rjust(12)+"to"+\
                (" %1.2e" % maxError[ij]).rjust(10)
            if (ij == nfilesOUT-1): print "All other files match."

# Print out if test data is non zero and reference data is zero
if (sum(noref) > 1):
    for ij in range(0,nfilesOUT):
        if (ij == 0): print "Test data exists, while reference data is missing for the following:"
        if (noref[ij] > 0):
            err  = valError[ij]*100.00
            print "  "+(varList[ij]+":").ljust(16)+(" %1.2f" % noref[ij]).rjust(10),"% of the time."
            if (ij == nfilesOUT-1): print "All other files match."


# Print summary for zero errors.
if (sum(fileError) == 0):print "All Files Match"

print "############################################################################################"

# END PROGRAM
