##########################################################################################
# Copyright (c) 2017, Regents of the University of Colorado
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
# Python module dependencies:
# numpy, netCDF4, argparse, warnings
#
# History
# Dec 2017 - D. Swales - Original version
#
# Comments:
# The purpose of this program is to compare COSP2 diagnostics.
# By default, errors are defined as relative differences greater than 1e-5. This threshold
# can be adjusted through command line arguments
# 
# Inputs:
# NAME                  TYPE            DESCRIPTION 
# fileRef               string          COSP output file (Reference)
# fileIN                string          COSP output file
# -zeroThresh           float           error threshold(optional)
#
##########################################################################################
import numpy,netCDF4,argparse,warnings

# Get command line arguments
parser = argparse.ArgumentParser(description='Read user specifications')
parser.add_argument('fileRef',metavar='COSP2_output_file_reference',type=str,nargs=1,
                    help='COSP2 output file (reference)')
parser.add_argument('fileIN',metavar='COSP2_output_file',type=str,nargs=1,
                    help='COSP2 output file')
parser.add_argument('-zeroThresh',nargs='?',type=float,
                    default=0.00001,help='Error tolerance threshold')
args    = parser.parse_args()
fileRef = ''.join(args.fileRef)
fileIN  = ''.join(args.fileIN)

# Supress Python warnings
warnings.filterwarnings("ignore")

print "############################################################################################"
print 'Treating relative differences less than {0:.10f}% as insignificant'.format(100*args.zeroThresh)

# Extract variable names from netCDF files.
dset    = netCDF4.Dataset(fileIN)
varlist = dset.variables.keys()   # Variable names
nvars   = len(varlist)            # Number of variables

# Loop over each variable name and perform comparison.
fracError = [0]*nvars             # Error fraction
minError  = [0]*nvars             # Minimum difference
maxError  = [0]*nvars             # Maximum difference
count = 0

for ij in range(0,nvars):

    # Read in COSP reference data
    fileID_ref = netCDF4.Dataset(fileRef,'r')
    var_ref    = fileID_ref.variables[varlist[ij]][:]

    # Read in COSP data
    fileID = netCDF4.Dataset(fileIN,'r')
    var    = fileID.variables[varlist[ij]][:]

    # Compute relative difference.
    diff  = var/var_ref-1
    diff[numpy.where(var_ref == 0.0)] = 0.0

    # Flag any point which exceeds error threshold. 
    error = numpy.argwhere(abs(diff) >= args.zeroThresh)

    # Error fraction
    fracError[ij] = float(len(error))/diff.size

    # If there are differences, gather range for differences.
    if (fracError[ij] > 0):
        # NOTE: The "argwhere" function only works on 1D arrays, so need to reshape first.
        diff2 = numpy.reshape(diff,diff.size)
        error2 = numpy.argwhere(abs(diff2) >= args.zeroThresh)
        minError[ij] = min(diff2[error2])
        maxError[ij] = max(diff2[error2])
        err  = fracError[ij]*100.00

        # Display differences to screen
        count = count + 1
        print "  "+(varlist[ij]+":").ljust(16)+(" %1.2f" % err).rjust(10),"% of values differ,"+\
            "relative range:".rjust(16)+(" %1.2e " % minError[ij]).rjust(12)+"to"+\
            (" %1.2e" % maxError[ij]).rjust(10)

if (count != 0): print "All other fields match."
if (count == 0): print "All fields match"
# END PROGRAM
