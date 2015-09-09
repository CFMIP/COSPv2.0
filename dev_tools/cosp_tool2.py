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
# The purpose of this program is to compare COSP simulator data. The user provides the 
# locations of the datasets to compare and the variables in which to plot.
#
# Coded by Dustin Swales 2014 CIRES/NOAA-ESRL-PSD
#
#                                INPUTS
# VARIABLE NAME         TYPE            DESCRIPTION 
# dirOUT                string          location of COSP simulator output data
# dirREF                string          location of COSP reference data
# dataset               string          either 1D or 2D
# fieldIN               string          list of variables to plot
# -zeroThresh           float           error threshold (optional)
# -logScale             logical         option for plot axis, default is linear (optional) 
#
##########################################################################################
import os,fnmatch,netCDF4,sys,numpy,warnings,math
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MaxNLocator,LogLocator,ScalarFormatter
from pylab import figure, show

# Get command line arguments
parser = argparse.ArgumentParser(description='Read user specifications')
parser.add_argument('dirOUT',metavar='Input Directory 1',type=str,nargs=1,
                    help='Directory of output data')
parser.add_argument('dirREF',metavar='Input Directory 2',type=str,nargs=1,
                    help='Directory of reference data')
parser.add_argument('dataset',metavar='Sample Dataset',type=str,nargs=1,
                    help='Reference dataset for comparison (1D or 2D)')
parser.add_argument('fields',nargs='+',
                    help='Variables to plot')
parser.add_argument('-zeroThresh',nargs='?',type=float,
                    default=0.00001,help='Error tolerance threshold')
parser.add_argument('-logScale',help='Plot axis option flag',action="store_true")

args       = parser.parse_args()
dirOUT     = ''.join(args.dirOUT)
dirREF     = ''.join(args.dirREF)
dataset    = ''.join(args.dataset)

# Supress Python warnings
warnings.filterwarnings("ignore")

print "#######################################################################################"
print "Using an error threshold of:",args.zeroThresh
print "Making plots for the following variable(s):"
for ii in range(0,len(args.fields)): 
    print "   ",args.fields[ii]

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
                result.append(os.path.join(root,name))

    return result

# Create list of input data files and reference data files
nfilesOUT = len(args.fields)
filesOUT  = ['']*nfilesOUT
filesREF  = ['']*nfilesOUT
fileflag = [0]*nfilesOUT
for ifile in range(0,nfilesOUT):
    varS = args.fields[ifile]
    filesOUT[ifile]  = find(varS+'*'+dataset+'*.nc',dirOUT);
    filesREF[ifile]  = find(varS+'*'+dataset+'*.nc',dirREF);
    fileflag[ifile] = 1
    if (''.join(filesOUT[ifile]) == ''):
        fileflag[ifile]=0
        print "The variable you requested, "+varS+\
            ", is not available in the output data directory you provided."
    if (''.join(filesREF[ifile]) == ''):
        fileflag[ifile]=0
        print "The variable you requested, "+varS+\
            ", is not available in the reference data directory you provided."

# HOOK 1) Quit if all variables requested are unavailable
if (sum(fileflag) == 0):
    print "ERROR: All variables you requested are not present."
    print "EXITING!"
    quit()

##############################################################################################
# 2) Loop over all output data files and compare each field within the file to the reference
#    data. Differences are determined by the user specified error threshold.
##############################################################################################

# Loop over each output data file.
valError  = [0]*nfilesOUT
minError  = [0]*nfilesOUT
maxError  = [0]*nfilesOUT
pcount    = 0
for ij in range(0,nfilesOUT):

    # Only perform comparison if there is both an input file and a reference file to compare to.
    if (fileflag[ij] == 1):

        # Open files
        fileID_out = netCDF4.Dataset(''.join(filesOUT[ij]),'r')
        fileID_ref = netCDF4.Dataset(''.join(filesREF[ij]),'r')

        # Read-in output data 
        var_out  = fileID_out.variables[args.fields[ij]][:]
        var_unit = fileID_out.variables[args.fields[ij]].units

        # Read-in reference data
        var_ref = fileID_ref.variables[args.fields[ij]][:]

        # Compare output data to reference data.
        diff         = var_out/var_ref-1
        diff[numpy.isinf(diff)]=1.0
        diff         = numpy.reshape(diff,diff.size)   
        a            = numpy.argwhere(abs(diff) >= args.zeroThresh)        
        valError[ij] = float(len(a))/diff.size

        # Plot variables that contain errors
        if (len(a) > 0):

            # Compute error range
            minError[ij] = min(diff[a])
            maxError[ij] = max(diff[a])

            # Set up figure size and properties
            pcount=pcount+1
#            plt.figure(pcount, figsize=(16,6), dpi=80, facecolor='w', edgecolor='k')
#
#            # 1) Create histograms of error distribution
#            plt.subplot(1,3,1)
#            diff           = abs(diff)
#            histBins       = [1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,\
#                                  1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
#            hist,bins      = numpy.histogram(diff[a],bins=histBins)
#            precntErrOccur = 100*float(len(a))/diff.size
#            widths         = numpy.diff(bins)
#            hist           = 100*hist/float(sum(hist))
#            plt.bar(bins[:-1],hist,widths)
#            plt.xscale("log")
#            plt.xlim([1e-9,1e9])
#            plt.xlabel("Relative Error")
#            plt.ylabel("Frequency of Occurence")
#            plt.title(args.fields[ij]+" Relative Errors")
#   
#            # 2) Create scatter plot of test values and reference values
#            plt.subplot(1,3,2)
#            y = numpy.reshape(var_out,var_out.size)
#            x = numpy.reshape(var_ref,var_ref.size)
#            if (args.logScale):
#                plt.xscale("log")
#                plt.yscale("log")
#            plt.plot(x,y,'*')
#            plt.xlabel("Test Value ("+var_unit+")")
#            plt.ylabel("Reference Value ("+var_unit+")")
#            plt.title(args.fields[ij]+" Comparison")
#                    
#            # 3) Create scatter plot of residuals
#            plt.subplot(1,3,3)
#            if (args.logScale):
#                plt.xscale("log")
#                plt.yscale("log")
#            plt.plot(x,y-x,'*')
#            plt.xlabel("Residual Value ("+var_unit+")")
#            plt.ylabel("Reference Value ("+var_unit+")")
#            plt.title(args.fields[ij]+" Residuals (Test-Reference)")

            # Set up figure size and properties
            fig = figure(pcount, figsize=(16,6), dpi=80, facecolor='w', edgecolor='k')

            # 1) Create histograms of error distribution
            diff           = abs(diff)
            histBins       = [1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,\
                                  1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
            hist,bins      = numpy.histogram(diff[a],bins=histBins)
            precntErrOccur = 100*float(len(a))/diff.size
            widths         = numpy.diff(bins)
            hist           = 100*hist/float(sum(hist))
            ax1 = fig.add_subplot(131)
            ax1.bar(bins[:-1],hist,widths)
            ax1.set_xscale("log")
            ax1.set_xlim([1e-9,1e9])
            ax1.set_xlabel("Relative Error")
            ax1.set_ylabel("Frequency of Occurence")
            ax1.set_title(args.fields[ij]+" Relative Errors")

   
            # 2) Create scatter plot of test values and reference values
            y = numpy.reshape(var_out,var_out.size)
            x = numpy.reshape(var_ref,var_ref.size)
            ax2 = fig.add_subplot(132)
            ax2.xaxis.set_major_locator(MaxNLocator(6))
            ax2.yaxis.set_major_locator(MaxNLocator(6))
            ax2.xaxis.get_major_formatter().set_powerlimits((0, 1))
            ax2.yaxis.get_major_formatter().set_powerlimits((0, 1))
            if (args.logScale):
                ax2.set_xscale('log')
                ax2.set_yscale('log')
                ax2.xaxis.set_major_locator(LogLocator(base = 10.0,numticks=6))
                ax2.yaxis.set_major_locator(LogLocator(base = 10.0,numticks=6))
            ax2.plot(x,y,'*')
            ax2.set_ylabel("Test Value ("+var_unit+")")
            ax2.set_xlabel("Reference Value ("+var_unit+")")
            ax2.set_title(args.fields[ij]+" Comparison")
                    
            # 3) Create scatter plot of residuals
            ax3 = fig.add_subplot(133)
            ax3.xaxis.set_major_locator(MaxNLocator(6))
            ax3.yaxis.set_major_locator(MaxNLocator(6))
            ax3.xaxis.get_major_formatter().set_powerlimits((0, 1))
            ax3.yaxis.get_major_formatter().set_powerlimits((0, 1))
            if (args.logScale):
                ax3.set_xscale('log')
                ax3.set_yscale('log')
                ax3.xaxis.set_major_locator(LogLocator(base = 10.0,numticks=6))
                ax3.yaxis.set_major_locator(LogLocator(base = 10.0,numticks=6))
            ax3.plot(x,y-x,'*')
            ax3.set_ylabel("Residual Value ("+var_unit+")")
            ax3.set_xlabel("Reference Value ("+var_unit+")")
            ax3.set_title(args.fields[ij]+" Residuals (Test-Ref)") 

            # Panel all subplots
            show(block=False)

        if (len(a) == 0):
            fileflag[ij]=0
            print "The variable you requested, "+args.fields[ij]+", does not meet the error threshold."
    
# Print some simple statistics on errors
for ij in range(0,nfilesOUT):
    if (ij == 0): print "Differences exist in reference and test files for:"
    if (fileflag[ij] == 1):
        err  = valError[ij]*100.00
        print "  "+(args.fields[ij]+":").ljust(16)+(" %1.2f" % err).rjust(10),"% of values differ,"+\
            "relative range:".rjust(16)+(" %1.2e " % minError[ij]).rjust(12)+"to"+\
            (" %1.2e" % maxError[ij]).rjust(10)

print "Plots for ",sum(fileflag),"of the",len(args.fields),"requested variables were created."
print "#######################################################################################"
raw_input()

# END PROGRAM
