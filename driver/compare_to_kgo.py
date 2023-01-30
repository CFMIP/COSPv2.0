# (c) British Crown Copyright 2020, the Met Office.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#   * Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#   * Neither the name of the Met Office nor the names of its contributors may
#     be used to endorse or promote products derived from this softwarewithout
#     specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

import netCDF4,argparse,sys
import numpy as np

def get_var_list(ncfile):
    """
    Returns a list with the names of the variables in a NetCDF file.
    """
    dataset = netCDF4.Dataset(ncfile)
    var_list = dataset.variables.keys()
    dataset.close()
    return var_list
    

def read_var(fname, vname):
    """
    Reads a variable from a NetCDF file.
    
    Arguments:
        fname: path to NetCDF file.
        vname: variable name.
    """
    f_id = netCDF4.Dataset(fname,'r')
    return f_id.variables[vname][:]

def calculate_stats(tst, kgo, atol=0.0, rtol=None):
    """
    Returns a dictionary with some basic summary statistics of the differences
    between two numpy arrays. When a threshold for relative differences is
    passed, then the results refer to relative differences.
    The dictionary contains the following keys:
        'N': number of differences.
        'AvgDiff': average difference.
        'MinDiff': minimum difference.
        'MaxDiff': maximum difference.
        'StDev': standard deviation of the differences.
    
    Arguments:
        tst: test variable (numpy array).
        kgo: reference variable (numpy array).
    Keyword arguments:
        atol: absolute tolerance threshold. Smaller differences pass the test.
        rtol: relative tolerance threshold. Smaller differences pass the test.
    """
    summary_stats = {'N':0, 'AvgDiff':0.0, 'MinDiff':0.0, 'MaxDiff':0.0, 'StDev':0.0}
    # All differences
    d = tst - kgo
    # Mask for differences larger than absolute tolerance
    maskAllDiff = (np.absolute(d) > atol)
    NallDiff = maskAllDiff.sum()
    # If there are differences larger than atol,
    # then calculate summary statsitics
    if (NallDiff > 0):
        diffs = d[maskAllDiff]
        # Are relative differences requested?
        if rtol is not None:
            # Calculate relative differences. When KGO=0,
            # use test as reference (i.e. relative diff will
            # be 1 or -1)
            maskedKgo = kgo[maskAllDiff]
            rdiffs = diffs / maskedKgo
            rdiffs[maskedKgo == 0.0] = np.sign(diffs[maskedKgo == 0.0])
            # Keep only those diffs larger than relative tolerance
            diffs = rdiffs[np.absolute(rdiffs) > rtol]
            NallDiff = len(diffs)
        # Calculate summary stats
        summary_stats['N'] = NallDiff
        if NallDiff > 0:
            summary_stats['AvgDiff'] = diffs.mean()
            summary_stats['MinDiff'] = diffs.min()
            summary_stats['MaxDiff'] = diffs.max()
            summary_stats['StDev'] = diffs.std()
        
    return summary_stats

def print_stats_table(summary_stats, print_all=False, stats_file=None):
    """
    Print table of summary statistics.
    
    Arguments:
      summary_stats: dictionary with summary statistics. Output from call
                     to function calculate_stats.

    Keywords:
      print_all: by default, it only prints lines with differences,
                 unless print_all==True.
      stats_file: full path to output file.
    """
    # Open file if required
    if stats_file is not None: 
        stdout = sys.stdout
        f_id = open(stats_file, 'w')
        sys.stdout = f_id
    # Header and column names
    print(42*'=', ' Summary statistics ', 42*'=')
    line = ('{:>40s} {:>10s} {:>12s} {:>12s} {:>12s} {:>12s}').format('Variable',
            'N', 'AvgDiff', 'MinDiff', 'MaxDiff', 'StDev')
    print(line)
    # Main table.
    for vname, vstats in summary_stats.items():
        line = ('{vname:>40s} {N:10d} {AvgDiff:12.4e} {MinDiff:12.4e} '
                '{MaxDiff:12.4e} {StDev:12.4e}').format(vname=vname,**vstats)
        if (vstats['N'] > 0):
            print(line)
        else:
            if print_all: print(line)
    print(106*'=')
    # Close file and restore stdout
    if stats_file is not None:
        f_id.close()
        sys.stdout = stdout

#######################
# Main
#######################
if __name__ == '__main__':

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("kgo_file", help="File with known good outputs.")
    parser.add_argument("tst_file", help="Test output file.")
    parser.add_argument("--atol", type=float, default=0.0,
                        help="Absolute tolerance.")
    parser.add_argument("--rtol", type=float, default=None,
                        help="Relative tolerance.")
    parser.add_argument("--allvar", type=bool, default=False,
                        help="Print summary stats for all variables.")
    parser.add_argument("--noerror", type=bool, default=False,
                        help="Exit with no error even if differences are "
                             "bigger than tolerances.")
    parser.add_argument("--stats_file", default=None,
                        help="Output file for summary statistics.")
    args = parser.parse_args()

    # Get list of variables
    kgo_vars = get_var_list(args.kgo_file)
    tst_vars = get_var_list(args.tst_file)
    nkgo = len(kgo_vars)
    ntst = len(tst_vars)

    # Colours for error messages
    red_colour   = '\033[91m'
    green_colour = '\033[92m'
    std_colour   = '\033[0m'

    # Dictionary for summary statistics
    summary_stats = {}

    # Iterate over shortest list and calculate stats
    errored = False
    if (nkgo <= ntst):
        vlst = kgo_vars
    else:
        vlst = tst_vars
    for vname in vlst:
        kgo = read_var(args.kgo_file, vname) # KGO
        tst = read_var(args.tst_file, vname) # test
        summary_stats[vname] = calculate_stats(tst, kgo, 
                                               atol=args.atol, rtol=args.rtol)
        if summary_stats[vname]['N'] > 0: errored = True

    # Print summary stats
    if errored:
        print(red_colour + "===== ERROR: some of the differences are larger "
              "than the tolerances." + std_colour)
    print_stats_table(summary_stats, print_all=args.allvar,
                      stats_file=args.stats_file)

    # Check if large differences should be treated as errors
    if args.noerror: errored = False

    # Error if files have different number variables. If the number 
    # of variables is the same but they have different names, it will
    # fail in summary_stats.
    if (nkgo != ntst):
        errored = True
        print(red_colour +
              "===== ERROR: files differ in the number of variables." +
              std_colour)
        print("===== Variables in KGO: ", nkgo)
        print(kgo_vars)
        print("===== Variables in Test: ", ntst)
        print(tst_vars)

    # Exit with correct error condition
    if errored:
        print(red_colour + "===== ERROR: the test is exiting with an error, "
              "please review the output above." + std_colour)
        sys.exit(1)
    else:
        print(green_colour + "===== Test completed successfully: "
              "all outputs are within the expected tolerances." + std_colour)
        sys.exit()
