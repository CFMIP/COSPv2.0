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
    between two numpy arrays. The dictionary contains the following keys:
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
            diffs = diffs[np.absolute(rdiffs) > rtol]
            NallDiff = len(diffs)
        # Calculate summary stats
        summary_stats['N'] = NallDiff
        if NallDiff > 0:
            summary_stats['AvgDiff'] = diffs.mean()
            summary_stats['MinDiff'] = diffs.min()
            summary_stats['MaxDiff'] = diffs.max()
            summary_stats['StDev'] = diffs.std()
        
    return summary_stats

def print_stats_table(summary_stats):
    """
    Print table of summary statistics.
    """
    for key, s in summary_stats.items():
        print(key,s)

#######################
# Main
#######################
if __name__ == '__main__':

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("kgo_file", help="File with known good outputs.")
    parser.add_argument("tst_file", help="Test output file.")
    parser.add_argument("--atol",type=float,
                        default=0.0,help="Absolute tolerance.")
    parser.add_argument("--rtol",type=float,
                        default=None,help="Relative tolerance.")
    args = parser.parse_args()

    # Get list of variables
    kgo_vars = get_var_list(args.kgo_file)
    tst_vars = get_var_list(args.tst_file)
    nkgo = len(kgo_vars)
    ntst = len(tst_vars)

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
    print_stats_table(summary_stats)
    
    # Error if files have different number variables. If the number 
    # of variables is the same but they have different names, it will
    # fail in summary_stats.
    if (nkgo != ntst):
        errored = True
        print("=== Variables in KGO ===")
        print(kgo_vars)
        print("=== Variables in Test ===")
        print(tst_vars)

    # Exit with correct error condition
    if errored:
        sys.exit(1)
    else:
        sys.exit()
