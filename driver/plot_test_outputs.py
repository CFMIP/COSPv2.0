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
#     be used to endorse or promote products derived from this software without
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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os

def collapse_dimensions_for_plotting(longitude, latitude, vname, vx, vd, dims):
    """
    Pre-processing of COSP variable for plotting.

    Arguments:
        longitude: longitude from COSP output file [loc].
        latitude: latitude from COSP output file [loc].
        vname: variable name.
        vx: variable from COSP output file [..., loc]
        vd: dictionary with metadata about variable.
        dims: dictionary with additional dimensions.
    Return:
        x: x axis values.
        y: y axis values.
        z: data array for plotting.
        d: dictionary with plot configuration.
    """
    yflip = False
    xticks, yticks, xticks_labels, yticks_labels, xlabel, ylabel, vmax = (None,)*7
    if vd['plot_type'] == 'map':
        d = None
        x = longitude[0]
        y = latitude[:,0]
        z = vx
        if vname == 'parasolGrid_refl': z = vx[2]
        # Roll longitude if there are values > 180
        m = (x > 180.0)
        if m is not None:
            Nroll = longitude.shape[1] // 2
            x[m] = x[m] - 360.0
            x = np.roll(x, Nroll)
            z = np.roll(z,Nroll,axis=1)
        # Calculate latitude and longitude edge points.
        # Assume they are increasing monotonically.
        # Extend length to N+2 and calculate midpoints.
        x = midpoints_to_edges(x)
        y = midpoints_to_edges(y)
        xticks = np.arange(-180,181,60)
        yticks = np.arange(-90,91,30)
        xlabel = 'Longitude (deg)'
        ylabel = 'Latitude (deg)'
    if vd['plot_type'] == '2Dhist':
        weights = np.cos(latitude * np.pi / 180.0)
        weights = weights / weights.sum()
        z = np.sum(vx * weights, axis=2)
        x = np.arange(z.shape[1]+1)
        y = np.arange(z.shape[0]+1)
    if vd['plot_type'] == 'zonal_cross_section':
        z = np.average(vx, axis=2)
        x = midpoints_to_edges(latitude[:,0])
        y = np.arange(z.shape[0] + 1)
    if vd['plot_type'] in ('2Dhist','zonal_cross_section'):
        if vd['xaxis_type'] == 'tau7':
            xticks_labels = ('0', '0.3', '1.3', '3.6', '9.4', '23', '60', '')
            xticks = x
            xlabel = 'Cloud optical depth'
        if vd['xaxis_type'] == 'cloudsat_DBZE_BINS':
            x = np.arange(-50,26,5)
            xticks = x
            xticks_labels = None
            xlabel = 'Radar reflectivity (dBZ)'
        if vd['xaxis_type'] == 'SR_BINS':
            xticks_labels = ('0', '0.01', '1.2', '3', '5', '7', '10', '15', '20',
                             '25', '30', '40', '50', '60', '80', '')
            xticks = x
            xlabel = 'Lidar scattering ratio'
        if vd['xaxis_type'] == 'latitude':
            xticks_labels = None
            xticks = np.arange(-90,91,30)
            xlabel = 'Latitude (deg)'
        if vd['yaxis_type'] == 'pres7':
            yticks_labels = ('1000', '800', '680', '560', '440', '310', '180','')
            yticks = y
            ylabel = 'Cloud Top Pressure (hPa)'
        if vd['yaxis_type'] == 'hgt16':
            yticks_labels =  ('', '0', '500', '1000', '1500', '2000', '2500',
                              '3000', '4000', '5000', '7000', '9000', '11000',
                              '13000', '15000', '17000', '')
            yticks = y
            ylabel = 'Cloud Top Height (m)'
        if vd['yaxis_type'] == 'REICE_MODIS':
            yticks_labels =  ('0', '10', '20', '30', '40', '60', '90')
            yticks = y
            ylabel = 'Ice particle size (micron)'
        if vd['yaxis_type'] == 'RELIQ_MODIS':
            yticks_labels = ('0', '8', '10', '13', '15', '20', '30')
            yticks = y
            ylabel = 'Liquid particle size (micron)'
        if vd['yaxis_type'] == 'levStat':
            if not dims['levStat'].any():
                # For diagnostics on model levels, all elements in levStat
                # are set to zero. Keep vertical coordinate as level index.
                ylabel = 'Model level'
            else:
                y = np.concatenate(([0.0], dims['levStat'][::-1] + dims['levStat'][-1]))
                ylabel = 'Altitude (m)'
            yticks = y[0::4]
            yticks_labels = None
            yflip = True
        if vd['yaxis_type'] == 'lev':
            yticks = y[0::4]
            yticks_labels = None
            ylabel = 'Model level'
            yflip = True
    # Extra processing for specific variables
    vmax = None
    if vname == 'cfadLidarsr355': vmax = 0.03
    if vname == 'cfadLidarsr532': vmax = 0.03
    if vname == 'cfadLidarsr532gr': vmax = 0.03
    if vname == 'cfadDbze94': vmax = 0.05
    if vname == 'iwpmodis': vmax = 2.0
    if vname == 'lwpmodis': vmax = 1.0
    if vname == 'tauisccp': vmax = 100.0
    if vname == 'tautmodis': vmax = 100.0
    if vname == 'tauwmodis': vmax = 100.0
    d = {'xticks':xticks,
         'yticks':yticks,
         'xticks_labels':xticks_labels,
         'yticks_labels':yticks_labels,
         'xlabel':xlabel,
         'ylabel':ylabel,
         'vmax':vmax}
    # Flip y axis?
    if yflip: z = np.flip(z, axis=0)
    return x, y, z, d

def midpoints_to_edges(x):
    """
    Calculate edge points. Midpoints must increase monotonically.

    Arguments:
        x: vector with mid points. Dimension N.
    Return:
        y: numpy vector with edges. Dimension N+1.
    """
    y = np.append(np.append(2 * x[0] - x[1], x), 2 * x[-1] - x[-2])
    return 0.5 * (y[1:] + y[:-1])

def plot_pcolormesh(x, y, v, d, fig_name, title=None, coastlines=False):
    """
    Plot pcolormesh and write the output to a png file.

    Arguments:
        x: x axis values.
        y: y axis values.
        v: data array. Dimensions [Nx,Ny,Np]
        d: dictionary with plot configuration.
        fig_name: output file name.
    Keywords:
        title: plot title.
        coastlines: plot coast lines.
    """
    fig = plt.figure(figsize=(10,5))
    cmap = plt.get_cmap('YlOrRd', 20)
    cmap.set_bad('grey', 1)
    if coastlines:
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax.coastlines()
    h = plt.pcolormesh(x, y, v, cmap=cmap, vmax=d['vmax'])
    if d['xticks_labels']:
        plt.xticks(d['xticks'],d['xticks_labels'])
    else:
        plt.xticks(d['xticks'])
    if d['yticks_labels']:
        plt.yticks(d['yticks'],d['yticks_labels'])
    else:
        plt.yticks(d['yticks'])
    plt.xlabel(d['xlabel'])
    plt.ylabel(d['ylabel'])
    plt.colorbar(h,orientation='vertical')
    if title is not None: plt.title(title)
    plt.savefig(fig_name, dpi=200)
    plt.close()

def read_dimensions(fname):
    """
    Read useful dimensions from COSP output file.

    Arguments:
        fname: path to NetCDF file.
    Return:
        d: dictionary with the following dimensions:
           'cloudsat_DBZE_BINS', 'hgt16', 'REICE_MODIS',
           'RELIQ_MODIS', 'levStat', 'SR_BINS', 'lev'
    """
    dim_names = ['cloudsat_DBZE_BINS', 'hgt16', 'REICE_MODIS', 'RELIQ_MODIS',
                 'levStat', 'SR_BINS', 'lev']
    d = {}
    f_id = netCDF4.Dataset(fname, 'r')
    for dim in dim_names:
        d[dim] = f_id.variables[dim][:]
    f_id.close()
    return d

def read_var_to_masked_array(fname, vname, fill_value, Nlat_lon = None):
    """
    Reads a variable from a NetCDF file, and produces a masked array.

    Arguments:
        fname: path to NetCDF file.
        vname: variable name.
        fill_value: missing data value.
    Keywords:
        Nlat_lon: tuple (Nrows, Ncols). If defined, variable is
                  reshaped to a lat-lon grid.
    Return:
        x: variable data array.
        lon: longitude array.
        lat: latitude array
        units: units attribute.
        long_name: long name attribute.
    """
    f_id = netCDF4.Dataset(fname, 'r')
    x = np.ma.masked_equal(f_id.variables[vname][:], fill_value)
    lon = np.ma.masked_equal(f_id.variables['longitude'][:], fill_value)
    lat = np.ma.masked_equal(f_id.variables['latitude'][:], fill_value)
    units = f_id.variables[vname].getncattr('units')
    long_name = f_id.variables[vname].getncattr('long_name')
    f_id.close()

    if Nlat_lon is not None:
        x = np.reshape(x, x.shape[:-1]+Nlat_lon)
        lon = np.reshape(lon, lon.shape[:-1] + Nlat_lon)
        lat = np.reshape(lat, lat.shape[:-1] + Nlat_lon)
    return x, lon, lat, units, long_name

def produce_cosp_summary_plots(fname, variables, output_dir, Nlat_lon = None):
    """
    Wrapper function that iterates over a list of COSP variables and produces
    a PNG figure for each of them.

    Arguments:
        fname: COSP output filename.
        variables: list of variable names.
        output_dir: output directory.
    Keywords:
        Nlat_lon: tuple with (lat,lon) dimensions model's grid.
    """
    fill_value = -1.0e30
    dimensions = read_dimensions(fname)
    for vname, vd in variables.items():
        new_shape = None
        if vd['reshape']: new_shape = Nlat_lon
        vx, longitude, latitude, units, long_name = read_var_to_masked_array(fname, vname, fill_value, Nlat_lon = new_shape)
        x, y, z, pkw = collapse_dimensions_for_plotting(longitude, latitude, vname, vx, vd, dimensions)
        title = long_name + ' (' + units + ')'
        fig_name = os.path.join(output_dir, ".".join([os.path.basename(fname), vname, 'png']))
        coastlines = False
        if vd['plot_type'] == 'map': coastlines = True
        plot_pcolormesh(x, y, z, pkw, fig_name, title=title, coastlines=coastlines)

def variable2D_metadata(var_list, fname):
    """
    Return dictionary with metadata for each variable.

    Arguments:
        var_list: list of variable names.
        fname: COSP output filename.
    Return:
        d: dictionary of dictionaries with relevant metadata.
    """
    map_dims = (('loc',),('PARASOL_NREFL','loc'))
    hist2D_dims = (('pres7', 'tau7', 'loc'),
                  ('levStat', 'cloudsat_DBZE_BINS', 'loc'),
                  ('hgt16', 'tau7', 'loc'),
                  ('REICE_MODIS', 'tau7', 'loc'),
                  ('RELIQ_MODIS', 'tau7', 'loc'),
                  ('levStat', 'SR_BINS', 'loc'))
    zcs_dims = (('levStat','loc'), ('lev','loc'))
    f_id = netCDF4.Dataset(fname, 'r')
    vmeta = {}
    print("=== Processing variables in output file:\n {}".format(fname))
    for vname in var_list:
        try:
            x = f_id.variables[vname]
            # Standard map
            if x.dimensions in map_dims:
                vmeta[vname] = {'plot_type':'map', 'reshape':True}
            # 2D histograms
            if x.dimensions in hist2D_dims:
                vmeta[vname] = {'plot_type':'2Dhist', 'reshape':False,
                                'xaxis_type': x.dimensions[1],
                                'yaxis_type': x.dimensions[0]}
            # Zonal cross section
            if x.dimensions in zcs_dims:
                vmeta[vname] = {'plot_type':'zonal_cross_section', 'reshape':True,
                                'xaxis_type': 'latitude',
                                'yaxis_type': x.dimensions[0]}
        except:
            print("Skipping {}, not found in output file.".format(vname))
    f_id.close()
    return vmeta



#######################
# Main
#######################
if __name__ == '__main__':

    # Command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("--tst_file", default="./data/outputs/UKMO/cosp2_output.um_global.nc",
                        help="Test output file.")
    parser.add_argument("--out_dir", default="./data/outputs/UKMO",
                        help="Output directory.")
    parser.add_argument("--Nlat", type=int, default=36, help="Number of latitude points.")
    parser.add_argument("--Nlon", type=int, default=48, help="Number of longitude points.")
    args = parser.parse_args()

    # Dictonaries with list of variables.
    v2D_maps_names = ['cllcalipsoice', 'clmcalipsoice', 'clhcalipsoice', 'cltcalipsoice', 'cllcalipsoliq', 'clmcalipsoliq',
                      'clhcalipsoliq', 'cltcalipsoliq', 'cllcalipsoun', 'clmcalipsoun', 'clhcalipsoun', 'cltcalipsoun',
                      'cllcalipso', 'clmcalipso', 'clhcalipso', 'cltcalipso', 'clopaquecalipso', 'clthincalipso',
                      'clzopaquecalipso', 'clopaquetemp', 'clthintemp', 'clzopaquetemp', 'clopaquemeanz', 'clthinmeanz',
                      'clthinemis', 'clopaquemeanzse', 'clthinmeanzse', 'clzopaquecalipsose', 'cllgrLidar532',
                      'clmgrLidar532', 'clhgrLidar532', 'cltgrLidar532', 'cllatlid', 'clmatlid', 'clhatlid',
                      'cltatlid', 'cloudsatpia', 'cltisccp', 'meantbisccp', 'meantbclrisccp', 'pctisccp', 'tauisccp',
                      'albisccp', 'misr_meanztop', 'misr_cldarea', 'cltmodis', 'clwmodis', 'climodis', 'clhmodis',
                      'clmmodis', 'cllmodis', 'tautmodis', 'tauwmodis', 'tauimodis', 'tautlogmodis', 'tauwlogmodis',
                      'tauilogmodis', 'reffclwmodis', 'reffclimodis', 'pctmodis', 'lwpmodis', 'iwpmodis',
                      'cltlidarradar', 'cloudsat_tcc', 'cloudsat_tcc2','parasolGrid_refl']
    v2D_hists_names = ['clisccp', 'clmodis', 'cfadDbze94', 'clMISR',
                       'modis_Optical_Thickness_vs_ReffICE',
                       'modis_Optical_Thickness_vs_ReffLIQ',
                       'cfadLidarsr532', 'cfadLidarsr532gr', 'cfadLidarsr355']
    v2D_zcs_names = ['clcalipsoice','clcalipsoliq','clcalipsoun','clcalipsotmp','clcalipsotmpice','clcalipsotmpliq',
                     'clcalipsotmpun','clcalipso','clcalipsoopaque','clcalipsothin','clcalipsozopaque',
                     'clcalipsoopacity','clgrLidar532','clatlid','clcalipso2',
                     'lidarBetaMol532gr','lidarBetaMol532','lidarBetaMol355']
    v2D_all_names = v2D_maps_names + v2D_hists_names + v2D_zcs_names
    # Plots for these variables are not yet developed
    # atb532_perp(lev, cosp_scol, loc);
    # atb532(lev, cosp_scol, loc);
    # calipso_tau(lev, cosp_scol, loc);
    # atb532gr(lev, cosp_scol, loc);
    # atb355(lev, cosp_scol, loc);
    # dbze94(lev, cosp_scol, loc);
    # parasolPix_refl(PARASOL_NREFL, cosp_scol, loc);
    # boxtauisccp(cosp_scol, loc);
    # boxptopisccp(cosp_scol, loc);

    # Build dictionary with metadata and produce plots.
    vars2D  = variable2D_metadata(v2D_all_names, args.tst_file)
    produce_cosp_summary_plots(args.tst_file, vars2D, args.out_dir,
                               Nlat_lon=(args.Nlat, args.Nlon))

    print("===== Summary plots produced =====")
