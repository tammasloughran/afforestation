#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob
import os
import pdb
import warnings

import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import numpy as np
import netCDF4 as nc
from scipy import stats
import ipdb

from shapely.errors import ShapelyDeprecationWarning

if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filename
    from analysis.constants import (
            DATA_DIR,
            DPI,
            ENSEMBLES,
            NENS,
            PLOTS_DIR,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filename
    from constants import (
            DATA_DIR,
            DPI,
            ENSEMBLES,
            NENS,
            PLOTS_DIR,
            )

warnings.filterwarnings(action='ignore', category=ShapelyDeprecationWarning)
warnings.filterwarnings(action='ignore', category=stats.SpearmanRConstantInputWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)

files = glob.glob('./*')
if PLOTS_DIR not in files:
    os.mkdir(PLOTS_DIR)
if DATA_DIR not in files:
    os.mkdir(DATA_DIR)

# Control flag
files = glob.glob(f'{DATA_DIR}/cVeg_*_global.npy')
if any(['.npy' in f for f in files]):
    load_npy_files = True
else:
    load_npy_files = False
load_npy_files = True # Uncomment to override previous check.

NTIMES = 1032

# Load data for only afforested grid cells.
treeFrac = np.load(f'{DATA_DIR}/treeFrac_area_anomaly.npy')
NLAT = treeFrac.shape[0]
NLON = treeFrac.shape[1]


def yearly_mean_from_monthly(data:np.ndarray)->np.ndarray:
    """Calculate a yearly mean on a numpy array of monthly data.
    The 0th dimension must be time and divisible by 12.
    """
    if np.mod(data.shape[0], 12)!=0:
        raise ValueError("Not enough months in 0th dimension.")
    toshape = list(data.shape)
    toshape.pop(0)
    toshape.insert(0, 12)
    toshape.insert(0, int(data.shape[0]/12))
    fraction_of_year = np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.0
    return np.average(data.reshape(toshape), axis=1, weights=fraction_of_year)


def my_spearmanr(a:np.ndarray, b:np.ndarray)->tuple:
    """Calculate spearman correlation allong the time axis for a and b.
    a is an array of shape (e,t,x,y)
    b is an array of shape (t,x,y)

    Return
    ------
    r is an array of shape (e,x,y)
    p is an array of shape (e,x,y)
    """
    r = np.ones((NENS,NLAT,NLON))*np.nan
    p = np.ones((NENS,NLAT,NLON))*np.nan
    for e in range(NENS):
        print('    - Correlating ensemble: ', e+1)
        for i in range(NLAT):
            for j in range(NLON):
                r[e,i,j], p[e,i,j] = stats.spearmanr(a[e,:,i,j], b[:,i,j])
    return r, p


def my_spearmanr2(a:np.ndarray, b:np.ndarray)->tuple:
    """Calculate spearman correlation allong the time axis for a and b.
    a is an array of shape (t,x,y)
    b is an array of shape (t,x,y)

    Return
    ------
    r is an array of shape (x,y)
    p is an array of shape (x,y)
    """
    r = np.ones((NLAT,NLON))*np.nan
    p = np.ones((NLAT,NLON))*np.nan
    for i in range(NLAT):
        for j in range(NLON):
            r[i,j], p[i,j] = stats.spearmanr(a[:,i,j], b[:,i,j])
    return r, p


def make_correlation_plot()->None:
    if not load_npy_files:
        # Load tree frac
        print('Loading treeFrac')
        frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', 'Lmon', 'treeFrac')[0]
        ncfile = nc.Dataset(frac_file, 'r')
        for_tree = ncfile.variables['treeFrac'][:]
        frac_file = get_filename('C4MIP', 'esm-ssp585', '1', 'Lmon', 'treeFrac')[0]
        ncfile = nc.Dataset(frac_file, 'r')
        ssp585_tree = ncfile.variables['treeFrac'][:]

        # Calculate difference in tree frac
        tree_diff = for_tree - ssp585_tree

        # Load temperature data
        for_temp = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
        ssp585_temp = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
        print('Loading temperature')
        for e,ens in enumerate(ENSEMBLES):
            print('    - ensemble: ', ens)
            temp_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Amon', 'tas')[0]
            for_temp[e,...] = nc.Dataset(temp_file, 'r').variables['tas'][:]
            temp_file = get_filename('C4MIP', 'esm-ssp585', ens, 'Amon', 'tas')[0]
            ssp585_temp[e,...] = nc.Dataset(temp_file, 'r').variables['tas'][:]

        # Calculate difference in surface temperature
        temp_diff = for_temp - ssp585_temp

        # Save np files
        np.save(f'{DATA_DIR}/temp_diff_monthly.npy', temp_diff.data)
        np.save(f'{DATA_DIR}/tree_diff_monthly.npy', tree_diff.data)
    else:
        frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', 'Lmon', 'treeFrac')[0]
        ncfile = nc.Dataset(frac_file, 'r')
        lats = ncfile.variables['lat'][:]
        lons = ncfile.variables['lon'][:]
        temp_diff = np.load(f'{DATA_DIR}/temp_diff_monthly.npy')
        tree_diff = np.load(f'{DATA_DIR}/tree_diff_monthly.npy')

    # Correlate.
    print('Correlating')
    correlation, pvalue = my_spearmanr(albedo_diff, tree_diff)
    correlation[pvalue>0.05] = 0

    # Plot
    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    absmax = max(abs(np.nanmin(correlation)), np.nanmax(correlation))
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-55, 65, 10)/100.0, ncolors=256)
    plt.pcolormesh(lons, lats, np.nanmean(correlation, axis=0),
            norm=discrete_bins,
            #vmin=-absmax,
            #vmax=absmax,
            cmap='bwr',
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(
            ticks=np.arange(-5, 6, 1)/10,
            orientation='horizontal',
            pad=0.05,
            )
    plt.title('Correlation between albedo and treeFrac')
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/correlation_tree_albedo.png', dpi=DPI)

    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    absmax = max(abs(np.nanmin(correlation)), np.nanmax(correlation))
    plt.pcolormesh(lons, lats, np.nansum(pvalue<0.05, axis=0),
            cmap='viridis',
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(orientation='horizontal', pad=0.05)
    plt.title('Number of significant correlations')
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/correlation_tree_tas_significance.png', dpi=DPI)


def make_ens_mean_first_correlation():
    """Same corrlation as make_correlation_plot() except it does the ensemble mean first.
    """
    if not load_npy_files:
        # Load tree frac
        print('Loading treeFrac')
        frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', 'Lmon', 'treeFrac')[0]
        ncfile = nc.Dataset(frac_file, 'r')
        for_tree = ncfile.variables['treeFrac'][:]
        frac_file = get_filename('C4MIP', 'esm-ssp585', '1', 'Lmon', 'treeFrac')[0]
        ncfile = nc.Dataset(frac_file, 'r')
        ssp585_tree = ncfile.variables['treeFrac'][:]

        # Calculate difference in tree frac
        tree_diff = for_tree - ssp585_tree

        print('Loading latent heat flux')
        # Load albedo correlation
        ncfile = nc.Dataset('/g/data/p66/rml599/ESM-ensemble/esm-ssp585/evspsbl_Amon_ACCESS-ESM1-5_esm-ssp585_ens_gn_201501-210012.ann.nc', 'r')
        albedo_ssp585 = ncfile.variables['evspsbl'][:]
        ncfile = nc.Dataset('/g/data/p66/rml599/ESM-ensemble/esm-ssp585-ssp126Lu/evspsbl_Amon_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ens_gn_201501-210012.ann.nc', 'r')
        albedo_for = ncfile.variables['evspsbl'][:]
        ncfile.close()

        # Calculate difference in surface temperature
        albedo_diff = albedo_for - albedo_ssp585

        # Save np files
        np.save(f'{DATA_DIR}/temp_diff_monthly.npy', temp_diff.data)
        np.save(f'{DATA_DIR}/tree_diff_monthly.npy', tree_diff.data)
        lats = ncfile.variables['lat'][:]
        lons = ncfile.variables['lon'][:]
        np.save(f'{DATA_DIR}/access_lats.npy', np.array(lats))
        np.save(f'{DATA_DIR}/access_lons.npy', np.array(lons))
    else:
        frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', 'Lmon', 'treeFrac')[0]
        ncfile = nc.Dataset(frac_file, 'r')
        lats = ncfile.variables['lat'][:]
        lons = ncfile.variables['lon'][:]

        # Load albedo correlation
        ncfile = nc.Dataset('/g/data/p66/rml599/ESM-ensemble/esm-ssp585/evspsbl_Amon_ACCESS-ESM1-5_esm-ssp585_ens_gn_201501-210012.ann.nc', 'r')
        albedo_ssp585 = ncfile.variables['evspsbl'][:]
        ncfile = nc.Dataset('/g/data/p66/rml599/ESM-ensemble/esm-ssp585-ssp126Lu/evspsbl_Amon_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ens_gn_201501-210012.ann.nc', 'r')
        albedo_for = ncfile.variables['evspsbl'][:]
        ncfile.close()

        albedo_diff = albedo_for - albedo_ssp585

        tree_diff = np.load(f'{DATA_DIR}/tree_diff_monthly.npy')

        lats = np.load(f'{DATA_DIR}/access_lats.npy')
        lons = np.load(f'{DATA_DIR}/access_lons.npy')

    # Correlate.
    print('Correlating')
    correlation, pvalue = my_spearmanr2(albedo_diff, yearly_mean_from_monthly(tree_diff))
    correlation[pvalue>0.05] = 0 # Show only values significant at the 5% level

    # Plot
    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    absmax = max(abs(np.nanmin(correlation)), np.nanmax(correlation))
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-95, 105, 10)/100.0, ncolors=256)
    plt.pcolormesh(lons, lats, correlation,
            norm=discrete_bins,
            #vmin=-absmax,
            #vmax=absmax,
            cmap='bwr',
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(
            ticks=np.arange(-9, 10, 1)/10,
            orientation='horizontal',
            pad=0.05,
            )
    plt.title('Correlation between evapotranspiration and tree fraction')
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/correlation_tree_evapotranspiration_ens_mean_first.png', dpi=DPI)

    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    absmax = max(abs(np.nanmin(correlation)), np.nanmax(correlation))
    plt.pcolormesh(lons, lats, pvalue<0.05,
            cmap='viridis',
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(orientation='horizontal', pad=0.05)
    plt.title('Is significant')
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/correlation_tree_evapostranspiration_significance_ens_mean_first.png', dpi=DPI)
    plt.show()


if __name__=='__main__':
    #make_correlation_plot()
    make_ens_mean_first_correlation()

