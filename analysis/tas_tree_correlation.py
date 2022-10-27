#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob
import os
import pdb
import sys
import warnings

from cdo import Cdo
import cdo_decorators as cdod
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import netCDF4 as nc
import scipy import stats

from shapely.errors import ShapelyDeprecationWarning

if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cdo_calc_load import cdo_cover_area_load, cdo_area_diff_load
    from analysis.constants import FRAC_VARIABLES, M2_IN_MILKM2, PLOTS_DIR, DATA_DIR, DPI
    from analysis.cmip_files import get_filename, LAND_FRAC_FILE
    from analysis.jaisnb import jaisnb
    from analysis.constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            M2_IN_MILKM2,
            NENS,
            NTIMES,
            PLOTS_DIR,
            SEC_IN_DAY,
            SEC_IN_YEAR,
            TABLES,
            VARIABLES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cdo_calc_load import cdo_cover_area_load, cdo_area_diff_load
    from constants import FRAC_VARIABLES, M2_IN_MILKM2, PLOTS_DIR, DATA_DIR, DPI
    from cmip_files import get_filename, LAND_FRAC_FILE
    from jaisnb import jaisnb
    from constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            M2_IN_MILKM2,
            NENS,
            NTIMES,
            PLOTS_DIR,
            SEC_IN_DAY,
            SEC_IN_YEAR,
            TABLES,
            VARIABLES,
            )

warnings.filterwarnings(action='ignore', category=ShapelyDeprecationWarning)
cdo = Cdo()
cdo.debug = False

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
load_npy_files = False # Uncomment to override previous check.

MIPS = {
        'esm-ssp585':'C4MIP',
        'esm-ssp585-ssp126Lu':'LUMIP',
        }

LABELS = {
        'AFF':
                {'treeFrac': 'tree AFF',
                'cropFrac': 'crop AFF',
                'shrubFrac': 'shrub AFF',
                'grassFrac': 'grass AFF'},
        'esm-ssp585':
                {'treeFrac': 'tree SSP5-8.5',
                'cropFrac': 'crop SSP5-8.5',
                'shrubFrac': 'shrub SSP5-8.5',
                'grassFrac': 'grass SSP5-8.5'}}

# Load data for only afforested grid cells.
treeFrac = np.load(f'{DATA_DIR}/treeFrac_area_anomaly.npy')/M2_IN_MILKM2
NLAT = treeFrac.shape[0]
NLON = treeFrac.shape[1]


def my_spearmanr(a:np.ndarray, b:np.ndarray)->np.ndarray:
    """Calculate spearman correlation allong the time axis for a and b.
    a is an array of shape (e,t,x,y)
    b is an array of shape (t,x,y)

    Return
    ------
    r is an array of shape (e,x,y)
    p is an array of shape (e,x,y)
    """
    r = np.ones((NENS,NLAT,NLON))*.np.nan
    p = np.ones((NENS,NLAT,NLON))*.np.nan
    for e in range(NENS):
        for i in range(NLAT):
            for j in range(NLON):
                r[e,i,j], p[e,i,j] = stats.spearmanr(a[e,:,i,j], b[:,i,j])
    return r, p


def make_correlation_plot()->None:
    # Load tree frac
    frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', 'Lmon', 'treeFrac')
    ncfile = nc.Dataset(frac_file, 'r')
    for_tree = ncfile.variables['treeFrac'][:]
    frac_file = get_filename('C4MIP', 'esm-ssp585', '1', 'Lmon', 'treeFrac')
    ncfile = nc.Dataset(frac_file, 'r')
    ssp585_tree = ncfile.variables['treeFrac'][:]

    # Calculate difference in tree frac
    tree_diff = for_tree - ssp585_tree

    lats = ncfile.variables['lat'][:]
    lons = ncfile.variables['lon'][:]

    # Load temperature data
    for_temp = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
    ssp585_temp = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
    for e,ens in enumerate(ENSEMBLES):
        temp_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Amon', 'tas')
        for_temp[e,...] = nc.Dataset(temp_file, 'r').variables['tas'][:]
        temp_file = get_filename('C4MIP', 'esm-ssp585', ens, 'Amon', 'tas')
        ssp585_temp[e,...] = nc.Dataset(temp_file, 'r').variables['tas'][:]

    # Calculate difference in surface temperature
    temp_diff = for_temp - ssp585_temp

    # Correlate, then ensemble mean.
    correlation, pvalue = my_spearmanr(temp_diff, tree_diff).mean(axis=0)

    # Plot
    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    plt.pcolormesh(
            lat,
            lon,
            correlation,
            vmin=-1,
            vmax=1,
            cmap='bwr',
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(orientation='horizontal', pad=0.05)
    plt.title('Correlation between TAS and treeFrac')
    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    make_correlation_plot()

