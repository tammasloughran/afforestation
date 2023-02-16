#!/usr/bin/env python3
# Plot the tree fraction and areas for all models participating in the LUMIP
# forestation experiment esm-ssp585-ssp126Lu.
import datetime as dt
import glob
import pdb
import warnings
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cftime
import ipdb
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4 as nc
import numpy as np
from cdo import Cdo
import shapely
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

if __name__ == 'analysis.plot_model_trees':
    from analysis.cmip_files import get_filenames
    from analysis.constants import DATA_DIR, PLOTS_DIR, DPI
else:
    from cmip_files import get_filenames
    from constants import DATA_DIR, PLOTS_DIR, DPI

warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')
warnings.filterwarnings('ignore', category=shapely.errors.ShapelyDeprecationWarning)
warnings.filterwarnings('ignore', category=mpl.MatplotlibDeprecationWarning)

cdo = Cdo()

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Constants
# Some models are missing treeFrac variables on the ESGF.
MODELS = [
        'ACCESS-ESM1-5',
        #'BCC-CSM2-MR',
        'CanESM5',
        'CESM2',
        'GFDL-ESM4',
        #'MIROC-ESL2',
        'MPI-ESM1-2-LR',
        #'NorESM2-LM',
        'UKESM1-0-LL',
        ]
INSTIT = {
        'ACCESS-ESM1-5':'CSIRO',
        'BCC-CSM2-MR':'BCC',
        'CanESM5':'CCma',
        'GFDL-ESM4':'NOAA-GFDL',
        'MIROC-ES2L':'MIROC',
        'MPI-ESM1-2-LR':'MPI-M',
        'NorESM2-LM':'NCC',
        'CESM2':'NCAR',
        'UKESM1-0-LL':'MOHC',
        }
COLORS = {
        'CSIRO':color_cycle[0],
        'BCC':color_cycle[1],
        'CCma':color_cycle[2],
        'NOAA-GFDL':color_cycle[3],
        'MIROC':color_cycle[4],
        'MPI-M':color_cycle[5],
        #'NCC':color_cycle[6],
        'NCAR':color_cycle[6],
        'MOHC':color_cycle[8],
        }
ENSEMBLES = {
        'ACCESS-ESM1-5':'r1i1p1f1',
        'BCC-CSM2-MR':'r1i1p1f1',
        'CanESM5':'r1i1p2f1',
        'MIROC-ES2L':'r1i1p1f2',
        'UKESM1-0-LL':'r1i1p1f2',
        'MPI-ESM1-2-LR':'r1i1p1f1',
        'NorESM2-LM':'r1i1p1f1',
        'CESM2':'r1i1p1f1',
        'GFDL-ESM4':'r1i1p1f1',
        }
M2_TOMILKM2 = 1000000000000
#cdict = {'red':   [[0.0,  84/255, 84/255],
#                   [0.5,  1.0, 1.0],
#                   [1.0,  0.0, 0.0]],
#         'green': [[0.0, 48/255, 48/255],
#                   [0.5, 1.0, 1.0],
#                   [1.0,  68/255, 68/255]],
#         'blue':  [[0.0,  5/255, 5/255],
#                   [0.5,  1.0, 1.0],
#                   [1.0,  27/255, 27/255]]}
#CMAP = LinearSegmentedColormap('diverging_trees', segmentdata=cdict, N=256)
CMAP = 'bwr'


def my_num2date(times:np.ndarray, units:str, calendar:str)->list:
    dates = nc.num2date(times, units, calendar=calendar)
    return [dt.datetime(d.year, d.month, d.day, d.hour) for d in dates]


def plot_global_sum_area()->None:
    plt.figure()
    for model in MODELS:
        cell_area = nc.Dataset(f'{DATA_DIR}/gridarea_{model}.nc', 'r').variables['cell_area'][:]
        if model=='CESM2':
            table = 'Eyr'
        else:
            table = 'Lmon'
        files = sorted(get_filenames(
                'LUMIP',
                INSTIT[model],
                model,
                'esm-ssp585-ssp126Lu',
                ENSEMBLES[model],
                table,
                'treeFrac',
                ))
        if len(files)>1:
            nc_tree_frac = nc.MFDataset(files, 'r')
        else:
            nc_tree_frac = nc.Dataset(files[0], 'r')
        tree_frac = nc_tree_frac.variables['treeFrac'][:].squeeze()/100
        tree_area = tree_frac*cell_area/M2_TOMILKM2
        global_tree_area = tree_area.sum(axis=(1,2))
        times = nc_tree_frac.variables['time']
        dates = my_num2date(times[:], times.units, times.calendar)
        plt.plot(dates, global_tree_area, color=COLORS[INSTIT[model]], label=model)
    plt.xlim(left=dates[0], right=dates[-1])
    plt.ylabel('Area (million km$^2$)')
    plt.title('Tree area')
    plt.legend()
    plt.savefig(f'{PLOTS_DIR}/tree_area_ssp126_all_models_global_sum.png', dpi=DPI)


def plot_global_mean_frac()->None:
    plt.figure()
    for model in MODELS:
        if model=='CESM2':
            table = 'Eyr'
        else:
            table = 'Lmon'
        files = sorted(get_filenames(
                'LUMIP',
                INSTIT[model],
                model,
                'esm-ssp585-ssp126Lu',
                ENSEMBLES[model],
                table,
                'treeFrac',
                ))
        if len(files)>1:
            nc_tree_frac = nc.MFDataset(files, 'r')
        else:
            nc_tree_frac = nc.Dataset(files[0], 'r')
        tree_frac = nc_tree_frac.variables['treeFrac'][:].squeeze()
        lats = nc_tree_frac.variables['lat'][:]
        coslats = np.cos(np.deg2rad(lats))[None,:,None]*np.ones(tree_frac.shape)
        global_tree_frac = np.ma.average(tree_frac, axis=(1,2), weights=coslats)
        times = nc_tree_frac.variables['time']
        dates = my_num2date(times[:], times.units, times.calendar)
        plt.plot(dates, global_tree_frac, color=COLORS[INSTIT[model]], label=model)
    plt.xlim(left=dates[0], right=dates[-1])
    plt.ylabel('Fraction %')
    plt.title('Tree fraction')
    plt.legend()
    plt.savefig(f'{PLOTS_DIR}/tree_fraction_ssp126_all_models_globalmean.png', dpi=DPI)


def plot_area_by_lat()->None:
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        cell_area = nc.Dataset(f'{DATA_DIR}/gridarea_{model}.nc', 'r').variables['cell_area'][:]
        if model=='CESM2':
            table = 'Eyr'
        else:
            table = 'Lmon'
        files = sorted(get_filenames(
                'LUMIP',
                INSTIT[model],
                model,
                'esm-ssp585-ssp126Lu',
                ENSEMBLES[model],
                table,
                'treeFrac',
                ))
        if len(files)>1:
            nc_tree_frac = nc.MFDataset(files, 'r')
        else:
            nc_tree_frac = nc.Dataset(files[0], 'r')
        tree_frac = nc_tree_frac.variables['treeFrac'][:].squeeze()/100
        tree_area = tree_frac*cell_area/M2_TOMILKM2
        lat_tree_area = tree_area.sum(axis=-1)
        times = nc_tree_frac.variables['time']
        dates = my_num2date(times[:], times.units, times.calendar)
        lats = nc_tree_frac.variables['lat'][:]
        plt.pcolormesh(
                dates,
                lats,
                (lat_tree_area - lat_tree_area[0]).T,
                label=model,
                vmin=-0.5,
                vmax=0.5,
                cmap=CMAP,
                shading='auto',
                )
        plt.title(model, fontsize=8)
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.2, left=0.1, right=0.95, hspace=0.2, wspace=0.2)
    fig.add_subplot(111, frameon=False)
    plt.tick_params(
            labelcolor='none',
            which='both',
            top=False,
            bottom=False,
            left=False,
            right=False,
            )
    plt.xlabel('Year')
    plt.ylabel('Latitude °N')
    cbar_ax = fig.add_axes([0.1, 0.085, 0.85, 0.04])
    plt.colorbar(cax=cbar_ax, label='Area (million km$^2$)', orientation='horizontal', pad=0.05)
    plt.savefig(f'{PLOTS_DIR}/tree_area_ssp126_all_models_zonsum.png', dpi=DPI)


def plot_frac_by_lat()->None:
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        if model=='CESM2':
            table = 'Eyr'
        else:
            table = 'Lmon'
        files = sorted(get_filenames(
                'LUMIP',
                INSTIT[model],
                model,
                'esm-ssp585-ssp126Lu',
                ENSEMBLES[model],
                table,
                'treeFrac',
                ))
        if len(files)>1:
            nc_tree_frac = nc.MFDataset(files, 'r')
        else:
            nc_tree_frac = nc.Dataset(files[0], 'r')
        tree_frac = nc_tree_frac.variables['treeFrac'][:].squeeze()
        lat_tree_frac = np.ma.average(tree_frac, axis=-1)
        times = nc_tree_frac.variables['time']
        dates = my_num2date(times[:], times.units, times.calendar)
        lats = nc_tree_frac.variables['lat'][:]
        plt.pcolormesh(
                dates,
                lats,
                (lat_tree_frac - lat_tree_frac[0]).T,
                label=model,
                vmin=-40,
                vmax=40,
                cmap=CMAP,
                shading='auto',
                )
        plt.title(model, fontsize=8)
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.2, left=0.1, right=0.95, hspace=0.2, wspace=0.1)
    fig.add_subplot(111, frameon=False)
    plt.tick_params(
            labelcolor='none',
            which='both',
            top=False,
            bottom=False,
            left=False,
            right=False,
            )
    plt.xlabel('Year')
    plt.ylabel('Latitude °N')
    cbar_ax = fig.add_axes([0.1, 0.085, 0.85, 0.04])
    plt.colorbar(cax=cbar_ax, label='Fraction %', orientation='horizontal', pad=0.05)
    plt.savefig(f'{PLOTS_DIR}/tree_area_ssp126_all_models_zonmean.png', dpi=DPI)


def plot_map_tree_area_change()->None:
    fig, axes = plt.subplots(
            nrows=3,
            ncols=2,
            sharex=True,
            sharey=True,
            subplot_kw={'projection': ccrs.EckertIV()},
            )
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        if model=='CESM2':
            table = 'Eyr'
        else:
            table = 'Lmon'
        files = sorted(get_filenames(
                'LUMIP',
                INSTIT[model],
                model,
                'esm-ssp585-ssp126Lu',
                ENSEMBLES[model],
                table,
                'treeFrac',
                ))
        if len(files)>1:
            nc_tree_frac = nc.MFDataset(files, 'r')
        else:
            nc_tree_frac = nc.Dataset(files[0], 'r')
        tree_frac = nc_tree_frac.variables['treeFrac'][:].squeeze()/100
        cell_area = nc.Dataset(f'{DATA_DIR}/gridarea_{model}.nc', 'r').variables['cell_area'][:]
        tree_area = tree_frac*cell_area/M2_TOMILKM2
        lats = nc_tree_frac.variables['lat'][:]
        lons = nc_tree_frac.variables['lon'][:]
        lons = np.append(lons, [lons[-1] + lons[1] - lons[0]]) # pcolormesh expects +1 lon.
        lons = lons - (lons[1] - lons[0])/2 # Shift plotting to centre grid, align with coastlines.
        lats = lats - (lats[1] - lats[0])/4
        data = plt.pcolormesh(
                lons,
                lats,
                tree_area[-1]-tree_area[0],
                cmap=CMAP,
                vmin=-0.025,
                vmax=0.025,
                transform=ccrs.PlateCarree(),
                )
        axes[i].coastlines()
        plt.title(model, fontsize=8)
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, hspace=0.12, wspace=0.05)
    cbar_ax = fig.add_axes([0.1, 0.085, 0.85, 0.04])
    plt.colorbar(
            data,
            cax=cbar_ax,
            label='Area (million km$^2$)',
            orientation='horizontal',
            pad=0.05,
            )
    plt.savefig(f'{PLOTS_DIR}/tree_area_map_2015-2100_delta.png', dpi=DPI)


def make_tree_frac_plots()->None:
    # Calculate model grid cell areas.
    for model in MODELS:
        if model=='CESM2':
            table = 'Eyr'
        else:
            table = 'Lmon'
        files = sorted(get_filenames(
            'LUMIP',
            INSTIT[model],
            model,
            'esm-ssp585-ssp126Lu',
            ENSEMBLES[model],
            table,
            'treeFrac',
            ))
        if f'gridarea_{model}.nc' not in os.listdir(DATA_DIR):
            cdo.gridarea(input=files[0], output=f'{DATA_DIR}/gridarea_{model}.nc')

    # Line plots of global sum of tree frac area
    plot_global_sum_area()

    # Plot global mean fraction.
    plot_global_mean_frac()

    # Plot area by latitude.
    plot_area_by_lat()

    # Plot frac by latitude.
    plot_frac_by_lat()

    # Plot global maps of tree area change.
    plot_map_tree_area_change()


if __name__ != 'analysis.plot_model_trees':
    make_tree_frac_plots()

