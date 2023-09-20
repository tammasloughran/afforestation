#!/usr/bin/env python3
"""Plot the tree fraction and areas for all models participating in the LUMIP
forestation experiment esm-ssp585-ssp126Lu.
"""
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
        #'BCC-CSM2-MR', # BCC has been excluded.
        #'MIROC-ESL2', # MIROC has missing data.
        #'NorESM2-LM', # NorESM has been excluded.
        'ACCESS-ESM1-5',
        'CESM2',
        'CanESM5',
        'GFDL-ESM4',
        'MPI-ESM1-2-LR',
        'UKESM1-0-LL',
        ]
INSTIT = {
        #'BCC-CSM2-MR':'BCC',
        #'NorESM2-LM':'NCC',
        'ACCESS-ESM1-5':'CSIRO',
        'CESM2':'NCAR',
        'CanESM5':'CCma',
        'GFDL-ESM4':'NOAA-GFDL',
        'MIROC-ES2L':'MIROC',
        'MPI-ESM1-2-LR':'MPI-M',
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
        #'BCC-CSM2-MR':'r1i1p1f1',
        #'NorESM2-LM':'r1i1p1f1',
        'ACCESS-ESM1-5':'r1i1p1f1',
        'CESM2':'r1i1p1f1',
        'CanESM5':'r1i1p2f1',
        'GFDL-ESM4':'r1i1p1f1',
        'MIROC-ES2L':'r1i1p1f2',
        'MPI-ESM1-2-LR':'r1i1p1f1',
        'UKESM1-0-LL':'r1i1p1f2',
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


import numpy as _np
def yearly_mean_from_monthly(data:_np.ndarray)->_np.ndarray:
    """Calculate a yearly mean on a numpy array of monthly data.
    The 0th dimension must be time and divisible by 12.
    """
    if _np.mod(data.shape[0], 12)!=0:
        raise ValueError("Not enough months in 0th dimension.")
    toshape = list(data.shape)
    toshape.pop(0)
    toshape.insert(0, 12)
    toshape.insert(0, int(data.shape[0]/12))
    fraction_of_year = _np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.0
    return _np.average(data.reshape(toshape), axis=1, weights=fraction_of_year)


def my_num2date(times:np.ndarray, units:str, calendar:str)->list:
    """Convert numbers to a list of datetime objects.
    """
    dates = nc.num2date(times, units, calendar=calendar)
    return [dt.datetime(d.year, d.month, d.day, d.hour) for d in dates]


def global_mean(data:np.ndarray, lats:np.ndarray)->np.ndarray:
    """Calculate an area weighted global mean. weights are the cosine of lats.
    """
    coslats = np.cos(np.deg2rad(lats))[None,:,None]*np.ones(data.shape)
    return np.ma.average(data, axis=(1,2), weights=coslats)


def load_frac(file_names:list)->tuple:
    """Load all required data.

    Returns:
    - (tree_frac,dates,lats,lons)
    """
    if len(file_names)>1:
        nc_tree_frac = nc.MFDataset(file_names, 'r')
    else:
        nc_tree_frac = nc.Dataset(file_names[0], 'r')
    tree_frac = nc_tree_frac.variables['treeFrac'][:].squeeze()
    tree_frac[tree_frac.mask] = 0 # Remove missing values
    lats = nc_tree_frac.variables['lat'][:]
    lons = nc_tree_frac.variables['lon'][:]
    times = nc_tree_frac.variables['time']
    dates = my_num2date(times[:], times.units, times.calendar)
    return tree_frac, dates, lats, lons


def calc_tree_area(tree_frac:np.ndarray, area_file:str)->np.ndarray:
    """Load gridcell area in units million km^2.
    """
    areas = nc.Dataset(area_file, 'r').variables['cell_area'][:]
    return (tree_frac/100)*areas/M2_TOMILKM2


def plot_global_sum_area()->None:
    """Plot the global sum of tree areas.
    """
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    for model in MODELS:
        if model=='UKESM1-0-LL' or model=='GFDL-ESM4':
            # Convert from monthly to yearly
            for_dates_plot = for_dates[model][5:-1:12]
            ssp585_dates_plot = ssp585_dates[model][5:-1:12]
            for_global_sum_plot = yearly_mean_from_monthly(for_global_sum[model])
            ssp585_global_sum_plot = yearly_mean_from_monthly(ssp585_global_sum[model])
            for_global_mean_plot = yearly_mean_from_monthly(for_global_mean[model])
            ssp585_global_mean_plot = yearly_mean_from_monthly(ssp585_global_mean[model])
        else:
            for_dates_plot = for_dates[model]
            ssp585_dates_plot = ssp585_dates[model]
            for_global_sum_plot = for_global_sum[model]
            ssp585_global_sum_plot = ssp585_global_sum[model]
            for_global_mean_plot = for_global_mean[model]
            ssp585_global_mean_plot = ssp585_global_mean[model]
        ax1.plot(
                for_dates_plot,
                for_global_sum_plot,
                color=COLORS[INSTIT[model]],
                label=model,
                )
        ax1.plot(
                ssp585_dates_plot,
                ssp585_global_sum_plot,
                color=COLORS[INSTIT[model]],
                linestyle='dashed',
                )
        ax2.plot(
                for_dates_plot,
                for_global_mean_plot,
                color=COLORS[INSTIT[model]],
                )
        ax2.plot(
                for_dates_plot,
                for_global_mean_plot,
                color=COLORS[INSTIT[model]],
                linestyle='dashed',
                )
        print(model, "Delta F:", for_global_sum[model][-1] - ssp585_global_sum[model][-1])
    plt.xlim(left=for_dates[model][0], right=for_dates[model][-1])
    ax1.set_ylabel('Area (million km$^2$)')
    ax2.set_ylabel('Fraction (%)')
    plt.title('Tree area')
    ax1.legend(frameon=False)
    plt.savefig(f'{PLOTS_DIR}/tree_area_ssp126_all_models_global_sum.png', dpi=DPI)


def plot_global_mean_frac()->None:
    """Plot the global mean of tree fractions.
    """
    plt.figure()

    for model in MODELS:
        plt.plot(
                for_dates[model],
                for_global_mean[model],
                color=COLORS[INSTIT[model]],
                label=model,
                )
        plt.plot(
                ssp585_dates[model],
                ssp585_global_mean[model],
                color=COLORS[INSTIT[model]],
                linestyle='dashed',
                )
    plt.xlim(left=for_dates[model][0], right=for_dates[model][-1])
    plt.ylabel('Fraction %')
    plt.title('Tree fraction')
    plt.legend()
    plt.savefig(f'{PLOTS_DIR}/tree_fraction_ssp126_all_models_globalmean.png', dpi=DPI)


def plot_area_by_lat()->None:
    """Plot the tree area by latitude.
    """
    global for_dates, for_tree_area
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        lat_tree_area = for_tree_area[model].sum(axis=-1)
        plt.pcolormesh(
                for_dates[model],
                lats[model],
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
    cbar_ax = fig.add_axes([0.1,0.085,0.85,0.04])
    plt.colorbar(cax=cbar_ax, label='Area (million km$^2$)', orientation='horizontal', pad=0.05)
    plt.savefig(f'{PLOTS_DIR}/tree_area_ssp126_all_models_zonsum.png', dpi=DPI)


def plot_frac_by_lat()->None:
    """Plot the tree fraction by latitude.
    """
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        lat_tree_frac = np.ma.average(for_tree_frac[model], axis=-1)
        plt.pcolormesh(
                for_dates[model],
                lats[model],
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
    cbar_ax = fig.add_axes([0.1,0.085,0.85,0.04])
    plt.colorbar(cax=cbar_ax, label='Fraction %', orientation='horizontal', pad=0.05)
    plt.savefig(f'{PLOTS_DIR}/tree_area_ssp126_all_models_zonmean.png', dpi=DPI)


def plot_map_tree_area_change(indata:dict, name:str)->None:
    """Plot maps of the tree area change.
    """
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
        # pcolormesh expects +1 lon.
        llons = np.append(lons[model], [lons[model][-1] + lons[model][1] - lons[model][0]])
        # Shift plotting to centre grid, align with coastlines.
        llons = llons - (llons[1] - llons[0])/2
        llats = lats[model] - (lats[model][1] - lats[model][0])/4
        data = plt.pcolormesh(
                llons,
                llats,
                indata[model][-1] - indata[model][0],
                cmap=CMAP,
                vmin=-0.025,
                vmax=0.025,
                transform=ccrs.PlateCarree(),
                )
        axes[i].coastlines()
        plt.title(model, fontsize=8)
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, hspace=0.12, wspace=0.05)
    cbar_ax = fig.add_axes([0.1,0.085,0.85,0.04])
    plt.colorbar(
            data,
            ticks=np.arange(-90, 100, 10),
            cax=cbar_ax,
            label='Area (million km$^2$)',
            orientation='horizontal',
            pad=0.05,
            )
    plt.savefig(f'{PLOTS_DIR}/tree_area_map_2015-2100_{name}_delta.png', dpi=DPI)


def plot_map_tree_fraction_change(indata:dict, name:str)->None:
    """Plot maps of the tree fraction change.
    """
    fig, axes = plt.subplots(
            nrows=3,
            ncols=2,
            sharex=True,
            sharey=True,
            subplot_kw={'projection': ccrs.EckertIV()},
            )
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-95, 105, 10), ncolors=256)
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        # pcolormesh expects +1 lon.
        llons = np.append(lons[model], [lons[model][-1] + lons[model][1] - lons[model][0]])
        # Shift plotting to centre grid, align with coastlines.
        llons = llons - (llons[1] - llons[0])/2
        llats = lats[model] - (lats[model][1] - lats[model][0])/4
        data = plt.pcolormesh(
                llons,
                llats,
                indata[model][-1] - indata[model][0],
                norm=discrete_bins,
                cmap=CMAP,
                #vmin=-100,
                #vmax=100,
                transform=ccrs.PlateCarree(),
                )
        axes[i].coastlines()
        plt.title(model, fontsize=8)
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, hspace=0.12, wspace=0.05)
    cbar_ax = fig.add_axes([0.1,0.085,0.85,0.04])
    plt.colorbar(data,
            ticks=np.arange(-90, 100, 10),
            cax=cbar_ax,
            label='Tree fraction (%)',
            orientation='horizontal',
            pad=0.05,
            )
    plt.savefig(f'{PLOTS_DIR}/tree_frac_map_2015-2100_{name}_delta.png', dpi=DPI)


def plot_map_tree_fraction_difference(indata:dict, indata2:dict)->None:
    """Plot maps of the tree fraction change.
    """
    fig, axes = plt.subplots(
            nrows=3,
            ncols=2,
            sharex=True,
            sharey=True,
            subplot_kw={'projection': ccrs.EckertIV()},
            )
    fig.set_size_inches(9, 7)
    axes = axes.flatten()
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-95, 105, 10), ncolors=256)
    for i,model in enumerate(MODELS):
        plt.sca(axes[i])
        # pcolormesh expects +1 lon.
        llons = np.append(lons[model], [lons[model][-1] + lons[model][1] - lons[model][0]])
        # Shift plotting to centre grid, align with coastlines.
        llons = llons - (llons[1] - llons[0])/2
        llats = lats[model] - (lats[model][1] - lats[model][0])/4
        data = plt.pcolormesh(llons, llats, indata[model][-1] - indata2[model][-1],
                norm=discrete_bins,
                cmap=CMAP,
                #vmin=-100,
                #vmax=100,
                transform=ccrs.PlateCarree(),
                )
        axes[i].coastlines()
        plt.title(model, fontsize=8)
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, hspace=0.12, wspace=0.05)
    cbar_ax = fig.add_axes([0.1,0.085,0.85,0.04])
    plt.colorbar(data,
            ticks=np.arange(-90, 100, 10),
            cax=cbar_ax,
            label='Tree fraction (%)',
            orientation='horizontal',
            pad=0.05,
            )
    plt.savefig(f'{PLOTS_DIR}/tree_frac_map_2015-2100_diff.png', dpi=DPI)


def make_tree_frac_plots()->None:
    """Main function to plot all tree area and fraction plots.
    """
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
    #plot_global_mean_frac() # This is now plotted together in Figure 1

    # Plot area by latitude.
    plot_area_by_lat()

    # Plot frac by latitude.
    plot_frac_by_lat()

    # Plot global maps of tree area change.
    plot_map_tree_area_change(for_tree_area, 'for')
    plot_map_tree_area_change(ssp585_tree_area, 'ssp585')
    plot_map_tree_fraction_change(for_tree_frac, 'for')
    plot_map_tree_fraction_change(ssp585_tree_frac, 'ssp585')

    # Plot the maps of differences between simulations.
    plot_map_tree_fraction_difference(for_tree_frac, ssp585_tree_frac)


if __name__ != 'analysis.plot_model_trees':
    print("Creating tree fraction plots.")
    # Load data:
    for_tree_frac = {}
    for_tree_area = {}
    for_global_mean = {}
    for_global_sum = {}
    for_dates = {}
    ssp585_tree_frac = {}
    ssp585_tree_area = {}
    ssp585_global_mean = {}
    ssp585_global_sum = {}
    ssp585_dates = {}
    lons = {}
    lats = {}
    print("Loading:")
    for model in MODELS:
        print(f"    - {model}")
        print(f"        - Fractions")
        # Load the forestation experiment data.
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
        for_tree_frac[model], for_dates[model], _, _, = load_frac(files)

        # Load the SSP585 data.
        if model=='CanESM5':
            e = 'r1i1p1f1'
        else:
            e = ENSEMBLES[model]
        files = sorted(get_filenames(
                'C4MIP',
                INSTIT[model],
                model,
                'esm-ssp585',
                e,
                'Lmon',
                'treeFrac',
                ))
        ssp585_tree_frac[model], ssp585_dates[model], lats[model], lons[model] = load_frac(files)

        # Calculate global mean cover fractions.
        for_global_mean[model] = global_mean(for_tree_frac[model], lats[model])
        ssp585_global_mean[model] = global_mean(ssp585_tree_frac[model], lats[model])

        # Cacluate areas
        print("        - Areas")
        for_tree_area[model] = calc_tree_area(
                for_tree_frac[model],
                f'{DATA_DIR}/gridarea_{model}.nc',
                )
        ssp585_tree_area[model] = calc_tree_area(
                ssp585_tree_frac[model],
                f'{DATA_DIR}/gridarea_{model}.nc',
                )

        # Calculate global sum areas.
        print("        - Sums")
        for_global_sum[model] = np.nansum(for_tree_area[model], axis=(1,2))
        ssp585_global_sum[model] = np.nansum(ssp585_tree_area[model], axis=(1,2))

    make_tree_frac_plots()
    plt.show()
