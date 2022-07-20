#!/usr/bin/env python3
import glob
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

if __name__ == 'analysis.plot_australia':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cdo_calc_load import (load_aus_base_flux, load_aus_flux,
                                        load_aus_pool, load_aus_clim)
    from analysis.cmip_files import get_filename
    from analysis.constants import CLIM_VARIABLES, ENSEMBLES, TABLES, VARIABLES, SEC_IN_DAY
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cdo_calc_load import load_aus_base_flux, load_aus_flux, load_aus_pool, load_aus_clim
    from cmip_files import get_filename
    from constants import CLIM_VARIABLES, ENSEMBLES, TABLES, VARIABLES, SEC_IN_DAY

NTIMES = 86
NENS = 10
COLORS = {'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black',
        'tas':'black',
        'pr':'blue'}
REGIONS = { # 'Name': ([lat1,lat2],[lon1,lon2]), # afforestation/deforesation
        'Amazonia': ([-19.63,12.70],[-81.81,-31.31]), # reforrestation
        'Eastern North America': ([24.94, 48.85],[-96.75,-51.87]), # afforest & deforestaion
        'Boreal North America': ([49.05,71.35],[-167.77,-53.81]), # reforrestation
        'Central Africa': ([-16.79,12.87],[-17.65,53.25]), # low reforrestation
        'Western Eruasia': ([46.21,60.23],[25.42,49.55]), # deforrestation
        'Boreal Eurasia': ([49.34,77.09],[50.9,175]), # afforestation
        'East Asia': ([8.34,45.87],[96.25,148.87])} # afforestation and deforestation
PLOTS_DIR = './plots'
DATA_DIR = './data'

files = glob.glob('./*')
if PLOTS_DIR not in files:
    os.mkdir(PLOTS_DIR)
if DATA_DIR not in files:
    os.mkdir(DATA_DIR)

# Control flag
files = glob.glob('./data/*')
if any(['.npy' in f for f in files]):
    load_npy_files = True
else:
    load_npy_files = False
#load_npy_files = False # Uncomment to override previous check.

years = list(range(2015, 2101))


def plot_regions_map():
    """Plot the box regions on a map.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for region,box in REGIONS.items():
        plt.plot([box[1][0],box[1][1],box[1][1],box[1][0],box[1][0]],
                [box[0][1],box[0][1],box[0][0],box[0][0],box[0][1]],
                transform=ccrs.PlateCarree())
        plt.annotate(region, (box[1][0],box[0][0]))
    ax.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    plt.xlabel('Longitude (°E)')
    plt.ylabel('Latitude (°N)')
    plt.title('Regions')
    plt.tight_layout()
    plt.savefig('regional_analysis_map.png')
    plt.show()


def plot_veg_region(years, data, data_mean, data_std, var, label=''):
    plt.figure()
    for i,ens in enumerate(ENSEMBLES):
        plt.plot(years, data[i,...], color='grey', alpha=0.4)
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color=COLORS[var], alpha=0.4)
    plt.plot(years, data_mean, color=COLORS[var])
    plt.hlines(0, years[0], years[-1], linestyle='dotted', color='black')
    plt.xlim(left=years[0], right=years[-1])
    plt.ylabel('$\Delta$ Pg(C)')
    plt.xlabel('Time (Year)')
    plt.title(f"ACCES-ESM1.5 Australia {var}")
    plt.savefig('plots/'+var+'_aus_'+label+'.png')
    plt.show()


def plot_clim_region(years, data, data_mean, data_std, var, label=''):
    plt.figure()
    for i,ens in enumerate(ENSEMBLES):
        plt.plot(years, data[i,...], color='grey', alpha=0.4)
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color=COLORS[var], alpha=0.4)
    plt.plot(years, data_mean, color=COLORS[var])
    plt.xlim(left=years[0], right=years[-1])
    if var=='pr':
        plt.ylabel('mm/day')
    elif var=='tas':
        plt.ylabel('$^\circ$C')
    plt.xlabel('Time (Year)')
    plt.title(f"ACCES-ESM1.5 Australia {var}")
    plt.savefig('plots/'+var+'_aus_'+label+'.png')


def make_australia_plots():
    for table in TABLES:
        for var in VARIABLES[table]:
            # Get the reference period value.
            ref_file = '[ data/'+var+'_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc ]'
            if var in ['cVeg', 'cLitter', 'cSoil']:
                reference = load_aus_pool(input=ref_file, var=var)
            elif var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                reference = load_aus_base_flux(input=ref_file, var=var)

            # Load the afforestation data.
            if load_npy_files:
                aff_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy')
            else:
                aff_data = np.ones((NENS,NTIMES))*np.nan
                for i, ens in enumerate(ENSEMBLES):
                    filenames = ' '.join(
                            get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var))
                    filenames = '[ '+filenames+' ]'
                    if var in ['cVeg', 'cLitter', 'cSoil']:
                        aff_data[i,:] = load_aus_pool(input=filenames, var=var)
                    elif var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                        aff_data[i,:] = load_aus_flux(input=filenames, var=var)
                np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy', aff_data)

            # Calculate anomaly.
            anom_data = aff_data - reference

            # Load the ssp585 data.
            if load_npy_files:
                ssp585_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy')
            else:
                ssp585_data = np.ones((NENS,NTIMES))*np.nan
                for i, ens in enumerate(ENSEMBLES):
                    filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, table, var))
                    filenames = '[ '+filenames+' ]'
                    if var in ['cVeg', 'cLitter', 'cSoil']:
                        ssp585_data[i,:] = load_aus_pool(input=filenames, var=var)
                    elif var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                        ssp585_data[i,:] = load_aus_flux(input=filenames, var=var)
                np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy', ssp585_data)

            # Clalculate the difference.
            diff_data = aff_data - ssp585_data

            # Plot
            anom_data_mean = anom_data.mean(axis=0)
            anom_data_std = anom_data.std(axis=0, ddof=1)
            plot_veg_region(years, anom_data, anom_data_mean, anom_data_std, var, label='anom')
            diff_data_mean = diff_data.mean(axis=0)
            diff_data_std = diff_data.std(axis=0, ddof=1)
            plot_veg_region(years, diff_data, diff_data_mean, diff_data_std, var, label='diff')

    for var in ['tas', 'pr']:
        # Load the afforestation data
        if load_npy_files:
            aff_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy')
            ssp585_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy')
        else:
            aff_data = np.ones((NENS,NTIMES))*np.nan
            ssp585_data = np.ones((NENS,NTIMES))*np.nan
            for i, ens in enumerate(ENSEMBLES):
                filenames = ' '.join(get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Amon', var))
                filenames = '[ '+filenames+' ]'
                aff_data[i,:] = load_aus_clim(input=filenames, var=var)
                filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, 'Amon', var))
                filenames = '[ '+filenames+' ]'
                ssp585_data[i,:] = load_aus_clim(input=filenames, var=var)
            np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy', aff_data)
            np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy', ssp585_data)

        if var=='tas':
            aff_data -= 273.15
            ssp585_data -= 273.15
        elif var=='pr':
            aff_data *= SEC_IN_DAY
            ssp585_data *= SEC_IN_DAY

        # Calculate the difference
        diff_data = aff_data - ssp585_data

        # Plot
        plot_clim_region(years, aff_data, aff_data.mean(axis=0), aff_data.std(axis=0, ddof=1),
                var, label='aff_'+var)
        plot_clim_region(years, ssp585_data, ssp585_data.mean(axis=0),
                ssp585_data.std(axis=0, ddof=1), var, label='ssp585_'+var)
        plot_clim_region(years, diff_data, diff_data.mean(axis=0), diff_data.std(axis=0, ddof=1),
                var, label='diff_'+var)

    plt.show()


if __name__ != 'analysis.plot_australia':
    plot_regions_map()
    make_australia_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

