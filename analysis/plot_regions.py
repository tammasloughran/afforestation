#!/usr/bin/env python3
import functools
import glob
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cdo_decorators as cdod
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pymannkendall as pmk
from cdo import Cdo
import pdb

if __name__ == 'analysis.plot_regions':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import LAND_FRAC_FILE, get_filename
    from analysis.cdo_calc_load import cdo_mul_land_area
    from analysis.constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
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
    from cmip_files import LAND_FRAC_FILE, get_filename
    from cdo_calc_load import cdo_mul_land_area
    from constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            NENS,
            NTIMES,
            PLOTS_DIR,
            SEC_IN_DAY,
            SEC_IN_YEAR,
            TABLES,
            VARIABLES,
            )

cdo = Cdo()
cdo.debug = False

COLORS = {
        'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black',
        'cLand':'maroon',
        'tas':'black',
        'pr':'blue',
        }
REGIONS = { # 'Name': ([lat1,lat2],[lon1,lon2]), # Forestation/deforesation
        #'Australia': ([-45,-10],[110,155]), # Neutral
        'Amazonia': ([-19.63,12.70],[-81.81,-31.31]), # Forestation
        'Eastern North America': ([24.94, 48.85],[-96.75,-51.87]), # Forestation & deforestaion
        'Boreal North America': ([49.05,71.35],[-167.77,-53.81]), # Forestation
        'Central Africa': ([-16.79,12.87],[-17.65,53.25]), # Low forestation
        #'Western Eruasia': ([46.21,60.23],[25.42,49.55]), # Deforrestation
        'Boreal Eurasia': ([49.34,77.09],[50.9,175]), # Forestation
        'East Asia': ([8.34,45.87],[96.25,148.87]), # Forestation and deforestation
        #'Boreal Eurasia Gridpoint': ([63.74,63.76],[78.74,78.76]),
        #'Central Africa Gridpoint': ([-7.6,-7.4],[18.74,18.76]),
        'Amazon Gridpoint': ([-11.26,-11.24],[309.374,309.376]), # Forestation ONLY gridpoint
        'Asia Gridpoint': ([29.75,30.25],[99,100]), # Forestation in the latter half of century.
        'North America Gridpoint': ([37,38],[273,274]), # Eastern North america, forestation but cooling.
        }

LAND_GT_50 = f'{DATA_DIR}/land_gt_50.nc'

files = glob.glob('./*')
if PLOTS_DIR not in files:
    os.mkdir(PLOTS_DIR)
if DATA_DIR not in files:
    os.mkdir(DATA_DIR)

# Control flag.
files = glob.glob('./data/*')
if any(['.npy' in f for f in files]):
    load_npy_files = True
else:
    load_npy_files = False
load_npy_files = True # Uncomment to override previous check.

years = list(range(2015, 2101))


def plot_regions_map()->None:
    """Plot the box regions on a map.
    """
    if not os.path.exists(f'{PLOTS_DIR}/regional'): os.mkdir(f'{PLOTS_DIR}/regional')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for region,box in REGIONS.items():
        if 'point' in region:
            plt.scatter(box[1][0], box[0][1])
        else:
            plt.plot(
                    [box[1][0],box[1][1],box[1][1],box[1][0],box[1][0]],
                    [box[0][1],box[0][1],box[0][0],box[0][0],box[0][1]],
                    transform=ccrs.PlateCarree(),
                    )
            plt.annotate(region, (box[1][0],box[0][0]))
    ax.set_xticks([-180,-120,-60,0,60,120,180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90,-60,-30,0,30,60,90], crs=ccrs.PlateCarree())
    plt.xlabel('Longitude (°E)')
    plt.ylabel('Latitude (°N)')
    plt.title('Regions')
    #plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/regional/regional_analysis_map.png', dpi=DPI)
    plt.close()


def plot_veg_region(
            years:np.ndarray,
            data:np.ndarray,
            var:str,
            region:str,
            label:str='')->None:
    """Plot vegetation data for regional analysis.
    """
    if not os.path.exists(f'{PLOTS_DIR}/regional'): os.mkdir(f'{PLOTS_DIR}/regional')
    plt.figure()
    #for i,ens in enumerate(ENSEMBLES):
    #    plt.plot(years, data[i,...], color='grey', alpha=0.4)
    plt.fill_between(years, data.min(axis=0), data.max(axis=0), color=COLORS[var], alpha=0.5)
    plt.plot(years, data.mean(axis=0), color=COLORS[var])
    plt.hlines(0, years[0], years[-1], linestyle='dotted', color='black')
    plt.xlim(left=years[0], right=years[-1])
    plt.ylabel('$\Delta$ Pg(C)')
    plt.xlabel('Time (Year)')
    plt.title(f"ACCESE-ESM1.5 {region} {var}")
    reg = region.replace(' ', '').lower()
    plt.savefig(f'{PLOTS_DIR}/regional/{var}_{reg}_{label}.png', dpi=DPI)
    plt.close()


def plot_clim_region(
            years:np.ndarray,
            data:np.ndarray,
            var:str,
            region:str,
            label:str='')->None:
    """Plot climate data for regional analysis.
    """
    plt.figure()
    #for i,ens in enumerate(ENSEMBLES):
    #    plt.plot(years, data[i,...], color='grey', alpha=0.4)
    plt.fill_between(years, data.min(axis=0), data.max(axis=0), color=COLORS[var], alpha=0.5)
    plt.plot(years, data.mean(axis=0), color=COLORS[var])
    if label=='diff':
        plt.hlines(0, years[0], years[-1], linewidth=0.5, color='black')
    plt.xlim(left=years[0], right=years[-1])
    if var=='pr':
        plt.ylabel('mm/day')
    elif var=='tas':
        plt.ylabel('$^\circ$C')
    plt.xlabel('Time (Year)')
    plt.title(f"ACCESS-ESM1.5 {region} {var}")
    reg = region.replace(' ', '').lower()
    plt.savefig(f'{PLOTS_DIR}/regional/{var}_{reg}_{label}.png', dpi=DPI)
    plt.close()


def make_regional_plots()->None:
    """Load data and generate regional analysis plots.
    """
    if not os.path.exists(f'{PLOTS_DIR}/regional'): os.mkdir(f'{PLOTS_DIR}/regional')
    stats_file = open('regions_stats.md', 'w')
    model = 'ACCESS-ESM1-5'
    for region,box in REGIONS.items():
        # Create loader functions
        @cdod.cdo_cat(input2='') # Concatenate all files in input1.
        @cdo_mul_land_area
        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        @cdod.cdo_fldsum
        @cdod.cdo_yearmonmean
        @cdod.cdo_divc(str(KG_IN_PG))
        def load_region_pool(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, returnCdf=True, options='-L').variables[var][:].squeeze()


        @cdod.cdo_cat(input2='')
        @cdo_mul_land_area
        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        @cdod.cdo_fldsum
        @cdod.cdo_mulc(str(SEC_IN_YEAR))
        @cdod.cdo_divc(str(KG_IN_PG))
        def load_region_base_flux(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


        @cdod.cdo_cat(input2='')
        @cdo_mul_land_area
        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        @cdod.cdo_fldsum
        @cdod.cdo_mulc(str(SEC_IN_DAY))
        @cdod.cdo_muldpm
        @cdod.cdo_yearsum
        @cdod.cdo_divc(str(KG_IN_PG))
        def load_region_flux(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


        @cdod.cdo_cat(input2='') # Concatenate all files in input1.
        @cdod.cdo_ifthen(input1=LAND_GT_50) # Mask for climate over land only. land frac > 50%
        # Alternatively I can multiply by the land fraction and divide by 100. Doesnt work for temp.
        #@cdod.cdo_mul(input2=LAND_FRAC_FILE) # Mask for climate over land only.
        #@cdod.cdo_divc(str(100)) # LAND_FRAC_FILE is in %.
        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        @cdod.cdo_fldmean(weights='TRUE') # Area weighted spatial mean.
        @cdod.cdo_yearmonmean # Annual mean
        def load_region_clim(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


        print("Plotting", region)
        rname = region.replace(' ', '').lower()

        # Vegetation variables
        for table in TABLES:
            for var in VARIABLES[table]:
                print("    - Processing", var)
                # Get the reference period value.
                ref_file = f'[ {DATA_DIR}/{var}_{model}_esm-hist-aff_ensmean_200501-202412mean.nc ]'
                if var in ['cVeg','cLitter','cSoil','cLand']:
                    reference = load_region_pool(input=ref_file, var=var)
                elif var in ['gpp','npp','ra','rh','nbp']:
                    reference = load_region_base_flux(input=ref_file, var=var)

                # Load the afforestation data.
                if load_npy_files:
                    aff_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_{rname}.npy')
                else:
                    aff_data = np.ones((NENS,NTIMES))*np.nan
                    for i, ens in enumerate(ENSEMBLES):
                        filenames = ' '.join(
                                get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var))
                        filenames = '[ '+filenames+' ]'
                        if var in ['cVeg','cLitter','cSoil','cLand']:
                            aff_data[i,:] = load_region_pool(input=filenames, var=var)
                        elif var in ['gpp','npp','ra', 'rh','nbp']:
                            aff_data[i,:] = load_region_flux(input=filenames, var=var)
                    np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_{rname}.npy', aff_data)

                # Calculate anomaly.
                anom_data = aff_data - reference

                # Load the ssp585 data.
                if load_npy_files:
                    ssp585_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585_{rname}.npy')
                else:
                    ssp585_data = np.ones((NENS,NTIMES))*np.nan
                    for i, ens in enumerate(ENSEMBLES):
                        filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, table, var))
                        filenames = '[ '+filenames+' ]'
                        if var in ['cVeg','cLitter','cSoil','cLand']:
                            ssp585_data[i,:] = load_region_pool(input=filenames, var=var)
                        elif var in ['gpp','npp','ra','rh','nbp']:
                            ssp585_data[i,:] = load_region_flux(input=filenames, var=var)
                    np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585_{rname}.npy', ssp585_data)

                # Clalculate the difference.
                diff_data = aff_data - ssp585_data

                # Plot.
                plot_veg_region(
                        years,
                        anom_data,
                        var,
                        region,
                        label='anom',
                        )
                plot_veg_region(
                        years,
                        diff_data,
                        var,
                        region,
                        label='diff',
                        )

        # Climate Variables
        table = 'Amon'
        for var in CLIM_VARIABLES[table]:
            print("    - Processing", var)
            # Load the afforestation data.
            if load_npy_files:
                aff_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_{rname}.npy')
                ssp585_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585_{rname}.npy')
            else:
                aff_data = np.ones((NENS,NTIMES))*np.nan
                ssp585_data = np.ones((NENS,NTIMES))*np.nan
                for i, ens in enumerate(ENSEMBLES):
                    filenames = ' '.join(get_filename(
                        'LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var))
                    filenames = '[ '+filenames+' ]'
                    aff_data[i,:] = load_region_clim(input=filenames, var=var)
                    filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, table, var))
                    filenames = '[ '+filenames+' ]'
                    ssp585_data[i,:] = load_region_clim(input=filenames, var=var)
                np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_{rname}.npy', aff_data)
                np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585_{rname}.npy', ssp585_data)

            if var=='tas':
                aff_data -= 273.15
                ssp585_data -= 273.15
            elif var=='pr':
                aff_data *= SEC_IN_DAY
                ssp585_data *= SEC_IN_DAY

            # Calculate the difference.
            diff_data = aff_data - ssp585_data

            # Print statistics to file.
            diff_mean = diff_data.mean(axis=0)
            stats_file.write(f'# {var} for {region}\n')
            stats_file.write('## Autocorrelation at lag 1\n')
            auto_corr, p = stats.spearmanr(diff_mean[:-1], diff_mean[1:])
            stats_file.write(f'r={auto_corr}, p={p}\n')
            stats_file.write('## Trend\n')
            # trend: tells the trend (increasing, decreasing or no trend)
            # h: True (if trend is present) or False (if trend is absent)
            # p: p-value of the significance test
            # z: normalized test statistics
            # Tau: Kendall Tau
            # s: Mann-Kendal's score
            # var_s: Variance S
            # slope: Theil-Sen estimator/slope
            # intercept: intercept of Kendall-Theil Robust Line
            trend, h, p, z, tau, s, var_s, slope, intercept = pmk.original_test(
                    diff_mean,
                    alpha=0.05,
                    )
            stats_file.write(f'trend={trend}, p={p}, h={h}\n')
            trend_line = slope*np.arange(NTIMES) + intercept

            # Plot.
            plot_clim_region(
                    years,
                    aff_data,
                    var,
                    region,
                    label='aff',
                    )
            plot_clim_region(
                    years,
                    diff_data,
                    var,
                    region,
                    label='diff',
                    )
            # I don't really need to plot the ssp585 simulation, but it's here just in case.
            #plot_clim_region(
            #        years,
            #        ssp585_data,
            #        var,
            #        region,
            #        label='ssp585',
            #        )
    stats_file.close()


if __name__ != 'analysis.plot_regions':
    plot_regions_map()
    #make_regional_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()
