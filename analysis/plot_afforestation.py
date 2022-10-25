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
from shapely.errors import ShapelyDeprecationWarning

if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cdo_calc_load import cdo_fetch_ensembles, cdo_load_anomaly_map
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
            PLOTS_DIR,
            SEC_IN_DAY,
            SEC_IN_YEAR,
            TABLES,
            VARIABLES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cdo_calc_load import cdo_fetch_ensembles, cdo_load_anomaly_map
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
            PLOTS_DIR,
            SEC_IN_DAY,
            SEC_IN_YEAR,
            TABLES,
            VARIABLES,
            )

warnings.filterwarnings(action='ignore', category=ShapelyDeprecationWarning)
cdo = Cdo()
cdo.debug = False

# Local constants
COLORS = {'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cLand':'maroon',
        'cSoil':'black',
        'tas':'black',
        'pr':'blue'}

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


# Load data for only afforested grid cells.
treeFrac = np.load(f'{DATA_DIR}/treeFrac_area_anomaly.npy')/M2_IN_MILKM2
NLAT = treeFrac.shape[0]
NLON = treeFrac.shape[1]


def reference_period(infile1:str, infile2:str, outfile:str, pyear:list=[2005, 2024]):
    """Extract a map of the mean over the reference period using CDO.
    The reference period spans 20 years centred on the start of the future simulation (2015), so
    [2005, 2024].
    """
    cdo.timmean(input=f'-selyear,{pyear[0]}/{pyear[1]} -cat '+infile1+' '+infile2, output=outfile)


def plot_ensembles(years:np.ndarray, data:np.ndarray, var:str, hlines:bool=True)->None:
    """Plot all ensemble members with ensemble mean and standard deviation.
    """
    model = 'ACCESS-ESM1.5'
    plt.figure()
    #for ens in ENSEMBLES:
    #    plt.plot(years, data[int(ens)-1], color='lightgray', linewidth=0.6, alpha=0.4)
    plt.plot(years, data.mean(axis=0), color=COLORS[var], label="Ensemble mean")
    plt.fill_between(
            years,
            data.min(axis=0),
            data.max(axis=0),
            color=COLORS[var],
            label="Ensemble range",
            alpha=0.4,
            )
    if hlines: plt.hlines(0, years[0], years[-1], color='black', linewidth=0.5)
    plt.xlim(left=years[0], right=years[-1])
    plt.xlabel('Year')
    if var[0]=='c': plt.ylabel(f'{var.upper()} anomaly [Pg(C)]')
    elif var=='tas': plt.ylabel(var.upper()+' anomaly ($^{\circ}$C)')
    elif var=='pr': plt.ylabel(f'{var.upper()} (mm/day)')
    else: plt.ylabel(f'{var.upper()} anomaly [Pg(C)/year]')
    plt.title(f"{model} {var.upper()}")


def plot_map(lons:np.ndarray, lats:np.ndarray, data:np.ndarray, var:str, label='')->None:
    """Plot a map of the difference between the start of the future period and the last year.
    """
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    max_abs = np.abs(data).max()
    plt.pcolormesh(lons, lats, data, cmap=jaisnb, vmax=max_abs, vmin=-max_abs,
            transform=ccrs.PlateCarree())
    ax.coastlines()
    if var in ['cVeg','cLitter','cSoil','cLand']:
        plt.colorbar(label='Pg(C)', orientation='horizontal', pad=0.05)
    elif var in ['pr']:
        plt.colorbar(label='mm/day', orientation='horizontal', pad=0.05)
    elif var in ['tas']:
        plt.colorbar(label='°C', orientation='horizontal', pad=0.05)
    else:
        plt.colorbar(label='Pg(C)/year', orientation='horizontal', pad=0.05)
    plt.title(var+' '+label)
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/global/{var}_ACCESS-ESM1.5_aff-esm-ssp585_{label}.png', dpi=DPI)


def make_veg_plots()->None:
    """Load vegetation data and run plotting routine.
    """
    if not os.path.exists(f'{PLOTS_DIR}/global'): os.mkdir(f'{PLOTS_DIR}/global')

    # Options.
    ## Recalculate the ensemble means from the CMORized CMIP6 data. If not, assume it's been done.
    recalculate_ens_mean = False

    # Local variables.
    global_sum_baselines = {}
    # Calculate the maps of base period means for each variable and ensemble member.
    print("\rCalculating baseline values.")
    for table in TABLES:
        for var in VARIABLES[table]:
            if recalculate_ens_mean:
                for ens in ENSEMBLES:
                    esmhist_files = get_filename('CMIP', 'esm-hist', ens, table, var)
                    aff_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
                    reference_period(
                            esmhist_files[-1],
                            aff_files[0],
                            f'{var}_ACCESS-ESM1-5_esm-hist-aff_r{ens}i1p1f1_200501-202412mean.nc'
                            )
                # The analysis for all ensemble members should be with respect to the same baseline
                # value. So I calculate the ensemble mean here.
                cdo.ensmean(
                        input=f'{var}_ACCESS-ESM1-5_esm-hist-aff_r*i1p1f1_200501-202412mean.nc',
                        output=f'{DATA_DIR}/' \
                                f'{var}_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc',
                        )
                # Clean up unneeded ensemble files.
                efiles = glob.glob(f'{var}_ACCESS-ESM1-5_esm-hist-aff_r*i1p1f1_200501-202412mean.nc')
                for f in efiles:
                    os.remove(f)
            # Calculate the global sum of each variable's ensemble mean. If a flux, convert to /year.
            # Need to also multiply by the grid cell area and land fraction.
            # LAND_FRAC_FILE is in % so must be divided by 100 first.
            # kg m-2 s-2 -> *land_frac*landarea*SEC_IN_YEAR/KG_IN_PG -> Pg/year
            if var in ['cVeg','cLitter','cSoil','cLand']:
                time_units = 1
            else:
                time_units = SEC_IN_YEAR
            global_sum_baselines[var] = cdo.divc(str(KG_IN_PG),
                    input=f'-mulc,{time_units} -fldsum -mul -mul ' \
                            f'{DATA_DIR}/' \
                            f'{var}_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc ' \
                            f'-divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
                    options='-L',
                    returnCdf=True).variables[var][:][0,0,0]

    print("Global sum baseline values:")
    print(global_sum_baselines)

    model = 'ACCESS-ESM1.5'
    for table in TABLES:
        for var in VARIABLES[table]:
            print(f"Processing {var}")
            # Load data
            if load_npy_files:
                aff_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy')
                ssp585_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy')
            else:
                aff_data = cdo_fetch_ensembles('LUMIP', 'esm-ssp585-ssp126Lu', table, var)
                ssp585_data = cdo_fetch_ensembles('C4MIP', 'esm-ssp585', table, var)
                np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy', aff_data)
                np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy', ssp585_data)

            # Anomaly, mean and standard deviation relative to 2015 baseline. Demonstrates overall
            # impact of climate, CO2 forcing and afforestation on pool/flux.
            data_anomaly = aff_data - global_sum_baselines[var]

            # Anomaly relative to the esm-ssp585 scenario from C4MIP. Demonstrates only the changes
            # due to afforestation.
            data_aff_diff = aff_data - ssp585_data

            # Plot the graphs for anomalies relative to 2015.
            years = list(range(2015, 2101))
            plot_ensembles(years, data_anomaly, var)
            plt.savefig(
                    f'{PLOTS_DIR}/global/{var}_{model}_esm-ssp585-ssp126Lu_ensembles_anomalies.png',
                    dpi=DPI,
                    )
            plt.close()

            # Plot the graphs for difference relative to the esm-ssp585
            plot_ensembles(years, data_aff_diff, var)
            plt.savefig(
                    f'{PLOTS_DIR}/global/{var}_{model}_esm-ssp585-ssp126Lu_ensembles_diff.png',
                    dpi=DPI,
                    )
            plt.close()


@cdod.cdo_cat(input2='')
@cdod.cdo_yearmonmean
def cdo_load_clim_map(var:str, input:str)->np.ma.MaskedArray:
    """Use cdo to load the first and last years of the future scenario and calculate the
    temporal anomaly between them.
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def make_clim_maps()->None:
    """Load climate data and run plotting routine.
    """
    if not os.path.exists(f'{PLOTS_DIR}/global'): os.mkdir(f'{PLOTS_DIR}/global')
    model = 'ACCESS-ESM1.5'
    for table in ['Amon',]:
        for var in CLIM_VARIABLES[table]:
            print(f"Processing {var}")
            # Load data
            if load_npy_files:
                diff_data = np.load(f'{DATA_DIR}/{var}_{model}_diff.npy')
                lats = np.load(f'{DATA_DIR}/lats.npy')
                lons = np.load(f'{DATA_DIR}/lons.npy')
            else:
                NTIMES = 86
                aff_data = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
                ssp585_data = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
                for e,ens in enumerate(ENSEMBLES):
                    aff_file = '[ '+ \
                            get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)[0]+' ]'
                    ssp585_file = '[ '+get_filename('C4MIP', 'esm-ssp585', ens, table, var)[0]+' ]'
                    aff_data[e,...] = cdo_load_clim_map(input=aff_file, var=var)
                    ssp585_data[e,...] = cdo_load_clim_map(input=ssp585_file, var=var)
                diff_data = aff_data - ssp585_data
                aff_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)[0]
                lats = nc.Dataset(aff_file, 'r').variables['lat'][:]
                lons = nc.Dataset(aff_file, 'r').variables['lon'][:]
                np.save(f'{DATA_DIR}/{var}_{model}_diff.npy',diff_data.data)
                np.save(f'{DATA_DIR}/lats.npy', lats.data)
                np.save(f'{DATA_DIR}/lons.npy', lons.data)

            # Calculate ensemble mean for the last 20 years
            diff_ens_mean = diff_data[:,-20:,...].mean(axis=(0,1))

            # Plot
            plot_map(lons, lats, diff_ens_mean, var, label='difference')


def make_clim_plots()->None:
    """Load climate data run plotting routine.
    """
    if not os.path.exists(f'{PLOTS_DIR}/global'): os.mkdir(f'{PLOTS_DIR}/global')
    model = 'ACCESS-ESM1.5'
    table = 'Amon'
    for var in CLIM_VARIABLES[table]:
        print(f"Processing {var}")
        if load_npy_files:
            aff_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy')
            ssp585_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy')
        else:
            aff_data = cdo_fetch_ensembles('LUMIP', 'esm-ssp585-ssp126Lu', table, var)
            ssp585_data = cdo_fetch_ensembles('C4MIP', 'esm-ssp585', table, var)
            np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy', aff_data)
            np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy', ssp585_data)

        # Calculate mean and standar deviation
        if var == 'tas':
            aff_data -= 273.15
            ssp585_data -= 273.15
        elif var=='pr':
            aff_data *= SEC_IN_DAY
            ssp585_data *= SEC_IN_DAY

        # Plot
        years = list(range(2015, 2101))
        plot_ensembles(years, aff_data, var, hlines=False)
        plt.savefig(f'{PLOTS_DIR}/global/{var}_{model}_esm-ssp585-ssp126Lu_ensembles.png', dpi=DPI)
        plt.close()
        plot_ensembles(years, ssp585_data, var, hlines=False)
        plt.savefig(f'{PLOTS_DIR}/global/{var}_{model}_esm-ssp585_ensembles.png', dpi=DPI)
        plt.close()
        plot_ensembles(years, ssp585_data-aff_data, var, hlines=True)
        plt.savefig(f'{PLOTS_DIR}/global/{var}_{model}_esm-ssp585_ensembles_diff.png', dpi=DPI)
        plt.close()


def make_veg_maps()->None:
    """Create maps of ensemble mean anomaly for the difference of the afforestation experiment
    and the esm-ssp585 scenario.
    """
    if not os.path.exists(f'{PLOTS_DIR}/global'): os.mkdir(f'{PLOTS_DIR}/global')
    for table in TABLES:
        for var in VARIABLES[table]:
            # Load the data.
            if load_npy_files:
                aff_data = np.load(f'{DATA_DIR}/{var}_aff_anomaly_maps.npy')
                ssp585_data = np.load(f'{DATA_DIR}/{var}_ssp585_anomaly_maps.npy')
                lats = np.load(f'{DATA_DIR}/lats.npy')
                lons = np.load(f'{DATA_DIR}/lons.npy')
            else:
                aff_data = np.ones((NENS,NLAT,NLON))*np.nan
                ssp585_data = np.ones((NENS,NLAT,NLON))*np.nan
                for e,ens in enumerate(ENSEMBLES):
                    aff_file = '[ ' + \
                            get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)[0]+' ]'
                    ssp585_file = '[ '+get_filename('C4MIP', 'esm-ssp585', ens, table, var)[0]+' ]'
                    aff_data[e,...] = cdo_load_anomaly_map(input=aff_file, var=var)
                    ssp585_data[e,...] = cdo_load_anomaly_map(input=ssp585_file, var=var)
                aff_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)[0]
                lats = nc.Dataset(aff_file, 'r').variables['lat'][:]
                lons = nc.Dataset(aff_file, 'r').variables['lon'][:]
                np.save(f'{DATA_DIR}/{var}_aff_anomaly_maps.npy', aff_data.data)
                np.save(f'{DATA_DIR}/{var}_ssp585_anomaly_maps.npy', ssp585_data.data)
                np.save(f'{DATA_DIR}/lats.npy', lats.data)
                np.save(f'{DATA_DIR}/lons.npy', lons.data)

            # Calculate ensemble mean of difference
            difference = aff_data - ssp585_data
            ensemble_mean_difference = difference.mean(axis=0)

            # Plot.
            plot_map(lons, lats, ensemble_mean_difference, var, label='difference')


@cdod.cdo_cat(input2='')
@cdod.cdo_yearmonmean
def cdo_clim_map_load(var:str, input:str)->np.ma.MaskedArray:
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def make_clim_aff_only()->None:
    """Make plots of climate variables for where there are afforesed gridcells only.
    """
    if not os.path.exists(f'{PLOTS_DIR}/global'): os.mkdir(f'{PLOTS_DIR}/global')
    # Load data for only afforested grid cells.
    treeFrac = np.load(f'{DATA_DIR}/treeFrac_area_anomaly.npy')/M2_IN_MILKM2
    NLAT = treeFrac.shape[0]
    NLON = treeFrac.shape[1]
    NENS = 10
    NTIMES = 86
    treeFrac = treeFrac*np.ones((NENS,NTIMES,NLAT,NLON))
    for var in CLIM_VARIABLES['Amon']:
        aff_data = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
        ssp585_data = np.ones((NENS,NTIMES,NLAT,NLON))*np.nan
        for e,ens in enumerate(ENSEMBLES):
            aff_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Amon', var)
            ssp585_file = get_filename('C4MIP', 'esm-ssp585', ens, 'Amon', var)
            aff_data[e,...] = cdo_clim_map_load(var=var, input=aff_file[0])
            ssp585_data[e,...] = cdo_clim_map_load(var=var, input=ssp585_file[0])

        # Calculate difference, mean and mask non afforested grid cells.
        clim_diff = aff_data - ssp585_data
        clim_diff[treeFrac<0.2] = np.nan
        lats = nc.Dataset(aff_file[0], 'r').variables['lat'][:]
        coslats = np.cos(lats[:,None]*np.ones((NLAT,NLON)))
        clim_diff = np.nansum(clim_diff*coslats, axis=(-1,-2))/np.sum(coslats)
        np.save(f'{DATA_DIR}/{var}_aff_only.np', clim_diff.data)

        # Plot
        years = list(range(2015, 2015 + NTIMES))
        plt.plot(years, clim_diff.mean(axis=0), color='black')
        for e in range(10):
            plt.plot(years, clim_diff[e,...], color='lightgray', alpha=0.5)
        plt.xlabel('Years')
        plt.ylabel('T')
        plt.show()


if __name__ != 'analysis.plot_afforestation':
    make_veg_plots()
    make_clim_plots()
    make_clim_maps()
    make_veg_maps()
    #make_clim_aff_only() # This didn't really work out.

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()

