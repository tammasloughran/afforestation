#!/usr/bin/env python3
"""Create histograms of temperature for particular regions in both experiments.
"""

import functools
import glob
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cdo_decorators as cdod
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from cdo import Cdo
import pdb

if __name__ == 'analysis.plot_regions':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import LAND_FRAC_FILE, get_filename
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
cdo.debug = True

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
        #'Amazonia': ([-19.63,12.70],[-81.81,-31.31]), # Forestation
        #'Eastern North America': ([24.94, 48.85],[-96.75,-51.87]), # Forestation & deforestaion
        #'Boreal North America': ([49.05,71.35],[-167.77,-53.81]), # Forestation
        #'Central Africa': ([-16.79,12.87],[-17.65,53.25]), # Low forestation
        #'Western Eruasia': ([46.21,60.23],[25.42,49.55]), # Deforrestation
        #'Boreal Eurasia': ([49.34,77.09],[50.9,175]), # Forestation
        #'East Asia': ([8.34,45.87],[96.25,148.87]), # Forestation and deforestation
        'Boreal Eurasia Gridpoint': ([63.74,63.76],[78.74,78.76]),
        'Central Africa Gridpoint': ([-7.6,-7.4],[18.74,18.76]),
        'Amazon Gridpoint': ([-11.26,-11.24],[309.374,309.376]), # Forestation ONLY gridpoint
        'Asia gridopint': ([29.75,30.25],[99,100]), # Forestation in the latter half of century.
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
load_npy_files = False # Uncomment to override previous check.

years = list(range(2015, 2101))


def is_not_nan(array):
    return np.logical_not(np.isnan(array))


def make_histogram_plots():
    model = 'ACCESS-ESM1-5'
    for region,box in REGIONS.items():


        # Create loader functions
        @cdod.cdo_selyear('2081/2100') # Select only the last 20 years of data.
        @cdod.cdo_selseason('JJA') # Select the summer season.
        @cdod.cdo_ifthen(input1=LAND_GT_50) # Mask for climate over land only. land frac > 50%
        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        def load_region_clim(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        @cdod.cdo_fldmean()
        @cdod.cdo_yearmonmean
        def load_region_frac(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


        print("Plotting", region)
        rname = region.replace(' ', '').lower()

        for_treeFrac = np.ones((NENS,len(years)))
        filenames = sorted(
                get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ENSEMBLES[0], 'Lmon', 'treeFrac'))
        for_treeFrac = load_region_frac(input=filenames[0], var='treeFrac')

        # Daily temperature variable
        table = 'day'
        for var in ['tas','tasmax','tasmin']:
            fig, ax = plt.subplots()
            print("Processing", var)
            # Load the data.
            if load_npy_files:
                for_data = np.load(f'{DATA_DIR}/{var}_daily_{model}_esm-ssp585-ssp126Lu_{rname}.npy')
                ssp585_data = np.load(f'{DATA_DIR}/{var}_daily_{model}_esm-ssp585_{rname}.npy')
            else:
                filenames = sorted(
                            get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ENSEMBLES[0], table, var))
                example_data = load_region_clim(input=filenames[-1], var=var)
                for_data = np.ones((NENS,)+example_data.shape)*np.nan
                ssp585_data = np.ones((NENS,)+example_data.shape)*np.nan
                for i, ens in enumerate(ENSEMBLES):
                    filenames = sorted(
                            get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var))
                    for_data[i,:] = load_region_clim(input=filenames[-1], var=var)
                    filenames = sorted(
                            get_filename('C4MIP', 'esm-ssp585', ens, table, var))
                    ssp585_data[i,:] = load_region_clim(input=filenames[-1], var=var)
                np.save(f'{DATA_DIR}/{var}_daily_{model}_esm-ssp585-ssp126Lu_{rname}.npy', for_data)
                np.save(f'{DATA_DIR}/{var}_daily_{model}_esm-ssp585_{rname}.npy', ssp585_data)

            for_data -= 273.15
            ssp585_data -= 273.15
            for_data[for_data>100] = np.nan
            ssp585_data[ssp585_data>100] = np.nan
            ks_test = stats.ks_2samp(for_data.flatten(), ssp585_data.flatten())
            print(f"{region} ks test p value is {ks_test.pvalue}")
            if ks_test.pvalue<0.05:
                print("The distributions are likely from different populations.")
            else:
                print("The distributions are likely the same.")

            n, bins, patches = plt.hist(
                    for_data.flatten()[is_not_nan(for_data.flatten())],
                    color='green',
                    bins=50,
                    density=True,
                    histtype='step',
                    label='esm-ssp585-ssp126Lu',
                    )
            plt.hist(
                    ssp585_data.flatten()[is_not_nan(ssp585_data.flatten())],
                    color='orange',
                    bins=bins,
                    density=True,
                    histtype='step',
                    label='esm-ssp585',
                    )
            #plt.annotate(f'p={ks_test.pvalue}', (0.1,0.9), xycoords='figure fraction')
            axin = ax.inset_axes([0.075,0.72,0.20,0.20])
            axin.plot(years,for_treeFrac)
            axin.set_title('Tree fraction')
            plt.xlabel('Temperature (°C)')
            plt.ylabel('Probablility')
            plt.title(region+f' {var}')
            plt.legend()
            plt.savefig(f'{PLOTS_DIR}/{var}_histogram_{rname}.png', dpi=DPI)


if __name__ != 'analysis.plot_regions':
    make_histogram_plots()

    plt.show()
