#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob
import os
import pdb
import sys

import cdo as cdo_module
import matplotlib.pyplot as plt
import numpy as np
import cdo_decorators as cdod
import pymannkendall as pmk
import ipdb

if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filenames
    from analysis.constants import (
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filenames
    from constants import (
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            )

cdo = cdo_module.Cdo()
cdo.debug = False

# Local constants
COLORS = {
        'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black',
        'tas':'black',
        'pr':'blue',
        }

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
#load_npy_files = True # Uncomment to override previous check.

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

MODELS = {
        'BCC':'BCC-CSM2-MR',
        'CCma':'CanESM5',
        'MIROC':'MIROC-ES2L',
        'MPI-M':'MPI-ESM1-2-LR',
        'NCC':'NorESM2-LM',
        'MOHC':'UKESM1-0-LL',
        }
COLORS = {
        'CSIRO':color_cycle[0],
        'BCC':color_cycle[1],
        'CCma':color_cycle[2],
        'MIROC':color_cycle[3],
        'MPI-M':color_cycle[4],
        'MOHC':color_cycle[5],
        'NCC':color_cycle[6],
        }
ENSEMBLES = {
        'BCC':'r1i1p1f1',
        'CCma':'r1i1p2f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCC':'r1i1p1f1',
        }
SSP585_ENSEMBLES = {
        'BCC':'r1i1p1f1',
        'CCma':'r1i1p1f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCC':'r1i1p1f1',
        }
VARIABLES = {
        #'Lmon':[
        #    #'gpp',
        #    #'npp',
        #    #'ra',
        #    #'rh',
        #    #'nbp',
        #    'cVeg',
        #    'cLitter',
        #    ],
        #'Emon':[
        #    'cSoil',
        #    'cLand',
        #    ],
        'Amon':[
            #'pr',
            'tas',
            ],
        }


def make_model_plots():
    if not os.path.exists(f'{PLOTS_DIR}/models'): os.mkdir(f'{PLOTS_DIR}/models')
    for table in VARIABLES.keys():
        for var in VARIABLES[table]:
            print(f'\nProcessing {var}:')
            aff_data = {}
            ssp585_data = {}
            for instit in MODELS.keys():
                if instit=='MOHC' and var=='cLitter': continue # Skip missing data.
                print('    -', instit, MODELS[instit])
                # Create the gridarea.nc file for this model.
                data_files = glob.glob(f'{DATA_DIR}/*')
                if f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc' not in data_files:
                    print('       - Making grid area file')
                    afile = get_filenames(
                            'LUMIP',
                            instit,
                            MODELS[instit],
                            'esm-ssp585-ssp126Lu',
                            ENSEMBLES[instit],
                            table,
                            var,
                            )[0]
                    cdo.gridarea(input=afile, output=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')


                # The loader function needs to be defined in this loop to account
                # for the model resolution.
                @cdod.cdo_cat(input2='')
                @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')
                @cdod.cdo_fldsum
                @cdod.cdo_yearmonmean
                @cdod.cdo_divc(str(KG_IN_PG))
                def cdo_pool_load_model(var:str, input:str)->np.ma.MaskedArray:
                    return cdo.copy(
                            input=input,
                            options='-L',
                            returnCdf=True,
                            ).variables[var][:].squeeze()


                model_land_frac = f'/g/data/p66/tfl561/CMIP6/C4MIP/{instit}/{MODELS[instit]}' \
                        f'/esm-ssp585/{SSP585_ENSEMBLES[instit]}/fx/sftlf/gn/latest' \
                        f'/sftlf_fx_{MODELS[instit]}_esm-ssp585_{SSP585_ENSEMBLES[instit]}_gn.nc'


                @cdod.cdo_cat(input2='')
                @cdod.cdo_mul(input2=model_land_frac)
                @cdod.cdo_fldmean()
                @cdod.cdo_yearmonmean
                def cdo_clim_load_model(var:str, input:str)->np.ma.MaskedArray:
                    return cdo.copy(
                            input=input,
                            options='-L',
                            returnCdf=True,
                            ).variables[var][:].squeeze()


                # Select pool or climate loader for this variable.
                if var in ['pr', 'tas']:
                    loader = cdo_clim_load_model
                else:
                    loader = cdo_pool_load_model

                # Load the CMIP6 models.
                file_list = sorted(get_filenames(
                        'LUMIP',
                        instit,
                        MODELS[instit],
                        'esm-ssp585-ssp126Lu',
                        ENSEMBLES[instit],
                        table,
                        var,
                        ))
                filenames = '[ ' + ' '.join(file_list) + ' ]'
                aff_data[instit] = loader(input=filenames, var=var)
                file_list = sorted(get_filenames(
                        'C4MIP',
                        instit,
                        MODELS[instit],
                        'esm-ssp585',
                        SSP585_ENSEMBLES[instit],
                        table,
                        var,
                        ))
                filenames = '[ ' + ' '.join(file_list) + ' ]'
                ssp585_data[instit] = loader(input=filenames, var=var)

            # Load the ACCESS-ESM1-5 data.
            access_aff = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
            access_ssp585 = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_global.npy')
            diff = access_aff - access_ssp585
            ens_mean = np.mean(diff, axis=0)

            # Some models have a slightly different time period.
            for m in MODELS.keys():
                if m=='MOHC' and var=='cLitter':
                    continue # Skip missing data.
                if m=='MPI-M':
                    aff_data[m] = aff_data[m][:-1]
                if m=='NCC':
                    ssp585_data[m] = ssp585_data[m][:-1]

            # Plot the difference between simulations for cpools.
            if not var=='tas' and not var=='pr':
                plt.figure()
                years = list(range(2015, 2015 + ens_mean.shape[0]))
                plt.plot(years, ens_mean, label='ACCESS-ESM1-5')
                for m in MODELS.keys():
                    diff_model = aff_data[m] - ssp585_data[m]
                    if m=='CCma' and var=='cVeg':
                        # CanESM5 cVeg has a large initial bias in cVeg. Remove this bias.
                        diff_model = diff_model - diff_model[0]
                    years = list(range(2015, 2015 + aff_data[m].shape[0]))
                    plt.plot(years, diff_model, color=COLORS[m], label=MODELS[m])
                for e in range(10):
                    plt.plot(years, diff[e,:], color='gray', alpha=0.4)
                if var=='tas':
                    plt.ylabel('Temperature (°C)')
                elif var=='pr':
                    plt.ylabel('Precipitation (mm/day)')
                else:
                    plt.ylabel('Pg(C)')
                plt.xlabel('Year')
                plt.title(f"{var} esm-ssp585-ssp126Lu - esm-ssp585")
                plt.legend()
                plt.savefig(f'{PLOTS_DIR}/models/{var}_model_intercomparison_diff.png', dpi=DPI)

            # Plot the trends for tas and pr.
            if var=='tas' or var=='pr':
                fig, axes = plt.subplots(nrows=4, ncols=2, sharex=True, sharey=True)
                axes = axes.flatten()
                years = list(range(2015, 2015 + ens_mean.shape[0]))
                axes[0].plot(years, ens_mean, label='ACCESS-ESM1-5')
                trend, h, p, z, tau, s, var_s, slope, intercept = pmk.original_test(
                        diff.mean(axis=0), # ACCESS-ESM1.5
                        alpha=0.05,
                        )
                print(f'    - ACCESS-ESM1-5 trend={trend}, p={p}, h={h}')
                trend_line = slope*np.arange(len(years)) + intercept
                if h:
                    axes[0].plot(years, trend_line, color=COLORS['CSIRO'])
                else:
                    axes[0].plot(years, trend_line, color=COLORS['CSIRO'], linestyle='dotted')
                for i,m in enumerate(MODELS.keys()):
                    diff_model = aff_data[m] - ssp585_data[m]
                    years = list(range(2015, 2015 + aff_data[m].shape[0]))
                    axes[i+1].plot(years, diff_model, color=COLORS[m], label=MODELS[m])
                    trend, h, p, z, tau, s, var_s, slope, intercept = pmk.original_test(
                            diff_model,
                            alpha=0.05,
                            )
                    print(f'    - {MODELS[m]} trend={trend}, p={p}, h={h}')
                    trend_line = slope*np.arange(len(years)) + intercept
                    if h: # Hypothesis that there exists a trend is true.
                        axes[i+1].plot(years, trend_line, color=COLORS[m])
                    else:
                        axes[i+1].plot(years, trend_line, color=COLORS[m],
                                linestyle='dotted')
                axes[7].set_axis_off()
                fig.add_subplot(111, frameon=False)
                plt.tick_params(
                        labelcolor='none',
                        which='both',
                        top=False,
                        bottom=False,
                        left=False,
                        right=False)
                plt.plot([0], [0], color=COLORS['CSIRO'], label='ACCESS-ESM1-5')
                for m in MODELS.keys():
                    plt.plot([0], [0], color=COLORS[m], label=MODELS[m])
                plt.legend(loc='lower right')
                plt.xlabel('Year')
                if var=='tas':
                    plt.ylabel('Temperature (°C)')
                else:
                    plt.ylabel('Precipitation (mm/day)')
                fig.suptitle(f'{var} global mean trends')
                plt.tight_layout()
                plt.savefig(f'{PLOTS_DIR}/models/{var}_trends.png', dpi=DPI)

            # Plot only the afforestation scenario for absolute values.
            plt.figure()
            access_aff = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
            access_ssp585 = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_global.npy')
            diff = access_aff
            ens_mean = np.mean(diff, axis=0)
            years = list(range(2015, 2015 + ens_mean.shape[0]))
            plt.plot(years, ens_mean, color=COLORS['CSIRO'], label='ACCESS-ESM1-5')
            for m in MODELS.keys():
                if m=='MOHC' and var=='cLitter': continue # Skip missing data.
                years = list(range(2015, 2015 + aff_data[m].shape[0]))
                plt.plot(years, aff_data[m], color=COLORS[m], label=MODELS[m])
            for e in range(10):
                plt.plot(years, diff[e,:], color='gray', alpha=0.4)
            if var=='tas':
                plt.ylabel('Temperature (°C)')
            elif var=='pr':
                plt.ylabel('Precipitation (mm/day)')
            else:
                plt.ylabel('Pg(C)')
            plt.xlabel('Year')
            plt.title(f"{var} esm-ssp585-ssp126Lu")
            plt.legend()
            plt.savefig(
                    f'{PLOTS_DIR}/models/{var}_model_intercomparison_esm-ssp585-ssp126Lu.png',
                    dpi=DPI,
                    )


if __name__ != 'analysis.plot_models':
    make_model_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()

