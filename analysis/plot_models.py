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
cdo.debug = True

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

models = {
        'BCC':'BCC-CSM2-MR',
        'CCma':'CanESM5',
        'MIROC':'MIROC-ES2L',
        'MOHC':'UKESM1-0-LL',
        'MPI-M':'MPI-ESM1-2-LR',
        'NCC':'NorESM2-LM',
        }
ensembles = {
        'BCC':'r1i1p1f1',
        'CCma':'r1i1p2f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCC':'r1i1p1f1',
        }
ssp585_ensembles = {
        'BCC':'r1i1p1f1',
        'CCma':'r1i1p1f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCC':'r1i1p1f1',
        }


VARIABLES = {
        'Lmon':[
            #'gpp',
            #'npp',
            #'ra',
            #'rh',
            #'nbp',
            'cVeg',
            'cLitter',
            ],
        'Emon':[
            'cSoil',
            'cLand',
            ],
        'Amon':[
            'pr',
            'tas',
            ],
        }


def make_model_plots():
    for table in VARIABLES.keys():
        for var in VARIABLES[table]:
            aff_data = {}
            ssp585_data = {}
            for instit in models.keys():
                if instit=='MOHC' and var=='cLitter': continue # Skip missing data.
                print(instit, models[instit])
                # Create the gridarea.nc file for this model.
                data_files = glob.glob(DATA_DIR)
                if 'gridarea_{models[instit]}.nc' not in data_files:
                    print('making grid area file')
                    afile = get_filenames(
                            'LUMIP',
                            instit,
                            models[instit],
                            'esm-ssp585-ssp126Lu',
                            ensembles[instit],
                            table,
                            var,
                            )[0]
                    cdo.gridarea(input=afile, output=f'{DATA_DIR}/gridarea_{models[instit]}.nc')


                # The loader function needs to be defined in this loop to account
                # for the model resolution.
                @cdod.cdo_cat(input2='')
                @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{models[instit]}.nc')
                @cdod.cdo_fldsum
                @cdod.cdo_yearmonmean
                @cdod.cdo_divc(str(KG_IN_PG))
                def cdo_pool_load_model(var:str, input:str)->np.ma.MaskedArray:
                    return cdo.copy(
                            input=input,
                            options='-L',
                            returnCdf=True,
                            ).variables[var][:].squeeze()


                model_land_frac = f'/g/data/p66/tfl561/CMIP6/C4MIP/{instit}/{models[instit]}' \
                        f'/esm-ssp585/{ssp585_ensembles[instit]}/fx/sftlf/gn/latest' \
                        f'/sftlf_fx_{models[instit]}_esm-ssp585_{ssp585_ensembles[instit]}_gn.nc'


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
                        models[instit],
                        'esm-ssp585-ssp126Lu',
                        ensembles[instit],
                        table,
                        var,
                        ))
                filenames = '[ ' + ' '.join(file_list) + ' ]'
                aff_data[instit] = loader(input=filenames, var=var)
                file_list = sorted(get_filenames(
                        'C4MIP',
                        instit,
                        models[instit],
                        'esm-ssp585',
                        ssp585_ensembles[instit],
                        table,
                        var,
                        ))
                filenames = '[ ' + ' '.join(file_list) + ' ]'
                ssp585_data[instit] = loader(input=filenames, var=var)

            # Plot the difference between simulations.
            plt.figure()
            # Some models have a slightly different time period.
            for m in models.keys():
                if m=='MOHC' and var=='cLitter':
                    continue # Skip missing data.
                if m=='MPI-M':
                    aff_data[m] = aff_data[m][:-1]
                if m=='NCC':
                    ssp585_data[m] = ssp585_data[m][:-1]
                years = list(range(2015, 2015 + aff_data[m].shape[0]))
                plt.plot(years, aff_data[m] - ssp585_data[m], label=models[m])
            # Load the ACCESS-ESM1-5 data.
            access_aff = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
            access_ssp585 = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_global.npy')
            diff = access_aff - access_ssp585
            ens_mean = np.mean(diff, axis=0)
            years = list(range(2015, 2015 + ens_mean.shape[0]))
            plt.plot(years, ens_mean, label='ACCESS-ESM1-5')
            for e in range(10):
                plt.plot(years, diff[e,:], color='gray', alpha=0.4)
            plt.ylabel('Pg(C)')
            plt.xlabel('Year')
            plt.title(f"{var} esm-ssp585-ssp126Lu - esm-ssp585")
            plt.legend()
            plt.savefig(f'{PLOTS_DIR}/{var}_model_intercomparison_diff.png', pdi=DPI)

            # Plot only the afforestation scenario.
            plt.figure()
            for m in models.keys():
                if m=='MOHC' and var=='cLitter': continue # Skip missing data.
                years = list(range(2015, 2015 + aff_data[m].shape[0]))
                plt.plot(years, aff_data[m], label=models[m])
            access_aff = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
            access_ssp585 = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_global.npy')
            diff = access_aff
            ens_mean = np.mean(diff, axis=0)
            years = list(range(2015, 2015 + ens_mean.shape[0]))
            plt.plot(years, ens_mean, label='ACCESS-ESM1-5')
            for e in range(10):
                plt.plot(years, diff[e,:], color='gray', alpha=0.4)
            plt.ylabel('Pg(C)')
            plt.xlabel('Year')
            plt.title(f"{var} esm-ssp585-ssp126Lu")
            plt.legend()
            plt.savefig(f'{PLOTS_DIR}/{var}_model_intercomparison_esm-ssp585-ssp126Lu.png', dpi=DPI)


if __name__ != 'analysis.plot_models':
    make_model_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()

