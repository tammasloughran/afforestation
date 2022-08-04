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
            CLIM_VARIABLES,
            DATA_DIR,
            ENSEMBLES,
            KG_IN_PG,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            VARIABLES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filenames
    from constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            ENSEMBLES,
            KG_IN_PG,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            VARIABLES,
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
        'CCma':'CanESM5',
        'MIROC':'MIROC-ES2L',
        'MOHC':'UKESM1-0-LL',
        'MPI-M':'MPI-ESM1-2-LR',
        'NCC':'NorESM2-LM',
        }
ensembles = {
        'CCma':'r1i1p2f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCC':'r1i1p1f1',
        }
ssp585_ensembles = {
        'CCma':'r1i1p1f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCC':'r1i1p1f1',
        }

aff_data = {}
ssp585_data = {}
for inst in models.keys():
    print(inst, models[inst])
    # Create the gridarea.nc file.
    data_files = glob.glob(DATA_DIR)
    if 'gridarea_{models[inst]}.nc' not in data_files:
        print('making grid area file')
        afile = get_filenames(
                'LUMIP',
                inst,
                models[inst],
                'esm-ssp585-ssp126Lu',
                ensembles[inst],
                'Lmon',
                'cVeg',
                )[0]
        cdo.gridarea(input=afile, output=f'{DATA_DIR}/gridarea_{models[inst]}.nc')


    # The loader function needs to be defined in this loop to account
    # for the model resolution.
    @cdod.cdo_cat(input2='')
    @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{models[inst]}.nc')
    @cdod.cdo_fldsum
    @cdod.cdo_yearmonmean
    @cdod.cdo_divc(str(KG_IN_PG))
    def cdo_pool_load_model(var:str, input:str)->np.ma.MaskedArray:
        return cdo.copy(
                input=input,
                options='-L',
                returnCdf=True,
                ).variables[var][:].squeeze()


    file_list = sorted(get_filenames(
            'LUMIP',
            inst,
            models[inst],
            'esm-ssp585-ssp126Lu',
            ensembles[inst],
            'Lmon',
            'cVeg',
            ))
    filenames = '[ ' + ' '.join(file_list) + ' ]'
    aff_data[inst] = cdo_pool_load_model(input=filenames, var='cVeg')
    file_list = sorted(get_filenames(
            'C4MIP',
            inst,
            models[inst],
            'esm-ssp585',
            ssp585_ensembles[inst],
            'Lmon',
            'cVeg',
            ))
    filenames = '[ ' + ' '.join(file_list) + ' ]'
    ssp585_data[inst] = cdo_pool_load_model(input=filenames, var='cVeg')


# Plot.
plt.figure()
for m in models.keys():
    if m=='MPI-M':
        aff_data[m] = aff_data[m][:-1]
    if m=='NCC':
        ssp585_data[m] = ssp585_data[m][:-1]
    years = list(range(2015, 2015 + aff_data[m].shape[0]))
    plt.plot(years, aff_data[m] - ssp585_data[m], label=models[m])
access_aff = np.load(f'{DATA_DIR}/cVeg_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
access_ssp585 = np.load(f'{DATA_DIR}/cVeg_ACCESS-ESM1.5_esm-ssp585_global.npy')
diff = access_aff - access_ssp585
ens_mean = np.mean(diff, axis=0)
years = list(range(2015, 2015 + ens_mean.shape[0]))
plt.plot(years, ens_mean, label='ACCESS-ESM1-5')
for e in range(10):
    plt.plot(years, diff[e,:], color='gray', alpha=0.4)
plt.ylabel('Pg(C)')
plt.xlabel('Year')
plt.title("cVeg esm-ssp585-ssp126Lu - esm-ssp585")
plt.legend()
plt.savefig(f'{PLOTS_DIR}/cVeg_model_intercomparison.png')
plt.show()

