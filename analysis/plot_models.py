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

if __name__ != 'analysis.plot_afforestation':
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filenames
    from constants import (CLIM_VARIABLES, ENSEMBLES, SEC_IN_DAY, TABLES,
                           VARIABLES, KG_IN_PG)
else:
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filenames
    from analysis.constants import (CLIM_VARIABLES, ENSEMBLES, SEC_IN_DAY,
                                    TABLES, VARIABLES, KG_IN_PG)
cdo = cdo_module.Cdo()

# Local constants
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
#load_npy_files = True # Uncomment to override previous check.

models = {
        'CCma':'CanESM5',
        'MIROC':'MIROC-ES2L',
        'MPI-M':'MPI-ESM1-2-LR'}

ensembles = {
        'CCma':'r1i1p2f1',
        'MIROC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1'}

data = {}
for inst in models.keys():

    filenames = get_filenames('LUMIP', inst, models[inst], 'esm-ssp585-ssp126Lu', ensembles[inst],
            'Lmon', 'cVeg')[0]

    cdo.gridarea(input=filenames, output=DATA_DIR+'/grid.nc')

    @cdod.cdo_cat(input2='')
    @cdod.cdo_mul(input2=DATA_DIR+'/grid.nc')
    @cdod.cdo_fldsum
    @cdod.cdo_yearmonmean
    @cdod.cdo_divc(str(KG_IN_PG))
    def cdo_pool_load_bad(var:str, input:str)->np.ma.MaskedArray:
        return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()

    filenames = ' '.join(get_filenames('LUMIP', inst, models[inst], 'esm-ssp585-ssp126Lu',
        ensembles[inst], 'Lmon', 'cVeg'))
    filenames = '[ '+filenames+' ]'
    data[inst] = cdo_pool_load_bad(input=filenames, var='cVeg')

years = list(range(2015, 2101))
plt.figure()
for m in models.keys():
    plt.plot(years, data[m]-data[m][0], label=models[m])
plt.ylabel('PgC')
plt.legend()
plt.show()

