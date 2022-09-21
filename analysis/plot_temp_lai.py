#!/usr/bin/env python3
import glob
import os

import cdo_decorators as cdod
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from cdo import Cdo
import pdb

if __name__ == 'analysis.plot_regions':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.constants import (
            DATA_DIR,
            DPI,
            NTIMES,
            PLOTS_DIR,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from constants import (
            DATA_DIR,
            DPI,
            NTIMES,
            PLOTS_DIR,
            )

cdo = Cdo()
cdo.debug = True

REGIONS = { # 'Name': ([lat1,lat2],[lon1,lon2]), # Forestation/deforesation
        'Asia Gridopint': ([29.75,30.25],[99,100]), # Forestation in the latter half of century.
        'Central Africa Gridpoint': ([-7.6,-7.4],[18.74,18.76]),
        'Boreal Eurasia Gridpoint': ([63.74,63.76],[78.74,78.76]),
        'Amazon Gridpoint': ([-11.26,-11.24],[309.374,309.376]), # Forestation ONLY gridpoint
        }
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
PFTS = {
        1:'Evergreen Needleleaf',
        2:'Evergreen Broadleaf',
        3:'Deciduous Needleleaf',
        4:'Deciduous Broadleaf',
        5:'Shrub',
        6:'C3 Grass',
        7:'C4 Grass',
        8:'Tundra',
        9:'C3 Crop',
        10:'C4 Crop',
        11:'Wetland',
        12:'',
        13:'',
        14:'Barren',
        15:'Urban',
        16:'Lakes',
        17:'Ice',
        }
EXPERIMENTS = {
        'esm-ssp585':'SSP-EDC-585',
        'esm-ssp585-ssp126Lu':'SSP-585-126',
        }
ENSEMBLES = ['03','04','06','07','08','09','10','11','12','13']
NENS = len(ENSEMBLES)

years = np.arange(2015, 2101)

files = glob.glob('./*')
if PLOTS_DIR not in files:
    os.mkdir(PLOTS_DIR)
if DATA_DIR not in files:
    os.mkdir(DATA_DIR)

def make_temp_lai_plots()->None:
    # Load Data

    for region, box in REGIONS.items():
        @cdod.cdo_sellonlatbox(str(box[1][0]), str(box[1][1]), str(box[0][0]), str(box[0][1]))
        #@cdod.cdo_fldmean()
        def load_region_lai(var:str, input:str)->np.ma.MaskedArray:
            return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


        for ens in ENSEMBLES:
            for_lai = load_region_lai(
                    input=f'{ARCHIVE_DIR}/esm-ssp585-ssp126Lu/lai_pfts_SSP-585-126-{ens}.nc',
                    var='lai_pfts',
                    )
            ssp585_lai = load_region_lai(
                    input=f'{ARCHIVE_DIR}/esm-ssp585/lai_pfts_SSP-EDC-585-{ens}.nc',
                    var='lai_pfts',
                    )
            for pp in [0,1,3,5,7]:
                plt.figure()
                plt.plot(for_lai[:,pp], color='green', label='esm-ssp585-ssp126Lu')
                plt.plot(ssp585_lai[:,pp], color='orange', label='esm-ssp585')
                plt.title(f'{region} {PFTS[pp+1]} LAI')
                plt.legend()
            plt.show()
            pdb.set_trace()


    # Plot data

if __name__ != 'analysis.plot_regions':
    make_temp_lai_plots()

    plt.show()

