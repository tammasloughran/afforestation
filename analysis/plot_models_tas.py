
import glob
import os
import ipdb
import sys

import cdo as cdo_module
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import netCDF4 as nc
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
            NENS,
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
            NENS,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            )

cdo = cdo_module.Cdo()
cdo.debug = False

MODELS = { # ACCESS is excluded here. Needs separate plotting for ensembles.
        #'CSIRO':'ACCESS-ESM1-5',
        #'BCC':'BCC-CSM2-MR', # BCC has been excluded.
        #'NCC':'NorESM2-LM', # NorESM has been excluded.
        'CCma':'CanESM5',
        'NCAR':'CESM2',
        'NOAA-GFDL':'GFDL-ESM4',
        'MIROC':'MIROC-ES2L',
        'MPI-M':'MPI-ESM1-2-LR',
        'MOHC':'UKESM1-0-LL',
        }
ENSEMBLES = {
        #'BCC':'r1i1p1f1',
        #'NCC':'r1i1p1f1',
        'CCma':'r1i1p2f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCAR':'r1i1p1f1',
        'NOAA-GFDL':'r1i1p1f1',
        }
SSP585_ENSEMBLES = {
        #'BCC':'r1i1p1f1',
        #'NCC':'r1i1p1f1',
        'CCma':'r1i1p1f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCAR':'r1i1p1f1',
        'NOAA-GFDL':'r1i1p1f1',
        }


@cdod.cdo_yearmonmean
@cdod.cdo_selyear('2081/2100')
@cdod.cdo_timmean
def cdo_load_last(var:str, input:str)->np.ma.MaskedArray:
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def load_grid(f):
    ncfile = nc.Dataset(f, 'r')
    lons = ncfile.variables['lon'][:]
    lats = ncfile.variables['lat'][:]
    return lons,lats


def make_tas_plots():
    fig , axes = plt.subplots(
            nrows=3,
            ncols=2,
            #sharex=True,
            #sharey=True,
            subplot_kw={'projection': ccrs.EckertIV()},
            )
    fig.set_size_inches(9, 7)
    axes = axes.flatten()

    for_data = {}
    ssp585_data = {}
    i = 0
    for instit,model in MODELS.items():
        print(model)
        f = sorted(get_filenames(
                'LUMIP',
                instit,
                model,
                'esm-ssp585-ssp126Lu',
                ENSEMBLES[instit],
                'Amon',
                'tas',
                ))[-1]
        for_data[model] = cdo_load_last('tas', input=f)
        f = sorted(get_filenames(
                'C4MIP',
                instit,
                model,
                'esm-ssp585',
                SSP585_ENSEMBLES[instit],
                'Amon',
                'tas',
                ))[-1]
        ssp585_data[model] = cdo_load_last('tas', input=f)
        lons, lats = load_grid(f)

        diff = for_data[model] - ssp585_data[model]

        rng = 3
        plt.sca(axes[i])
        mappable = axes[i].pcolormesh(
                lons,
                lats,
                diff,
                vmin=-rng,
                vmax=rng,
                cmap='seismic',
                linewidth=0.2,
                edgecolors='face',
                transform=ccrs.PlateCarree(),
                )
        axes[i].coastlines()
        plt.title(model)

        i += 1
    plt.tight_layout()
    plt.sca(axes[-1])
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.1, right=0.95, hspace=0.12, wspace=0.05)
    cbar_ax = fig.add_axes([0.1,0.085,0.85,0.04])
    cbar = plt.colorbar(
            mappable,
            cax=cbar_ax,
            label='Temperature difference (°C)',
            orientation='horizontal',
            pad=0.05,
            )
    cbar.solids.set_edgecolor('face')
    plt.savefig(f'plots/models/tas_maps_models.png', dpi=DPI)


if __name__ != 'analysis.plot_models_tas':
    make_tas_plots()
    plt.show()
