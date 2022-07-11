#!/usr/bin/env python3
# plot_land_cover_fractions.py generates plots of the land cover fractions from the afforestation
# experiment and the esm-ssp585 experiment.
import glob
import os

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

if __name__ != 'analysis.plot_land_cover_fractions':
    from cdo_calc_load import cdo_cover_area_load
    from cmip_files import get_filename
    from constants import FRAC_VARIABLES, M2_IN_MILKM2
else:
    from analysis.cdo_calc_load import cdo_cover_area_load
    from analysis.cmip_files import get_filename
    from analysis.constants import FRAC_VARIABLES, M2_IN_MILKM2

COLORS = {
        'treeFrac':'green',
        'cropFrac': 'gold',
        'shrubFrac': 'olive',
        'grassFrac': 'lime'}
LABELS = {
        'AFF':
                {'treeFrac': 'tree AFF',
                'cropFrac': 'crop AFF',
                'shrubFrac': 'shrub AFF',
                'grassFrac': 'grass AFF'},
        'esm-ssp585':
                {'treeFrac': 'tree SSP5-8.5',
                'cropFrac': 'crop SSP5-8.5',
                'shrubFrac': 'shrub SSP5-8.5',
                'grassFrac': 'grass SSP5-8.5'}}
PLOTS_DIR = 'plots'


def make_land_cover_plot():
    """Create plot of the land cover for forest, grass and crop for esm-ssp585 and
    esm-ssp585-ssp126Lu.
    """
    table = 'Lmon'
    years = list(range(2015, 2101))

    aff_areas = {}
    ssp585_areas = {}
    for var in FRAC_VARIABLES[table]:
        frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', table, var)
        aff_areas[var] = cdo_cover_area_load(var, input=frac_file[0])/M2_IN_MILKM2
        plt.plot(years, aff_areas[var], color=COLORS[var], label=LABELS['AFF'][var])
        frac_file = get_filename('C4MIP', 'esm-ssp585', '1', table, var)
        ssp585_areas[var] = cdo_cover_area_load(var, input=frac_file[0])/M2_IN_MILKM2
        plt.plot(years, ssp585_areas[var], color=COLORS[var], linestyle='dashed',
                label=LABELS['esm-ssp585'][var])
    plt.xlim(left=years[0], right=years[-1])
    plt.ylabel("Area (million km$^2$)")
    plt.xlabel("Year")
    plt.legend()
    plt.savefig(f'{PLOTS_DIR}/land_areas.png')

    plt.figure()
    for var in FRAC_VARIABLES[table]:
        plt.plot(years, aff_areas[var]-ssp585_areas[var], color=COLORS[var],
                label=f'{LABELS["AFF"][var]}-SSP585')
    plt.xlim(left=years[0], right=years[-1])
    plt.ylabel("Area difference (million km$^2$)")
    plt.xlabel("Year")
    plt.hlines(y=0, xmin=years[0], xmax=years[-1], linestyle='dashed', color='black')
    plt.legend(loc='lower right')
    plt.savefig(f'{PLOTS_DIR}/land_areas_diff.png')

    print("By the end of the century:")
    print(f"Forrest expansion: {(aff_areas['treeFrac']-ssp585_areas['treeFrac'])[-1]}")
    print(f"Abandonment: {(aff_areas['cropFrac']-ssp585_areas['cropFrac'])[-1]}")
    print("")
    imin = np.argmin(aff_areas['cropFrac']-ssp585_areas['cropFrac'])
    print(f"Crop minimum: {(aff_areas['cropFrac']-ssp585_areas['cropFrac'])[imin]}")


def make_afforestation_pft_plot():
    """Plot the pfts for the afforestation scenario.
    """
    CABLE_LUC = '/g/data/p66/txz599/data/luc_ssp126/cableCMIP6_LC_2*.nc'
    CABLE_AREA = 'data/gridarea.nc'
    ncarea = nc.Dataset(CABLE_AREA, 'r')
    grid_area = ncarea.variables['cell_area'][:]
    files = glob.glob(CABLE_LUC)
    files.sort()
    nyears = len(files)
    fractions = np.ones((nyears,4,145,192))*np.nan
    for i,f in enumerate(files):
        ncfile = nc.Dataset(f, 'r')
        fractions[i,...] = ncfile.variables['fraction'][:,0:4,...]
    fractions[fractions>5] = np.nan # The input files do not have a missing value attribute.
    areas = fractions*grid_area/M2_IN_MILKM2
    global_sum = np.nansum(areas, axis=(2,3))
    years = list(range(2015,2102))
    plt.figure()
    plt.plot(years, global_sum[:,0,...], label='Evergreen needle leaf')
    plt.plot(years, global_sum[:,1,...], label='Evergreen broad leaf')
    plt.plot(years, global_sum[:,3,...], label='Deciduous broad leaf')
    plt.title("CABLE forests in SSP1-2.6")
    plt.xlabel('Time (year)')
    plt.ylabel('Area (million km$^2$)')
    plt.legend()
    plt.savefig('plots/CABLE_forests.png')
    plt.figure()
    plt.plot(years, global_sum[:,2,...], label='Deciduous needle leaf')
    plt.title("CABLE forests in SSP1-2.6")
    plt.xlabel('Time (year)')
    plt.ylabel('Area (million km$^2$)')
    plt.legend()
    plt.savefig('plots/CABLE_forests_deciduous_needle.png')


if __name__ != 'analysis.plot_land_cover_fractions':
    #make_land_cover_plot()
    make_land_cover_map()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()

