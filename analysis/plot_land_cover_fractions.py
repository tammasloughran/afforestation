#!/usr/bin/env python3
# plot_land_cover_fractions.py generates plots of the land cover fractions from the afforestation
# experiment and the esm-ssp585 experiment.
import glob
import os

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np

if __name__ == 'analysis.plot_land_cover_fractions':
    # imported as a module of the analysis package.
    from analysis.cdo_calc_load import cdo_cover_area_load, cdo_area_diff_load
    from analysis.cmip_files import get_filename, LAND_FRAC_FILE
    from analysis.constants import FRAC_VARIABLES, M2_IN_MILKM2, PLOTS_DIR, DATA_DIR
    from analysis.jaisnb import jaisnb
else:
    # is main program or imported as a module from another script.
    from cdo_calc_load import cdo_cover_area_load, cdo_area_diff_load
    from cmip_files import get_filename, LAND_FRAC_FILE
    from constants import FRAC_VARIABLES, M2_IN_MILKM2, PLOTS_DIR, DATA_DIR
    from jaisnb import jaisnb

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

# Control flag
files = glob.glob(f'{DATA_DIR}/*')
if any(['.npy' in f for f in files]):
    load_npy_files = True
else:
    load_npy_files = False
#load_npy_files = True # Uncomment to override previous check.


def make_land_cover_plot():
    """Create plot of the land cover for forest, grass and crop for esm-ssp585 and
    esm-ssp585-ssp126Lu.
    Plot the difference in land cover fractions/areas between 2015 and 2100 for the
    esm-ssp585-ssp126Lu experiment.
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
    CABLE_AREA = DATA_DIR+'/gridarea.nc'
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
    plt.savefig(f'{PLOTS_DIR}/CABLE_forests_deciduous_needle.png')


def make_area_anomaly_map():
    """Plot the change in area for trees, crop and grass in esm-ssp585-ssp126Lu between 2015
    and 2100.
    """
    for table in FRAC_VARIABLES.keys():
        for var in FRAC_VARIABLES[table]:
            frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', table, var)
            area_anomaly = cdo_area_diff_load(var=var, input=frac_file[0])/M2_IN_MILKM2
            ncfile = nc.Dataset(LAND_FRAC_FILE)
            mask = ncfile.variables['sftlf'][:]
            area_anomaly[mask==0] = np.nan
            lats = ncfile.variables['lat'][:]
            lons = ncfile.variables['lon'][:]

            np.save(f'{DATA_DIR}/{var}_area_anomaly.npy', area_anomaly.data)

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
            abs_max = np.nanmax(np.abs(area_anomaly))
            plt.pcolormesh(lons, lats, area_anomaly,
                    vmax=abs_max/2, vmin=-abs_max/2,
                    #cmap='nipy_spectral',#'PRGn', # nipy_spectral has a sharp color change at 0
                    cmap=jaisnb,
                    transform=ccrs.PlateCarree())
            ax.coastlines() # Drawing coastlines covers the coastal gridpoints.
            # Only needed for diverging colormaps that have white in the centre.
            plt.colorbar(label='Million km$^2$', orientation='horizontal', pad=0.05)
            plt.title(var.upper())
            plt.tight_layout()
            plt.savefig(f'{PLOTS_DIR}/'+var+'_anomaly.png')


if __name__ != 'analysis.plot_land_cover_fractions':
    make_land_cover_plot()
    make_afforestation_pft_plot()
    make_area_anomaly_map()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()

