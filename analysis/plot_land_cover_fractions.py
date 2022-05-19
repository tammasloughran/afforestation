#!/usr/bin/env python3
# plot_land_cover_fractions.py generates plots of the land cover fractions from the afforestation
# experiment and the esm-ssp585 experiment.
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

if __name__ != 'analysis.plot_afforestation':
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


def make_land_cover_plot():
    table = 'Lmon'
    years = list(range(2015, 2101))

    aff_areas = {}
    ssp585_areas = {}
    for var in FRAC_VARIABLES[table]:
        frac_file = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', table, var)
        aff_areas[var] = cdo_cover_area_load(frac_file[0], var)
        plt.plot(years, aff_areas[var]/M2_IN_MILKM2, color=COLORS[var], label=LABELS['AFF'][var])
        frac_file = get_filename('C4MIP', 'esm-ssp585', '1', table, var)
        ssp585_areas[var] = cdo_cover_area_load(frac_file[0], var)
        plt.plot(years, ssp585_areas[var]/M2_IN_MILKM2, color=COLORS[var], linestyle='dashed',
                label=LABELS['esm-ssp585'][var])
    plt.ylabel("Area (million km$^2$)")
    plt.xlabel("Year")
    plt.legend()


if __name__ != 'analysis.plot_land_cover_fractions':
    make_land_cover_plot()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()

