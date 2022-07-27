# Main entry point of ananlysis package
import os
import sys

from analysis.plot_afforestation import (
        make_clim_plots,
        make_veg_plots,
        )
from analysis.plot_land_cover_fractions import (
        make_land_cover_plot,
        make_afforestation_pft_plot,
        make_area_anomaly_map,
        )
from analysis.plot_regions import (
        make_regional_plots,
        plot_regions_map,
        )
from analysis.constants import PLOTS_DIR, DATA_DIR
from analysis.plot_atmoshpere import (
        make_co2_plot,
        )


files = os.listdir('.')
if 'data' not in files: od.mkdir(DATA_DIR)
if 'plots' not in files: od.mkdir(PLOTS_DIR)


def main():
    make_veg_plots()
    make_clim_plots()
    make_land_cover_plot()
    make_area_anomaly_map()
    make_afforestation_pft_plot()
    make_area_anomaly_map()
    plot_regions_map()
    make_regional_plots()
    make_co2_plot()


sys.exit(main())

