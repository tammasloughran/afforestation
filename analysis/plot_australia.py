#!/usr/bin/env python3
import matplotlib.pyplot as plt

if __name__ != 'analysis.plot_australia':
    # plot_afforestation.py is main program or imported as a module from another script.
    from constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES
    from cmip_files import get_filename
    from cdo_calc_load import load_aus_pool
else:
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES
    from analysis.cmip_files import get_filename
    from analysis.cdo_calc_load import load_aus_pool

for ens in ENSEMBLES:
    filename = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Lmon', 'cVeg')[0]
    data = load_aus_pool(input=filename, var='cVeg')
    plt.plot(data)

plt.show()

