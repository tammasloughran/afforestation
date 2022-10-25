#!/usr/bin/env python3
"""Create some stats of the carbon stocks.
"""
import numpy as np
from markdownTable import markdownTable
if __name__ == 'analysis.plot_regions':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.constants import (
            DATA_DIR,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from constants import (
            DATA_DIR,
            )

VARIABLES = (
        'cLand',
        'cVeg',
        'cLitter',
        'cSoil',
        )


def print_stats()->None:
    # Load data
    ssp585_data = {}
    for_data = {}
    diff_data = {}
    last_20_mean = {}
    diff_p_present = {}
    diff_p_ssp585 = {}
    for var in VARIABLES:
        ssp585_data[var] = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_global.npy')
        for_data[var] = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
        diff_data[var] = for_data[var] - ssp585_data[var]

        # Create stats
        # Mean last 20 years of difference
        last_20_mean[var] = diff_data[var][:,-20:].mean(axis=1)

        # The last 20 years difference as a percentage compared to the first time step (present day)
        diff_p_present[var] = (100*diff_data[var]/ssp585_data[var][:,0,None])[:,-20:].mean(axis=1)

        # The last 20 years difference as a percentage compared to the last 20 years of ssp585.
        diff_p_ssp585[var] = (
                100*diff_data[var][:,-20:].mean(axis=1) / ssp585_data[var][:,-20:].mean(axis=1)
                )

    # Prinnt as a table
    table = [ {
            '':var,
            '$\Delta$ Carbon (Pg)':f'{last_20_mean[var].mean():.1f}±{last_20_mean[var].std():.1f}',
            '% present day':f'{diff_p_present[var].mean():.1f}±{diff_p_present[var].std():.1f}',
            '% of 2080-2100 ssp585':f'{diff_p_ssp585[var].mean():.1f}±{diff_p_ssp585[var].std():.1f}',
            } for var in VARIABLES]
    print(markdownTable(table).setParams(row_sep='markdown').getMarkdown())


if __name__ != 'analysis.cstock_stats':
    print_stats()

