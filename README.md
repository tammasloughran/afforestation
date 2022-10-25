# Afforestation

# Project structure
```
.
├── analysis
│   ├── __init__.py             Analysis package init file.
│   ├── __main__.py             Main entry point program.
│   ├── cdo_calc_load.py        CDO calculations and data loading.
│   ├── cmip_files.py           Locations and function for CMIP files.
│   ├── constants.py            Constants definitions.
│   ├── cstock_stats.py         Table of C stock stats.
│   ├── jaisnb.py               A custom colormap.
│   ├── plot_afforestation.py   Carbon plots script.
│   ├── plot_atmosphere.py      Atmosphere plots script.
│   ├── plot_histograms.py      Histograms plots.
│   ├── plot_land_cover_fractions.py Plotting script for land cover types.
│   ├── plot_models.py          Plotting script for model intercomparison.
│   ├── plot_ocean.py           Plotting for ocean.
│   ├── plot_regions.py         Plotting for regional analysis.
│   └── plot_temp_lai.py        Temperature and LAI correlation.
├── data
│   └── README.md
├── doc                         Notes and paper draft go here.
│   ├── afforestation_draft.tex Notes on afforestation simulations.
│   └── afforestation_references.bib
├── LICENCE                     CSIRO Non-Commercial Source Code Licence Agreement v1.0
├── plots                       Location of output plots.
│   └── README.md
├── pyproject.toml              Project packaging backend configuration.
├── README.md
└── setup.cfg                   Installation configuration.
```

## Setup and run
The following installs the analysis scripts as a package in a virtual environment and runs the
analysis.
```bash
virtualenv affenv
source ./affenv/bin/activate
pip install git+https://gitlab.com/tammasloughran/afforestation.git
python -m analysis
```

## Generate ctags

```bash
ctags --recurse --python-kinds=-i analysis
```

