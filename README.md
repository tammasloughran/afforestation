# Afforestation

# Project structure
```
.
├── analysis
│   ├── baseline.py             Calculate baseline values over a base period.
│   ├── cdo_calc_load.py        CDO calculations and data loading.
│   ├── cmip_files.py           Locations and function for CMIP files.
│   ├── constants.py            Constants definitions.
│   ├── __init__.py             Analysis package init file. (empty)
│   ├── __main__.py             Main entry point program.
│   ├── plot_afforestation.py   Carbon plots script.
│   ├── plot_australia.py       Regional Plots (Australia)
│   └── plot_land_cover_fractions.py Plotting script for land cover types.
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

