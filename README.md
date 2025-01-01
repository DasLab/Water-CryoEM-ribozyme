# Water-CryoEM-ribozyme

Repository accompanying the manuscript: 

## Maps and models

### SWIM
The code for the automated identification of water and ions used in this work to create the models used in this repository can be found at [https://github.com/gregdp/segger](https://github.com/gregdp/segger).

## Molecular dyanmics
Simulation input files, analysis, links for raw simulations, and data files used in analysis elsewhere in this respository can be found in [simulations](simulations).

## Analysis
All required data to recreate figures is included in this repository, in order to recreate these files, scripts in the [analysis](analysis).

## Figures
Code to create figures used in this paper can be found in [chimerax_code](chimerax_code) for figures of maps and models, and in [graphs_code](graphs_code) for graphical figures. See README in those respetive folders for more details.

#### Locations of graphing code:

- Fig 4A-B  [get_MD_density_values.ipynb](analysis/simulations/get_MD_density_values.ipynb)
- ED Fig 1H-I [Bfactor-plot.ipynb](analysis/bfactor/Bfactor-plot.ipynb)
- ED Fig 4I [analysis/water_consensus/swim_criteria_compare](analysis/water_consensus/swim_criteria_compare)
- ED Fig 9A,C,D,F [MD_summary_plots.ipynb](analysis/simulations/MD_summary_plots.ipynb)
- ED Fig 9E [get_MD_density_values.ipynb](analysis/simulations/get_MD_density_values.ipynb)

#### Revision response tables or figures
- Response #3.1 - running SWIM on half maps [analysis/water_consensus/swim_criteria_compare](analysis/water_consensus/swim_criteria_compare)