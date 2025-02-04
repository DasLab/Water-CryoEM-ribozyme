# Complex Water Networks Visualized by Cryogenic Electron Microscopy of RNA

Repository accompanying the manuscript: 

[preprint]

Rachael C. Kretsch, Shanshan Li, Grigore Pintilie, Michael Z. Palo, David A. Case, Rhiju Das, Kaiming Zhang, Wah Chiu. Complex Water Networks Visualized through 2.2-2.3 Å Cryogenic Electron Microscopy of RNA. bioRxiv 2025.01.23.634578; doi: https://doi.org/10.1101/2025.01.23.634578

## Maps and models

The atomic coordinates used in this study can be found in [models](models). The code for the automated identification of water and ions used in this work to create the models used in this repository can be found at, [SWIM](models/SWIM) they were based on [https://github.com/gregdp/segger](https://github.com/gregdp/segger).

## Molecular dynamics
Simulation input files, analysis, links for raw simulations, and data files used in analysis elsewhere in this repository can be found in [simulations](simulations).

## Analysis
All required data to recreate figures is included in this repository, in order to recreate these files, scripts in the [analysis](analysis). Code to create figures of atomic models are maps can be found in [chimerax_code](chimerax_code). See README in those respective folders for more details.

#### Locations of graphing code:

- Fig 2A-B [models_water_mg_by_consensus.cxc](chimerax_code/models_water_mg_by_consensus.cxc)
- Fig 2C-E [consensus_nonconsensus_compare.ipynb](analysis/water_consensus/consensus_nonconsensus_compare.ipynb)
- Fig 4A-B  [get_MD_density_values.ipynb](analysis/simulations/get_MD_density_values.ipynb)
- Fig 5 [Visualize_water_networks](chimerax_code/visualize_water_networks)
- ED Fig 1H-I [Bfactor-plot.ipynb](analysis/bfactor/Bfactor-plot.ipynb)
- ED Fig 2C [Qscore_22_23_31.ipynb](analysis/per_residue_comparison/Qscore_22_23_31.ipynb)
- ED Fig 3A-D [localres_rmsd_rmsf_q.cxc](chimerax_code/localres_rmsd_rmsf_q.cxc)
- ED Fig 3E [nucleotide_v_qscore_rmsf.ipynb](analysis/per_residue_comparison/nucleotide_v_qscore_rmsf.ipynb)
- ED Fig 3F-H [per_residue_correlations.ipynb](analysis/per_residue_comparison/per_residue_correlations.ipynb)
- ED Fig 4I [swim_criteria_compare](analysis/water_consensus/swim_criteria_compare)
- ED Fig 4L,M [consensus_nonconsensus_compare.ipynb](analysis/water_consensus/consensus_nonconsensus_compare.ipynb)
- ED Fig 4J,K,N,O [SWIM_summary_Qscore_distances.ipynb](analysis/binding_motifs/SWIM_summary_Qscore_distances.ipynb)
- ED Fig 5A-E [models_water_mg_by_agreement.cxc](chimerax_code/models_water_mg_by_agreement.cxc)
- ED Fig 5F-G [water_mg_consensus_tables.ipynb](analysis/water_consensus/water_mg_consensus_tables.ipynb)
- ED Fig 7B-D [interaction_types.ipynb](analysis/binding_motifs/interaction_types.ipynb)
- ED Fig 7E [avg_around_nuc](chimerax_code/avg_around_nuc)
- ED Fig 8A,B [categorize_water_ions.ipynb](analysis/water_consensus/categorize_water_ions.ipynb)
- ED Fig 9A,C,D,F [MD_summary_plots.ipynb](analysis/simulations/MD_summary_plots.ipynb)
- ED Fig 9E [get_MD_density_values.ipynb](analysis/simulations/get_MD_density_values.ipynb)
- Supp Table 1 [categorize_water_ions.ipynb](analysis/water_consensus/categorize_water_ions.ipynb)

#### Revision response tables or figures
- Response #3.1 - running SWIM on half maps [swim_criteria_compare](analysis/water_consensus/swim_criteria_compare)
- Response #3.5 - Q-score waters based on binding base vs backbone [consensus_nonconsensus_compare.ipynb](analysis/water_consensus/consensus_nonconsensus_compare.ipynb)