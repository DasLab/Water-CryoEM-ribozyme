# Models

This folder contains the various atomic coordinates used in this study. The main atomic models can be found at [2.2A_SWIM_consensus.pdb](2.2A_SWIM_consensus.pdb) and [2.3A_SWIM_consensus.pdb](2.3A_SWIM_consensus.pdb), with the all automatically identified waters and ions found in [2.2A_SWIM.pdb](2.2A_SWIM.pdb) and [2.3A_SWIM.pdb](2.3A_SWIM.pdb). The 2.3 Å model was aligned to the 2.2 Å model, these coordinates can be found in [consensus_23_aligned.pdb](consensus_23_aligned.pdb). How the water and ion positions were generated can be found in [SWIM](SWIM).

During the study waters and ions from one map where placed in the other map, accounting for shifts between the two maps. These water and ion placements can be found in [22in23.pdb](22in23.pdb) and [23in22.pdb](23in22.pdb).

To assess the stringency of the SWIM criteria, random positions were sampled in the solvent shell around the 2.2 Å and 2.3 Å model: [random_solvent_22.pdb](random_solvent_22.pdb) and [random_solvent_23.pdb](random_solvent_23.pdb). See [Qscore](Qscore) for how this was generated, as well as how Q-score were generated throughout this study.

This study also compared against previously structure of the ribozyme or subdomains of the ribozyme. These models, split by chain when more than one chain was found in the asymmetric unit, can be found in [other_models](other_models). These models aligned to the structure in the current study can be found in [aligned_models](aligned_models).

