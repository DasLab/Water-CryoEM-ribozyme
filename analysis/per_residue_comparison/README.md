## Extract per-residue information

Run `python get_per_residue_summary.py` to create [per_residue_summary.csv](per_residue_summary.csv) which contains the Q-score and local resolution of each residue in each map, the RMSD of that residue after all-heavy atom superposition between the two structures, RMSF of the residue in simulation after all-heavy atom superposition, and simQ, which is the average Q-score of the residue between simulated map and average structure from the MD simulations.
