# Calculate the Q-score of models

## Qscore

Note to enable no loss of significant figures in Q-score, a function was added to the Q-score scripts, [added_function_to_qscores.txt](added_function_to_qscores.txt) to save this information. This produces a csv files which is then easily readable in future scripts.

### Q-score of maps

To calculate the Q-score of water, magnesium, and RNA atoms, of our and previous models we used the simply command `python $MAPQ_CMD $CHIMERA_LOC map=$MRC $PDB np=24`, full commands used can be found in [Qscore.sh](Qscore.sh).

### Q-score of random solvent positions
To test to strigency of SWIM, we wanted to know what the Qscore of sampling of positions in the RNA solvent shell was (1.5-3.5A for any RNA heavy-atom). First get a pdb of potential solvent locations for both the 2.2A and 2.3A structure using [get_random_solvent.py](get_random_solvent.py). Note, the pdb included in this repository were generated for earlier models, results do not differ significantly, but pdbs prodcued with deposited structures may differ.

```
python get_random_solvent.py 0.6 PDB23_F ../models/random_solvent_23.pdb
python get_random_solvent.py 0.6 PDB22_F ../models/random_solvent_22.pdb
```

For reasonable processing times, we split these pdbs up using [split_rand_pdb_for_processing.py](split_rand_pdb_for_processing.py)

```
python split_rand_pdb_for_processing.py ../models/random_solvent_22.pdb ../models/random_solvent_22_part
python split_rand_pdb_for_processing.py ../models/random_solvent_23.pdb ../models/random_solvent_23_part
```

Then for each calculated the Q-scores, for full and half maps:
```
python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc map=../../maps/Con1_2.3A_sh.mrc ../models/random_solvent_23_part${i}.pdb np=16
python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc ../../maps/Con2-2.2A_sh.mrc ../models/random_solvent_22_part${i}.pdb np=16

python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc map=../../maps/Con1_2.3A_sh_half_A.mrc ../models/random_solvent_23_part${i}.pdb np=16
python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc map=../../maps/Con1_2.3A_sh_half_B.mrc ../models/random_solvent_23_part${i}.pdb np=16
python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc map=../../maps/Con2-2.2A_sh_half_A.mrc ../models/random_solvent_22_part${i}.pdb np=16
python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc map=../../maps/Con2-2.2A_sh_half_B.mrc ../models/random_solvent_22_part${i}.pdb np=16
```

Then all Q-scores were combined:

```
python combine_rand_pdb_Q.py ../models/random_solvent_23_part 2.3 ../models/random_solvent_22_part  2.2 random_sovlent_shell_Q.csv random_sovlent_shell_halfQ.csv
```

Creating the final [random_sovlent_shell_Q.csv](outputs/random_sovlent_shell_Q.csv) and [random_sovlent_shell_halfQ.csv](outputs/random_sovlent_shell_halfQ.csv).
