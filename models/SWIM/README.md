# Applying SWIM to the Cryo-EM maps

## SWIM installation

```
git clone https://github.com/gregdp/segger.git
cd segger/download
unzip Segger_2_9_7.zip
cd Segger
cp ${GIT_REPO}/models/SWIM/SWIM.py . # an updated version of SWIM.py is needed, for more flexibility with parameters, no substantial changes to algorithmn were made.
python install.py $CHIMERA_LOC
# found chimera folder with: chimera --root
```
## Running SWIM

A script for running SWIM from the command line was modified for this paper, [SWIM_script.py](SWIM_script.py). Please note it required the updated [SWIM.py](SWIM.py) found here and mentioned in install instructions above.

The commands used to run SWIM for the main models in the manuscript can be found in [run_swim.sh](run_swim.sh). The scripts for running SWIM while varying the desnity and Qscore thresholds can be found in [run_swim_metrics.sh](run_swim_metrics.sh) and the scirpts for running SWIM without half-maps and on half-maps in [run_swim_half.sh](run_swim_half.sh).

All outputs, models and outputs, can be found in [models](models). The models used throughout the manuscript were taken from [2.2A-G414_cdist3.2_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb](models/2.2A-G414_cdist3.2_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb) and [2.3A-G414_cdist3.2_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb](models/2.3A-G414_cdist3.2_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb) for the 2.2 and 2.3 Å maps respectively. The residue G414 which was added for the SWIM run was deleted before further analysis, [2.2A_SWIM.pdb](2.2A_SWIM.pdb) and [2.3A_SWIM.pdb](2.3A_SWIM.pdb).

