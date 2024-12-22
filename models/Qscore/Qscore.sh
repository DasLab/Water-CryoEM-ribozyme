#for pdb in Con2hf-4b_is2__Con2-2.2A_sh__thr2.000_Q0.70_0.60_rQ0.60__hA__hB__281-water__43-ion--.pdb; do for i in 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8 9; do echo ${i}; python $MAPQ_CMD $CHIMERA_LOC phase/2303_phase_simple_${i}.mrc ${pdb} np=24; done; done

# note outputs may end up in the maps folder

MAPQ_CMD=~/mapq/mapq_cmd.py
CHIMERA_LOC=~/.local/UCSF-Chimera64-1.14rc

PDB23=../2.3A_SWIM.pdb
MRC23=../../maps/Con1_2.3A_sh.mrc
MRC23A=../../maps/Con1_2.3A_sh_half_A.mrc
MRC23B=../../maps/Con1_2.3A_sh_half_B.mrc

PDB22=../2.2A_SWIM.pdb
MRC22=../../maps/Con2-2.2A_sh.mrc
MRC22A=../../maps/Con2-2.2A_sh_half_A.mrc
MRC22B=../../maps/Con2-2.2A_sh_half_B.mrc

# Need to generate these first run 
W23in22=../../models/23in22.pdb
W22in23=../../models/22in23.pdb

# the 2.3A map and half maps
python $MAPQ_CMD $CHIMERA_LOC map=$MRC23 $PDB23 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC23A $PDB23 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC23B $PDB23 np=24

# the 2.2A map and half maps
python $MAPQ_CMD $CHIMERA_LOC map=$MRC22 $PDB22 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC22A $PDB22 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC22B $PDB22 np=24

# the 2.2 A water and mg aligned to the 2.3 A model and  vice versa
python $MAPQ_CMD $CHIMERA_LOC map=$MRC22 $W23in22 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC22A $W23in22 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC22B $W23in22 np=24

python $MAPQ_CMD $CHIMERA_LOC map=$MRC23 $W22in23 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC23A $W22in23 np=24
python $MAPQ_CMD $CHIMERA_LOC map=$MRC23B $W22in23 np=24

# the previous 3.1 A map
python $MAPQ_CMD $CHIMERA_LOC map=./../maps/emd_31385.map ../other_models/7ez0.pdb np=24

