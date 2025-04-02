

Analyzing trajectories
RMSF+Bfactor CPPTRAJ
parm XXX_nowater.thr.prmtop
trajin XXX_md-nowater.nc 1 last # already every 1ns
autoimage

rms fit first (:1-387)&!@H=
average crdset MyAvg (:1-387)&!@H=
run
rms ref MyAvg (:1-387)&!@H=

atomicfluct rmsf out XXX_rmsf_all.dat :1-387&!@H= byres
atomicfluct bfactor out XXX_bfactor_all.dat :1-387&!@H= byres bfactor


import numpy as np
import glob
import pandas as pd

resnum = np.arange(22,409)

rmsf = []
for f in glob.glob("*rmsf_all*"):
    rmsf.append(np.loadtxt(f)[:,1])
rmsf = np.array(rmsf)
avg_rmsf = np.mean(rmsf,axis=0)

bfactor = []
for f in glob.glob("*bfactor_all*"):
    bfactor.append(np.loadtxt(f)[:,1])
bfactor = np.array(bfactor)
avg_bfactor = np.mean(bfactor,axis=0)

data = pd.DataFrame({"#Res":resnum,"RMSF":avg_rmsf,"Bfactor":avg_bfactor})
data.to_csv("cpptraj_avg_rmsf_qscore.csv",index=False)

Secondary and tertiary structure stability

Get PDBs
import numpy as np
import MDAnalysis

sims = ["n2ed1","n2ed2","n2ed3","n2ed4",
        "n2em1","n2em2","n2em3","n2em4",
        "n2en1","n2en2","n2en3","n2en4",
        "n2rd1","n2rd2","n2rd3","n2rd4","n2rd5","n2rd6",
        "n2rn1","n2rn2","n2rn3","n2rn4","n2rn5","n2rn6",
        "n2rm1","n2rm2","n2rm3","n2rm4","n2rm5","n2rm6"]

universes = []
total_frames = 0
for sim in sims:
    u = MDAnalysis.Universe(f"{sim}_nowater.thr.prmtop", f"{sim}_md-nowater.nc")
    total_frames += u.trajectory.n_frames
    universes.append(u)
RNA_selection = "not name H* and not name Na* and not name *M* and not name *Z* and not resname WAT and not (name O5' and resnum 22)"
for u,sim in zip(universes,sims):
    for ts in u.trajectory:
        rna = u.select_atoms(RNA_selection)
        rna.write(f"secstruct/{sim}_{ts.frame}ns.pdb")


Run RNA motif to find interactions
rna_motif.default.linuxgccrelease -in:file:s thr_MD.pdb
for f in secstruct/*pdb
do
  sed -i "s/ G3 / G  /g" $f
  sed -i "s/ G5 / G  /g" $f
  rna_motif.default.linuxgccrelease -in:file:s $f
done


Count interactions found in original and prevalence in simulation
import glob

bps = {}
with open('thr_MD.pdb.base_pairs.txt') as f:
    for line in f:
        bps[line] = 0


for file in glob.glob("secstruct_*/*base_pairs.txt"):
    with open(file) as f:
        for line_ in f:
            line = line_.replace("X:SYST:"," :")
            if line in bps:
                bps[line] += 1


print(bps)

only_WC = {}
for bp,count in bps.items():
    if "W W C" in bp:
        only_WC[bp] = count

print(only_WC)


RMSD plots
See Fig S9A section of notebook.
Simulated Q-score
Append all pdbs (made in secondary and tertiary structure section) of each simulation to make a large PDB with all 400 1ns frames
for f in secstruct_n2en/n2en4*.pdb; do tail -n +3 $f | head -n 8279 >> n2en4_all.pdb; done


Make map in ChimeraX
molmap #1 2.2


Use average structure calculated by RMSF section in and map and calculate Q-score, average structure and map can be found here
for f in n2rd*avg_struct.pdb;
do
  echo $f ${f:0:5}_all_2.2.mrc;
  python ~/mapq/mapq_cmd.py ~/.local/UCSF-Chimera64-1.14rc ${f:0:5}_all_2.2.mrc $f;
done


Correlation Qscore, RMSF, Simulated Q-score, Bfactor
Plot per residue Qscore, RMSF, simulated Q
See Fig 3G section of notebook
Plot correlation plots
See Fig S8 section of notebook.
Plot Qscore nucleotide, base, backbone
See Fig 1F section of notebook.
Compare Mg with previous structures
See table S4 and Fig S6A sections of the notebook.

Binding locations
Compare binding site in MD and EM
See Fig S6B,C in notebook.
Summary ion and water distance to RNA in MD
See Fig 6D section in notebook.


MD Movies - pymol
load n2XXX_nowater.thr.prmtop, XXX
load_traj n2XXX_md-nowater.nc, XXX
bg_color white
set cartoon_ladder_mode, 1
set cartoon_ring_finder, 1
set cartoon_ring_mode, 3
alter polymer.nucleic, resi=str(int(resi)+21)
color grey, all
select P9_2, resi 369-402
color 0xde9ccc, P9_2
select P9_1, resi 353-367+332-346
color 0x99fc99, P9_1
select P9, resi 318-331
color 0xe82421, P9
select P7, resi 307-314+262-268
color 0x2bbfbf, P7
select P3, resi 272-278+96-102 
color 0xf2ba45, P3
select P8, resi 279-299
color 0xab21ba, P8
select P13, resi 73+75-81+347-352
color 0xfa5424, P13
select P2, resi 31-42+46-56
color 0x172696, P2
select P2_1, resi 58-68+82-93
color 0x963636, P2_1
select P14, resi 43-45+169-172 
color 0xed40f5, P14
select P6, resi 215-258
color 0x0d800d, P6
select P4, resi 208-214+107-112
color 0xfa9687, P4
select P5, resi 116-121+200-205 
color 0xed8226, P5
select P9a, resi 315-316+405-406 
color 0xf2e8a1, P9a
select P5abc, resi 127-168+173-195 
color 0x002bf7, P5abc
select Mg, name MG or name nMg or name mMg
color magnesium, Mg
select Na, not polymer.nucleic and not (name MG or name nMg or name mMg)
color sodium, Na
alter Na, vdw=0.8
show spheres, Na

smooth
reset
morph XXX_morph, XXX, XXX, 400, 1, steps=20

join_states XXX_all, XXX*, 0
disable XXX_morph

reset
mview reset
mset 1 x20 1 -400 400 x60 400 -420 1 x20 1 -400 400 x60 400 -420 1 x20 1 -400 400 x60  400 -420 1 x20 1 -400 400 x60  400 -420 1 x20 1 -400 400 x60  400 -420 1 x20 1 -400 400 x60
mview store, 1
mview store, 420
turn y, 90
mview store, 480
mview store, 920
turn x, 90
mview store, 980
mview store, 1420
turn x, 90
mview store, 1480
mview store, 1920
turn x, 90
mview store, 1980
mview store, 2420
turn x, 90
mview store, 2480
mview store, 2920
turn y, 90
mview store, 2980

movie.produce XXX.mpg, quality=90


