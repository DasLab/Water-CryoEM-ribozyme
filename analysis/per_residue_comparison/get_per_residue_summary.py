
# import utils
import sys, os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0,parent_dir)
from utils import *

import MDAnalysis
from MDAnalysis.analysis.rms import RMSD

# inputs
outfile = "per_residue_summary.csv"


# Read rmsf values from running calculate_rmsf.py
rmsf_data = []
for sim in SIMULATIONS:
    pdb = PandasPdb().read_pdb(f'../../simulations/simulation_rmsfs/{sim}_rmsf.pdb')
    df = pd.concat([pdb.df["ATOM"],pdb.df["HETATM"]])
    df["simulation"] = sim
    df = df.rename(columns={"b_factor":"rmsf"})
    df = df[df.element_symbol!="H"]
    df["residue_number"] = df["residue_number"]+21
    df = df[~((df.atom_name=="O5'") &(df["residue_number"]==22))] # not in original model
    rmsf_data.append(df)
rmsf_data = pd.concat(rmsf_data)

# get per residue rmsf
per_res_rmsf = rmsf_data[~rmsf_data.residue_name.isin(["HOH","MG","Na+","mMg","nMg"])].groupby("residue_number").rmsf.mean().to_numpy()

# Align all RNA heavy atoms and calculate rmsd
u22 = MDAnalysis.Universe(PDB22_F)
u23 = MDAnalysis.Universe(PDB23_F)

R22_23 = RMSD(u22.select_atoms(RNA_SELECTION), u23.select_atoms(RNA_SELECTION),
              superposition=RNA_SELECTION,groupselections=[f"{RNA_SELECTION} and resnum {x}" for x in np.unique(u22.select_atoms(RNA_SELECTION).resnums)])
R22_23.run()
rmsd22_23 = R22_23.results.rmsd.T 
total_rmsd = rmsd22_23[2]
per_res_rmsd_22_23 = rmsd22_23[3:]

# read and format Qscores
col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','Qscore']

qscore_22 = pd.read_csv(QSCORE_22,names=col_names)
qscore_22["residue_number"] = qscore_22.residue_number.apply(lambda x: int(x.split(".")[0]))
qscore_23 = pd.read_csv(QSCORE_23,names=col_names)
qscore_23["residue_number"] = qscore_23.residue_number.apply(lambda x: int(x.split(".")[0]))

# ignore residue 414
qscore_22 = qscore_22[qscore_22.residue_number!=414]
qscore_23 = qscore_23[qscore_23.residue_number!=414]

# average over all RNA reidues
q_score_per_residue_22 = qscore_22[(~qscore_22.residue_name.isin(['MG', 'HOH']))&(qscore_22.atom_name.str[0]!="H")].groupby("residue_number").Qscore.mean()
q_score_per_residue_23 = qscore_23[(~qscore_23.residue_name.isin(['MG', 'HOH']))&(qscore_23.atom_name.str[0]!="H")].groupby("residue_number").Qscore.mean()

col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','Qscore']
simdfs=[]
for sim in SIMULATIONS:
    simdf = pd.read_csv(f'../../simulations/simulation_Q/{sim}_avg_struct.pdb__Q__{sim}_all_2.2.mrc.csv',names=col_names)
    simdf["residue_number"] = simdf.residue_number.apply(lambda x: int(x.split(".")[0]))
    simdfs.append(simdf)
simdfs = pd.concat(simdfs)

# RMSD, RMSF, Qscore 2.2, Qscore 2.3, resmap value, simulated Qscore
resmap23 = pd.read_csv(LOCALRES_23)
resmap22 = pd.read_csv(LOCALRES_22)
resmap22 = resmap22[resmap22.residue_name.isin(["A","C","G","U"])]
resmap23 = resmap23[resmap23.residue_name.isin(["A","C","G","U"])]
# ignore residue 414
resmap22 = resmap22[resmap22.residue_number!=414]
resmap23 = resmap23[resmap23.residue_number!=414]

# combine all 3 per residue values
df = q_score_per_residue_22.reset_index().rename(columns={"Qscore":"Qscore 2.2A"})
df["Qscore 2.3A"] = q_score_per_residue_23.reset_index().Qscore
df["local resolution 2.2A"] = resmap22.groupby("residue_number").mean().mapvalue.values
df["local resolution 2.3A"] = resmap23.groupby("residue_number").mean().mapvalue.values
df["RMSD"] = per_res_rmsd_22_23
df["RMSF"] = per_res_rmsf
df["simQ"] = simdfs.groupby("residue_number").mean().Qscore.values
# save 
df.to_csv(outfile,index=False)
