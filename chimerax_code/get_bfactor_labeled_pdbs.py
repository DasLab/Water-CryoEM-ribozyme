# import utils
import sys, os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../analysis'))
sys.path.insert(0,parent_dir)
from utils import *


import MDAnalysis
from MDAnalysis.analysis import distances


# inputs
per_res_file = "../analysis/per_residue_comparison/per_residue_summary.csv"
consensus_criteria = "both" #either
consensus_status_file = "../analysis/water_consensus/all_solvent_consensus_status_withconsensus.csv"

# outputs
rmsf_22_out = "bfactor_models/22_rmsf.pdb"
rmsd_22_out = "bfactor_models/22_rmsd.pdb"
q_22_out = "bfactor_models/22_q.pdb"
mg_wat_cons_22 = "bfactor_models/22_mg_water_consensus.pdb"
mg_wat_cons_23_align = "bfactor_models/23_aligned_mg_water_consensus.pdb"
mg_wat_cons_23_unalign = "bfactor_models/23_unaligned_mg_water_consensus.pdb"
consensus_3A_23 = "bfactor_models/23_colored_3A_agreement.pdb"
consensus_3A_22_q = "bfactor_models/22_colored_3A_agreement_Q.pdb"
consensus_3A_22_rmsd = "bfactor_models/22_colored_3A_agreement_RMSD_avg.pdb"
consensus_23_22_q = "bfactor_models/22_colored_2.3A_agreement_Q.pdb"
consensus_23_22_rmsd = "bfactor_models/22_colored_2.3A_agreement_RMSD.pdb"
consensus_X_23 = 'bfactor_models/23_colored_xray_agreement.pdb'
consensus_X_22_q = "bfactor_models/22_colored_xray_agreement_Q.pdb"
consensus_X_22_rmsd = "bfactor_models/22_colored_xray_agreement_num_structs.pdb"


# save pdbs for visualization of RMSD, RMSF, and Q-score
df = pd.read_csv(per_res_file)
save_pdb_with_bfactor_by_res(PDB22_F,RNA_SELECTION,df.RMSF,rmsf_22_out)
save_pdb_with_bfactor_by_res(PDB22_F,RNA_SELECTION,df.RMSD,rmsd_22_out)
save_pdb_with_bfactor_by_res(PDB22_F,RNA_SELECTION,df['Qscore 2.2A'],q_22_out)


# write to pdbs for visualization of consensus status
df = pd.read_csv(consensus_status_file)
matching_22_water,matching_22_mg,matching_23_water,matching_23_mg = get_consensus_wat_mg(df,consensus_criteria)
u22 = MDAnalysis.Universe(PDB22_F)
u23 = MDAnalysis.Universe(PDB23_ALIGNED)
    
bfactors = np.zeros(len(u22.atoms))
for resnum in matching_22_water:
    index = u22.select_atoms(f"{WAT_SELECTION} and resnum {resnum}").indices
    bfactors[index] = 1
for resnum in matching_22_mg:
    index = u22.select_atoms(f"{MG_SELECTION} and resnum {resnum}").indices
    bfactors[index] = 1
u22.atoms.tempfactors = bfactors
u22.atoms.write(mg_wat_cons_22)

bfactors = np.zeros(len(u23.atoms))
for resnum in matching_23_water:
    index = u23.select_atoms(f"{WAT_SELECTION} and resnum {resnum}").indices
    bfactors[index] = 1
for resnum in matching_23_mg:
    index = u23.select_atoms(f"{MG_SELECTION} and resnum {resnum}").indices
    bfactors[index] = 1
u23.atoms.tempfactors = bfactors
u23.atoms.write(mg_wat_cons_23_align)
u23_unaligned = MDAnalysis.Universe(PDB23_F)
u23_unaligned.atoms.tempfactors = bfactors
u23_unaligned.atoms.write(mg_wat_cons_23_unalign)

# print relevant counts
print(len(matching_22_mg))
print(len(df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG")])-len(matching_22_mg))
print(len(matching_22_water))
print(len(df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH")])-len(matching_22_water))
print(len(matching_23_mg))
print(len(df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG")])-len(matching_23_mg))
print(len(matching_23_water))
print(len(df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH")])-len(matching_23_water))


# save detailed on consensus with ion and water
# Now 5 colors
# 0 - no overlap
# 0.25 - bind at least 1 of same atom
# 0.5 - within 1A after alignment of RNA core
# 0.75 - bound to same RNA atoms
# 1.0 - consensus
per_res_info = pd.read_csv(per_res_file)
df = pd.read_csv(consensus_status_file)

to_draw = {"2.3":["2.3Å"],
   "3A":['7ez0', '7ez2', '7r6l',
       '7xd5', '7xd6', '7xd7', '7yg8', '7yg9',
       '7yga', '7ygb', '7ygc', '8i7n'],
             "xray":['1gid_A','1gid_B','1hr2_A','1hr2_B','1x8w_A','1x8w_B','1x8w_C','1x8w_D']}

our22df = df[df.model=="2.2Å"].copy()
our22df['2.3 bfactor'] = our22df.apply(lambda row: get_bfactor_based_on_overlap(row,to_draw["2.3"],"mg") if row.solvent=="MG" else get_bfactor_based_on_overlap(row,to_draw["2.3"],"wat"),axis=1)
our22df['3A bfactor'] = our22df.apply(lambda row: get_bfactor_based_on_overlap(row,to_draw["3A"],"mg") if row.solvent=="MG" else get_bfactor_based_on_overlap(row,to_draw["3A"],"wat"),axis=1)
our22df['xray bfactor'] = our22df.apply(lambda row: get_bfactor_based_on_overlap(row,to_draw["xray"],"mg") if row.solvent=="MG" else get_bfactor_based_on_overlap(row,to_draw["xray"],"wat"),axis=1)

our23df = df[df.model.isin(["2.2Å","2.3Å"])].copy()
our23df['3A bfactor'] = our23df.apply(lambda row: get_bfactor_based_on_overlap(row,to_draw["3A"],"mg") if row.solvent=="MG" else get_bfactor_based_on_overlap(row,to_draw["3A"],"wat"),axis=1)
our23df['xray bfactor'] = our23df.apply(lambda row: get_bfactor_based_on_overlap(row,to_draw["xray"],"mg") if row.solvent=="MG" else get_bfactor_based_on_overlap(row,to_draw["xray"],"wat"),axis=1)

u22 = PandasPdb().read_pdb(PDB22_F).df["ATOM"]
u23 = PandasPdb().read_pdb(PDB23_ALIGNED).df["ATOM"]

# get RMSD
cryoem_3_struct = [f'../models/aligned_models/{pdb}.pdb' for pdb in to_draw['3A']]
#cryoem_3_struct = ['../models/aligned_models/7ez0_aligned.pdb','../models/aligned_models/7r6n_aligned.pdb', 
#                   '../models/aligned_models/7r6m_aligned.pdb', '../models/aligned_models/7r6l_aligned.pdb',
#       '../models/aligned_models/7ez2_aligned.pdb']

rmsds22 = []
rmsds23 = []
for struct in cryoem_3_struct:
    u = PandasPdb().read_pdb(struct).df["ATOM"]
    u = u.merge(u22,on=["atom_name","residue_name","residue_number"]).reset_index()
    print(struct,"num atoms",len(u),"num_res",len(u.residue_number.unique()))
    diff = u[["x_coord_x","y_coord_x","z_coord_x"]].to_numpy()-u[["x_coord_y","y_coord_y","z_coord_y"]].to_numpy()#.astype(float)-coords.astype(float)
    dist = np.linalg.norm(diff,axis=1)
    u["dev"] = dist*dist
    rmsd = u.groupby("residue_number").dev.mean()
    rmsd = np.sqrt(rmsd)
    rmsd = rmsd.reset_index()
    rmsds22.append(rmsd)
    u = PandasPdb().read_pdb(struct).df["ATOM"]
    u = u.merge(u23,on=["atom_name","residue_name","residue_number"]).reset_index()
    print(struct,"num atoms",len(u),"num_res",len(u.residue_number.unique()))
    diff = u[["x_coord_x","y_coord_x","z_coord_x"]].to_numpy()-u[["x_coord_y","y_coord_y","z_coord_y"]].to_numpy()#.astype(float)-coords.astype(float)
    dist = np.linalg.norm(diff,axis=1)
    u["dev"] = dist*dist
    rmsd = u.groupby("residue_number").dev.mean()
    rmsd = np.sqrt(rmsd)
    rmsd = rmsd.reset_index()
    rmsds23.append(rmsd)
rmsds22 = pd.concat(rmsds22)
rmsds23 = pd.concat(rmsds23)

cryoem_avg_rmsd22 = rmsds22.groupby("residue_number").mean().reset_index()
cryoem_avg_rmsd23 = rmsds23.groupby("residue_number").mean().reset_index()
cryoem_avg_rmsd = (cryoem_avg_rmsd22+cryoem_avg_rmsd23)/2
print(cryoem_avg_rmsd)

# get xray count
xray_struct = [f'../models/aligned_models/{pdb}.pdb' for pdb in to_draw['xray']]

#xray_struct = [MODEL_RENAME_R[pdb] for pdb in to_draw['xray']]

#xray_struct = ['../models/aligned_models/1x8w_D_aligned.pdb',
# '../models/aligned_models/1gid_A_aligned.pdb',
# '../models/aligned_models/1x8w_C_aligned.pdb',
# '../models/aligned_models/1hr2_B_aligned.pdb',
# '../models/aligned_models/1x8w_B_aligned.pdb',
# '../models/aligned_models/1hr2_A_aligned.pdb',
# '../models/aligned_models/1x8w_A_aligned.pdb',
# '../models/aligned_models/1gid_B_aligned.pdb']

u22 = MDAnalysis.Universe(PDB22_F)
rna_atoms = u22.select_atoms(RNA_SELECTION)
residues_xray_count = {res:0 for res in np.unique(rna_atoms.resids)}
residues_xray_count_wat = {res:0 for res in np.unique(rna_atoms.resids)}
for struct in xray_struct:
    u = MDAnalysis.Universe(struct)
    for res in residues_xray_count.keys():
        closest_dist = distances.distance_array(u22.select_atoms(f"{RNA_SELECTION} and resid {res}").positions, # reference
                                    u.select_atoms(RNA_SELECTION).positions, # configuration
                                    box=u22.dimensions).min()
        if closest_dist < 1:
            residues_xray_count[res] = residues_xray_count[res] + 1
            if '1hr2' in struct or '1gid' in struct:
                residues_xray_count_wat[res] = residues_xray_count_wat[res] + 1

u22 = MDAnalysis.Universe(PDB22_F)
u23 = MDAnalysis.Universe(PDB23_ALIGNED)

bfactors = np.zeros(len(u22.atoms))
bfactors23 = np.zeros(len(u23.atoms))
for i,row in our22df.iterrows():
    i,val = change_value_bfactor(row,u22,'3A bfactor')
    bfactors[i] = val
for i,row in our23df.iterrows():
    i,val = change_value_bfactor(row,u23,'3A bfactor')
    bfactors23[i] = val
u23.atoms.tempfactors = bfactors23
u23.atoms.write(consensus_3A_23)
for i,row in per_res_info.iterrows():
    i,val = change_value_bfactor(row,u22,'Qscore 2.2A',"residue_number")
    bfactors[i] = val
u22.atoms.tempfactors = bfactors
u22.atoms.write(consensus_3A_22_q)
for i,row in cryoem_avg_rmsd.iterrows():
    i,val = change_value_bfactor(row,u22,'dev',"residue_number")
    bfactors[i] = val
u22.atoms.tempfactors = bfactors
u22.atoms.write(consensus_3A_22_rmsd)

bfactors = np.zeros(len(u22.atoms))
for i,row in our22df.iterrows():
    i,val = change_value_bfactor(row,u22,'2.3 bfactor')
    bfactors[i] = val
for i,row in per_res_info.iterrows():
    i,val = change_value_bfactor(row,u22,'Qscore 2.2A',"residue_number")
    bfactors[i] = val
u22.atoms.tempfactors = bfactors
u22.atoms.write(consensus_23_22_q)
for i,row in per_res_info.iterrows():
    i,val = change_value_bfactor(row,u22,'RMSD',"residue_number")
    bfactors[i] = val
u22.atoms.tempfactors = bfactors
u22.atoms.write(consensus_23_22_rmsd)

bfactors = np.zeros(len(u22.atoms))
for i,row in our22df.iterrows():
    i,val = change_value_bfactor(row,u22,'xray bfactor')
    bfactors[i] = val
for i,row in our23df.iterrows():
    i,val = change_value_bfactor(row,u23,'xray bfactor')
    bfactors23[i] = val
u23.atoms.tempfactors = bfactors23
u23.atoms.write(consensus_X_23)
for i,row in per_res_info.iterrows():
    i,val = change_value_bfactor(row,u22,'Qscore 2.2A',"residue_number")
    bfactors[i] = val
u22.atoms.tempfactors = bfactors
u22.atoms.write(consensus_X_22_q)
for i,row in per_res_info.iterrows():
    i,val = change_value_bfactor(row,u22,'Qscore 2.2A',"residue_number")
    bfactors[i] = residues_xray_count[row.residue_number]
u22.atoms.tempfactors = bfactors
u22.atoms.write(consensus_X_22_rmsd)


# print out region to hide for xray (not present in any xray structure)
print(",".join([str(x) for x in residues_xray_count if residues_xray_count[x]==0]))
print(",".join([str(x) for x in residues_xray_count_wat if residues_xray_count_wat[x]==0]))

'''
double check same
22,23,24,25,26,27,28,29,30,31,36,37,38,39,40,41,42,43,45,46,47,48,49,50,51,52,57,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,287,288,289,290,291,332,333,334,335,336,337,338,339,340,341,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401
22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,143,144,260,262,263,264,266,267,268,270,271,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408

'''

