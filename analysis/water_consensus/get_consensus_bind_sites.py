# import utils
import sys, os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0,parent_dir)
from utils import *

import pandas as pd

dfile = 'all_solvent_consensus_status.csv'
df = pd.read_csv(dfile)

def get_if_subset(siteA,siteB):
    # siteB subset siteA
    if siteA != siteA or siteA=="":
        return False
    siteA_list = siteA.split()
    for atom in siteB.split():
        if atom not in siteA_list:
            return False
    return True

def get_sym_overlap(cA,fB,cB,fA):
    if fB != fB or fB == "":
        return False
    if cB != cB or cB == "":
        cB = fB
    return get_if_subset(fB,cA) and get_if_subset(fA,cB)

def get_number_sym_overlap(cA,fBs,cBs,fA):
    if fA != fA or fA == "":
        return 0
    if cA != cA or cA == "":
        cA = fA
    overlap = [get_sym_overlap(cA,fB,cB,fA) for cB,fB in zip(cBs,fBs)]
    return sum(overlap)

def get_consensus(row,model,df,solv):
    if row[f'within 1A of {solv} in {model}'] == 0:
        return False
    cA = row["close_binding_site"]
    fA = row["binding_site"]
    if fA != fA or fA == "":
        return False
    if cA != cA or cA == "":
        cA = fA
    resnum = row[f'within 1A of {solv} in {model}'].split(":")[0].split("-")[1]
    compare = df[(df.residue_number == int(resnum)) & (df.model == model)]
    cB = compare["close_binding_site"].item()
    fB = compare["binding_site"].item()
    # this check is redudant
    #if compare[f'within 1A of {solv} in {row.model}'].item() == 0:
    #	return False
    if fB != fB or fB == "":
        return False
    if cB != cB or cB == "":
        cB = fB
    return get_if_subset(fB,cA) and get_if_subset(fA,cB)

from biopandas.pdb import PandasPdb

# remark the within 1A with what atom it is close to (imperfect)
new_df = []
all_columns = [col for col in df.columns if "within 1A" in col]
# load all PDBs into dict
pdbs = {}
for model in df['model'].unique():
    if model == '2.2Å':
        pdb_f = PDB22_F #'../models/Con2hf-4b_is2__Con2-2.2A_sh__thr2.000_Q0.70_0.60_rQ0.60__hA__hB__297-water__54-ion--.pdb'
    elif model == '2.3Å':
        pdb_f = PDB23_ALIGNED #f'../models/aligned_models/consensus_23.pdb'
    elif model in ['7r6m','7r6n','7ygd']: #,'8hd6','8hd7' '7xd3','7xd4',
        pdb_f = f'../../models/other_models/{model}.pdb'
    else:
        pdb_f = f'../../models/aligned_models/{model}.pdb'#f'../models/aligned_models/{model}.pdb'
    pdb = PandasPdb().read_pdb(pdb_f)
    print(pdb.df['HETATM'].residue_name.unique(),model)
    pdbs[model] = pdb
    
 # find closest het, these alignments should be good enough for closest
for model in df['model'].unique():
    tempdf = df[df.model==model].copy()
    for i,atom in tempdf.iterrows():
        coord = pdbs[model].df["HETATM"][pdbs[model].df["HETATM"].residue_number == atom.residue_number][['x_coord','y_coord','z_coord']].to_numpy()
        if len(coord)>1:
            print("ERROR",coord)
        for col in all_columns:
            pdbB = col.split()[-1]
            if atom[col]>0: #and pdbB!=model:
                hetdist = pdbs[pdbB].distance(xyz=coord, records=('HETATM',))
                min_dist = hetdist.min()
                minatom = pdbs[pdbB].df["HETATM"].iloc[hetdist.idxmin()]
                while not (("wat" in col and minatom.residue_name=="HOH") or ("mg" in col and minatom.residue_name=="MG")):
                    #print('minatom',minatom.atom_name,col)#minatom)
                    hetdist[hetdist.idxmin()] = 1000
                    min_dist = hetdist.min()
                    minatom = pdbs[pdbB].df["HETATM"].iloc[hetdist.idxmin()]
                # fixing error of slightly misaligned MgH206
                if model == "1hr2_A" and pdbB == "1hr2_B" and atom.residue_number==356:
                    name = f"{minatom.residue_name}-363:{minatom.atom_name}"
                elif model == "1hr2_B" and pdbB == "1hr2_A" and atom.residue_number==367:
                    name = f"{minatom.residue_name}-357:{minatom.atom_name}"
                else:
                    name = f"{minatom.residue_name}-{minatom.residue_number}:{minatom.atom_name}"
                if min_dist > 3:
                    print(min_dist,name,model,col)
                if (minatom.residue_name == "MG" and "mg" not in col) or (minatom.residue_name == "HOH" and "wat" not in col):
                    print("ERR",model,atom.residue_number,name)
                tempdf.loc[i,col] = name
    for col in all_columns:
        pdbB = col.split()[-1]
        if pdbB!=model:
            if tempdf[col][tempdf[col]!=0].duplicated().any():
                print(pdbB,tempdf[[col,'residue_number']][tempdf[col]!=0][tempdf[col][tempdf[col]!=0].duplicated()])
                print(tempdf[col].unique())
    #    check no repeats
    new_df.append(tempdf)
df = pd.concat(new_df)

for model in df['model'].unique():
    compare_wat = df[(df.model==model)&(df.solvent=="HOH")].copy()
    compare_wat_bind = compare_wat["binding_site"].to_numpy()
    compare_wat_close = compare_wat["close_binding_site"].to_numpy()
    compare_mg = df[(df.model==model)&(df.solvent=="MG")].copy()
    compare_mg_bind = compare_mg["binding_site"].to_numpy()
    compare_mg_close = compare_mg["close_binding_site"].to_numpy()
    
    
    df[f"binding spot of mg in {model}"] = df.apply(lambda row: get_number_sym_overlap(row["close_binding_site"],compare_mg_bind,compare_mg_close,row["binding_site"]) ,axis=1)
                                                    #if (row["close_binding_site"]!="" and row["close_binding_site"]==row["close_binding_site"]) else get_number_exact_overlap(compare_mg_bind,row["binding_site"]),axis=1)
    df[f"binding spot of wat in {model}"] = df.apply(lambda row: get_number_sym_overlap(row["close_binding_site"],compare_wat_bind,compare_wat_close,row["binding_site"]) ,axis=1)
                                                     #if (row["close_binding_site"]!="" and row["close_binding_site"]==row["close_binding_site"]) else get_number_exact_overlap(compare_wat_bind,row["binding_site"]),axis=1)
    df[f"consensus of mg in {model}"] = df.apply(lambda row: get_consensus(row,model,df,'mg') ,axis=1)
    df[f"consensus of wat in {model}"] = df.apply(lambda row: get_consensus(row,model,df,'wat') ,axis=1)


df.to_csv('all_solvent_consensus_status_withconsensus.csv',index=False)

df22wat = df[(df.model=='2.2Å') & (df.solvent=="HOH")].copy()
df22mg = df[(df.model=='2.2Å') & (df.solvent=="MG")].copy()
df23wat = df[(df.model=='2.3Å') & (df.solvent=="HOH")].copy()
df23mg = df[(df.model=='2.3Å') & (df.solvent=="MG")].copy()

# both
df22wat['Consensus'] = ((df22wat['within 1A of wat in 2.3Å']!=0) & (df22wat['exact binding spot of wat in 2.3Å']>0))
df23wat['Consensus'] = ((df23wat['within 1A of wat in 2.2Å']!=0) & (df23wat['exact binding spot of wat in 2.2Å']>0))
df22mg['Consensus'] = ((df22mg['within 1A of mg in 2.3Å']!=0) & (df22mg['exact binding spot of mg in 2.3Å']>0))
df23mg['Consensus'] = ((df23mg['within 1A of mg in 2.2Å']!=0) & (df23mg['exact binding spot of mg in 2.2Å']>0))
df22wat['Consensus2'] = ((df22wat['within 1A of wat in 2.3Å']!=0) & (df22wat['binding spot of wat in 2.3Å']>0))
df23wat['Consensus2'] = ((df23wat['within 1A of wat in 2.2Å']!=0) & (df23wat['binding spot of wat in 2.2Å']>0))
df22mg['Consensus2'] = ((df22mg['within 1A of mg in 2.3Å']!=0) & (df22mg['binding spot of mg in 2.3Å']>0))
df23mg['Consensus2'] = ((df23mg['within 1A of mg in 2.2Å']!=0) & (df23mg['binding spot of mg in 2.2Å']>0))

print(df22wat[f"consensus of wat in 2.3Å"].sum())
print(df23wat[f"consensus of wat in 2.2Å"].sum())
print(df22mg[f"consensus of mg in 2.3Å"].sum())
print(df23mg[f"consensus of mg in 2.2Å"].sum())
print(df22wat[f"Consensus"].sum())
print(df23wat[f"Consensus"].sum())
print(df22mg[f"Consensus"].sum())
print(df23mg[f"Consensus"].sum())
print(df22wat[f"Consensus2"].sum())
print(df23wat[f"Consensus2"].sum())
print(df22mg[f"Consensus2"].sum())
print(df23mg[f"Consensus2"].sum())
#print(df22wat[df22wat[f"Consensus2"] & ~df22wat[f"consensus of wat in 2.3Å"]])