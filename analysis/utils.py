from glob import glob
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
from math import floor,ceil
from tqdm import tqdm
from sys import argv


###############################################################################
# Models
###############################################################################

PDB23_F = '../models/2.3A-G414_cdist3.2_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb' #'../models/Con1_230717_all_SWIM.pdb'
PDB22_F = '../models/2.2A-G414_cdist3.2_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb' #'../models/Con2_230717_all_SWIM.pdb'
PDB23_ALIGNED = '../models/aligned_models/consensus_23.pdb'
OTHER_STRUCTS = glob('../models/other_models/*.pdb')
if "\\" in OTHER_STRUCTS[0]:
    MODEL_RENAME = {struct:struct.rsplit("/",1)[1].rsplit("\\",1)[1].rsplit(".",1)[0] for struct in OTHER_STRUCTS}
else:
    MODEL_RENAME = {struct:struct.rsplit("/",1)[1].rsplit(".",1)[0] for struct in OTHER_STRUCTS}
MODEL_RENAME[PDB23_F] = "2.3Å"
MODEL_RENAME[PDB22_F] = "2.2Å"
MODEL_RENAME_R = {y:x for x,y in MODEL_RENAME.items()}
SIMULATIONS = ["n2ed1","n2ed2","n2ed3","n2ed4",
        "n2em1","n2em2","n2em3","n2em4",
        "n2en1","n2en2","n2en3","n2en4",
        "n2rd1","n2rd2","n2rd3","n2rd4","n2rd5","n2rd6",
        "n2rn1","n2rn2","n2rn3","n2rn4","n2rn5","n2rn6",
        "n2rm1","n2rm2","n2rm3","n2rm4","n2rm5","n2rm6"]

MRC22_F = '../maps/Con2-2.2A_sh.mrc'
MRC23_F = '../maps/Con1_2.3A_sh.mrc'

# and Q of phase
QSCORE_22 = "../analysis/Qscore/Con2hf-4b_is2__Con2-2.2A_sh__thr2.000_Q0.70_0.60_rQ0.60__hA__hB__297-water__54-ion--.pdb__Q__Con2-2.2A_sh.mrc.csv"
QSCORE_23 = "../analysis/Qscore/Con3_is_4__Con1_2.3A_sh__thr2.000_Q0.70_0.60_rQ0.60__hA__hB__328-water__56-ion--.pdb__Q__Con1_2.3A_sh.ccp4.csv"
QSCORE_31 = "../analysis/Qscore/7ez0.pdb__Q__emd_31385.map.csv"

LOCALRES_22 = 'localres/Con2hf-4b_is2__Con2-2.2A_sh__thr2.000_Q0.70_0.60_rQ0.60__hA__hB__297-water__54-ion--.pdb__mapvalue__2.2A-cryosparc_P5_J1243_map_locres.mrc.csv'
LOCALRES_23 = 'localres/Con3_is_4__Con1_2.3A_sh__thr2.000_Q0.70_0.60_rQ0.60__hA__hB__328-water__56-ion--.pdb__mapvalue__2.3A-cryosparc_P5_J1242_map_locres.mrc.csv'

###############################################################################
# Geometric, simulation criteria
###############################################################################

MG_MAX_DIST = 2.5
MG_MED_DIST = 2.2
WAT_MAX_DIST = 3.5 
WAT_MED_DIST = 3.2
NA_MED_DIST = 2.5
SOLVENT_MIN_DIST = 1.5
TIME_FOR_BOUND = 10
ALLOW_SKIP_FOR_BOUND = 1
SIM_RMSD_CUTOFF = 3.4
SWIM_QSCORE_22 = 0.7
SWIM_QSCORE_23 = 0.7
SWIM_MAPVALUE = 5

###############################################################################
# Selections
###############################################################################

RNA_SELECTION = "not name H* and not name Na* and not name *M* and not name *Z* and not resname WAT and not resname HOH and not (name O5' and resnum 22) and not resnum 414"
MG_SELECTION = "name MG or name nMg or name mMg"
WAT_SELECTION = "name O and (resname WAT or resname HOH)"
NA_SELECTION = "name Na+ or name Na"

# define the stems and core
HELICES = {}
for i in range(369,402+1):
    HELICES[i] = "P9.2"
for i in list(range(353,367+1))+list(range(332,346+1)):
    HELICES[i] = "P9.1"
for i in range(318,331+1):
    HELICES[i] = "P9"
for i in list(range(315,316+1))+list(range(405,406+1)):
    HELICES[i] = "P9a"
for i in list(range(200,205+1))+list(range(116,121+1)):
    HELICES[i] = "P5"
for i in list(range(307,314+1))+list(range(262,268+1)):
    HELICES[i] = "P7"
for i in list(range(272,278+1))+list(range(96,102+1)):
    HELICES[i] = "P3"
for i in range(279,299+1):
    HELICES[i] = "P8"
for i in list(range(73,73+1))+list(range(75,81+1))+list(range(347,352+1)):
    HELICES[i] = "P13"
for i in list(range(31,42+1))+list(range(46,56+1)):
    HELICES[i] = "P2"
for i in list(range(58,68+1))+list(range(82,93+1)):
    HELICES[i] = "P2.1"
for i in range(22,27+1):
    HELICES[i] = "IGS"
for i in list(range(42,45+1))+list(range(169,172+1)):
    HELICES[i] = "P14"
for i in range(215,258+1):
    HELICES[i] = "P6"
for i in list(range(208,214+1))+list(range(107,112+1)):
    HELICES[i] = "P4"
for i in list(range(127,138+1))+list(range(178,195+1)):
    HELICES[i] = "P5a"
for i in range(139-164+1):
    HELICES[i] = "P5b"
for i in list(range(165,168+1))+list(range(173,177+1)):
    HELICES[i] = "P5c"
    
# define the core residues
CORE = ["P3", "P7", "P4", "P5", "P9", "P2", 'P9a']
non_base_paired = [210,262,313,314,324,325]
CORE_HELICES=[]
for i,h in HELICES.items():
    if h in CORE:
        CORE_HELICES.append(i)
CORE_HELICES.sort()
for i in non_base_paired:
    CORE_HELICES.remove(i)

str_selection = ""
for i in CORE_HELICES:
    str_selection += f"resnum {i} or "
CORE_SELECTION = f"{RNA_SELECTION} and ({str_selection[:-4]})"

CORE_HELICES_RENUMBER = [i-21 for i in CORE_HELICES]
str_selection = ""
for i in CORE_HELICES_RENUMBER:
    str_selection += f"resnum {i} or "
CORE_SELECTION_RENUMBERED = f"{RNA_SELECTION} and ({str_selection[:-4]})"


###############################################################################
# Formatting
###############################################################################

def join_RNA_atom_names(atoms,add_21=False):
    nums = []
    ats = []
    new_atoms = []

    for atom in atoms:
        at = atom.split("_")[2]
        nuc = atom.split("_")[1][0]
        num = atom.split("_")[0]
        if at == "OP1" or at == "OP2":
            at = "OP"
        if add_21:
            num = str(int(num)+21)
        nums.append(int(num))
        ats.append(at)
        new_atoms.append(nuc+"-"+num+":"+at)
    new_atoms = [x for _, __, x in sorted(zip(nums, ats, new_atoms),key=lambda pair: (pair[0],pair[1]))]
    return " ".join(new_atoms)

def get_number_within_1A(compare,coords):
    diff = compare.astype(float)-coords.astype(float)
    dist = np.linalg.norm(diff,axis=1)
    return sum(dist<=1)

def get_if_subset(siteA,siteB):
    if siteA != siteA or siteA=="":
        return False
    if siteB != siteB or siteB=="":
        return False
    siteA_list = siteA.split()
    for atom in siteB.split():
        if atom not in siteA_list:
            return False
    return True

def get_number_exact_overlap(compare,binding_site,compare_close):
    if binding_site != binding_site or binding_site == "":
        return 0
    overlap = [x for x,xc in zip(compare,compare_close) if get_if_subset(x,binding_site) and get_if_subset(binding_site,xc)]
    return len(overlap)

def get_number_exact_overlap2(compare,close_bind,bind):
    overlap = [x for x in compare if get_if_subset(x,close_bind) and get_if_subset(bind,x)]
    return len(overlap)

def get_if_one_agree(siteA,siteB):
    if siteA != siteA:
        return False
    siteA_list = siteA.split()
    for atom in siteB.split():
        if atom in siteA_list:
            return True
    return False

def get_number_partial_overlap(compare,binding_site):
    if binding_site != binding_site or binding_site == "":
        return 0
    overlap = [x for x in compare if get_if_one_agree(x,binding_site)]
    return len(overlap)

def replace_dist_with_rna_atom(row):
    row[~row.isnull()] = row.RNA_atom
    return row

# find fracs value that is subset the max frac
def frac_best(fracs,binding_site):
    if binding_site != binding_site or binding_site == "":
        return 0
    overlap = [x for x in fracs.keys() if get_if_subset(x,binding_site)]
    if len(overlap) == 0:
        return 0
    values = [fracs[x] for x in overlap]
    return max(values)

def get_category_both_one_none(row,modelA,modelB,name):
    if row[modelA+" exact match"] and row[modelB+" exact match"]:
        return f"{name} binding site found in both cryoEM maps."
    elif row[modelA+" exact match"] or row[modelB+" exact match"]:
        return f"{name} binding site found in only one cryoEM maps."
    else:
        return f"{name} binding site not found in any cryoEM maps."

def get_category_bind_overlap(row,model,mod_name,name):
    if row[f"exact binding spot of wat in {model}"]>0:
        return f"{name} with same bonded RNA atoms as {mod_name}"
    elif row[f"atom binding of wat in {model}"]>0:
        return f"{name} with at least one overlapping bonded RNA atoms as {mod_name}"
    else:
        return f"{name} with no overlap as {mod_name}"

def correct_for_time_limit_allow_skip(dists,time_lim,skip,ret_list=[],curr_list=[],curr_length=0,curr_skip=0):
    if dists==[]:
        if curr_length>=time_lim:
            return ret_list+curr_list
        else:
            return ret_list+([np.nan]*len(curr_list))
    else:
        if curr_list == []:
            if np.isnan(dists[0]):
                return correct_for_time_limit_allow_skip(dists[1:],time_lim,skip,
                                                         ret_list=ret_list+[np.nan],
                                                         curr_list=[],curr_length=0,curr_skip=0)
            else:
                return correct_for_time_limit_allow_skip(dists[1:],time_lim,skip,
                                                         ret_list=ret_list,
                                                         curr_list=[dists[0]],curr_length=1,curr_skip=0)
        else:
            if curr_skip>=skip and np.isnan(dists[0]):
                if curr_length>=time_lim:
                    return correct_for_time_limit_allow_skip(dists[1:],time_lim,skip,
                                                             ret_list=ret_list+curr_list+[np.nan],
                                                             curr_list=[],curr_length=0,curr_skip=0)
                else:
                    return correct_for_time_limit_allow_skip(dists[1:],time_lim,skip,
                                                             ret_list=ret_list+([np.nan]*(len(curr_list)+1)),
                                                             curr_list=[],curr_length=0,curr_skip=0)
            elif np.isnan(dists[0]):
                return correct_for_time_limit_allow_skip(dists[1:],time_lim,skip,
                                                         ret_list=ret_list,
                                                         curr_list=curr_list+[np.nan],curr_length=curr_length,curr_skip=curr_skip+1)
            else:
                return correct_for_time_limit_allow_skip(dists[1:],time_lim,skip,
                                                         ret_list=ret_list,
                                                         curr_list=curr_list+[dists[0]],curr_length=curr_length+1,curr_skip=0)
   

def get_resident_times(values, allow_skip, current_skip=0, current_time=0, times=[]):
    if values == []:
        # can happen because of rmsd cut
        if current_time >= TIME_FOR_BOUND:
            times.append(current_time)
        return times
    elif np.isnan(values[0]):
        if current_skip<allow_skip:
            return get_resident_times(values[1:], allow_skip=allow_skip, current_skip=current_skip+1, current_time=current_time, times=times)
        else:
            if current_time >= TIME_FOR_BOUND:
                times.append(current_time)
            return get_resident_times(values[1:], allow_skip=allow_skip, current_skip=0, current_time=0, times=times)
    else:
        return get_resident_times(values[1:], allow_skip=allow_skip, current_skip=0, current_time=current_time+1, times=times) 

def get_atom_category_md(atom):
    num,nuc,name = atom.split('_')
    nuc = nuc[0]
    if name in ["OP","O2'"]:
        return name 
    elif name in ["O3'","O5'","O4'"]:
        return "other ribose oxygen"
    elif f'{nuc}_{name}' in ['G_O6','C_O2','U_O4','U_O2']:
        return 'base carbonyl'
    elif f'{nuc}_{name}' in ['A_N7','G_N7','A_N1','A_N3','G_N3','C_N3']:
        return 'base secondary amine'
    elif f'{nuc}_{name}' in ['G_N1','U_N3']:
        return 'base primary ketimine'
    elif f'{nuc}_{name}' in ['A_N6','C_N4','G_N2']:
        return 'base primary amine'
    elif f'{nuc}_{name}' in ['C_N1','U_N1','G_N9','A_N9']:
        return 'base tertiary amine'
    elif 'P' in name or 'C' in name:
        return 'carbon or phosphorous' 
    else:
        return f'{nuc}_{name}'

def save_pdb_with_bfactor_by_res(pdb_file,selection,df,name):
    import MDAnalysis
    u22 = MDAnalysis.Universe(pdb_file)
    bfactors = np.zeros(len(u22.atoms))
    for resnum,res_quant in zip(np.unique(u22.select_atoms(selection).resnums),df):
        index = u22.select_atoms(f"resnum {resnum}").indices
        bfactors[index] = res_quant
    u22.atoms.tempfactors = bfactors
    u22.atoms.write(name)


def get_consensus_wat_mg(df,consensus_criteria):
    if True:
        matching_22_water = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH") & (df["consensus of wat in 2.3Å"])].residue_number.to_list()
        matching_22_mg = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG") & (df["consensus of mg in 2.3Å"])].residue_number.to_list()
        matching_23_water = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH") & (df["consensus of wat in 2.2Å"])].residue_number.to_list()
        matching_23_mg = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG") & (df["consensus of mg in 2.2Å"])].residue_number.to_list()
    else:
        if consensus_criteria == "both":
            matching_22_water = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH") & ((df["binding spot of wat in 2.3Å"]>0) & (df["within 1A of wat in 2.3Å"]>0))].residue_number.to_list()
            matching_22_mg = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG") & ((df["binding spot of mg in 2.3Å"]>0) & (df["within 1A of mg in 2.3Å"]>0))].residue_number.to_list()
            matching_23_water = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH") & ((df["binding spot of wat in 2.2Å"]>0) & (df["within 1A of wat in 2.2Å"]>0))].residue_number.to_list()
            matching_23_mg = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG") & ((df["binding spot of mg in 2.2Å"]>0) & (df["within 1A of mg in 2.2Å"]>0))].residue_number.to_list()
            #matching_22_water = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH") & ((df["exact binding spot of wat in 2.3Å"]>0) & (df["within 1A of wat in 2.3Å"]>0))].residue_number.to_list()
            #matching_22_mg = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG") & ((df["exact binding spot of mg in 2.3Å"]>0) & (df["within 1A of mg in 2.3Å"]>0))].residue_number.to_list()
            #matching_23_water = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH") & ((df["exact binding spot of wat in 2.2Å"]>0) & (df["within 1A of wat in 2.2Å"]>0))].residue_number.to_list()
            #matching_23_mg = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG") & ((df["exact binding spot of mg in 2.2Å"]>0) & (df["within 1A of mg in 2.2Å"]>0))].residue_number.to_list()
        elif consensus_criteria == "binding_site":
            matching_22_water = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH") & (df["exact binding spot of wat in 2.3Å"]>0)].residue_number.to_list()
            matching_22_mg = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG") & (df["exact binding spot of mg in 2.3Å"]>0)].residue_number.to_list()
            matching_23_water = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH") & (df["exact binding spot of wat in 2.2Å"]>0)].residue_number.to_list()
            matching_23_mg = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG") & (df["exact binding spot of mg in 2.2Å"]>0)].residue_number.to_list()
        elif consensus_criteria == "1A":
            matching_22_water = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH") & (df["within 1A of wat in 2.3Å"]>0)].residue_number.to_list()
            matching_22_mg = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG") & (df["within 1A of mg in 2.3Å"]>0)].residue_number.to_list()
            matching_23_water = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH") & (df["within 1A of wat in 2.2Å"]>0)].residue_number.to_list()
            matching_23_mg = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG") & (df["within 1A of mg in 2.2Å"]>0)].residue_number.to_list()
        elif consensus_criteria == "either":
            matching_22_water = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="HOH") & ((df["exact binding spot of wat in 2.3Å"]>0) | (df["within 1A of wat in 2.3Å"]>0))].residue_number.to_list()
            matching_22_mg = df[(df.model.astype(str)=="2.2Å") & (df.solvent=="MG") & ((df["exact binding spot of mg in 2.3Å"]>0) | (df["within 1A of mg in 2.3Å"]>0))].residue_number.to_list()
            matching_23_water = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="HOH") & ((df["exact binding spot of wat in 2.2Å"]>0) | (df["within 1A of wat in 2.2Å"]>0))].residue_number.to_list()
            matching_23_mg = df[(df.model.astype(str)=="2.3Å") & (df.solvent=="MG") & ((df["exact binding spot of mg in 2.2Å"]>0) | (df["within 1A of mg in 2.2Å"]>0))].residue_number.to_list()
    return matching_22_water,matching_22_mg,matching_23_water,matching_23_mg


def get_bfactor_based_on_overlap(row,structs,atom="mg"):
    # consensus 
    consensus_bool = row[f"consensus of {atom} in {structs[0]}"] #(row[f"exact binding spot of {atom} in {structs[0]}"]>0) & (row[f"within 1A of {atom} in {structs[0]}"]>0)
    for struct in structs[1:]:
        next_bool = row[f"consensus of {atom} in {struct}"] #(row[f"exact binding spot of {atom} in {struct}"]>0) & (row[f"within 1A of {atom} in {struct}"]>0)
        consensus_bool = consensus_bool | next_bool
    if consensus_bool:
        return 1.0
    # bound to same RNA atoms
    bound_bool = row[f"binding spot of {atom} in {structs[0]}"]>0 #row[f"exact binding spot of {atom} in {structs[0]}"]>0
    for struct in structs[1:]:
        next_bool = row[f"binding spot of {atom} in {struct}"]>0 #row[f"exact binding spot of {atom} in {struct}"]>0
        bound_bool = bound_bool | next_bool
    if bound_bool:
        return 0.75
    # near 
    near_bool = not row[f"within 1A of {atom} in {structs[0]}"] in [0,"0"]
    for struct in structs[1:]:
        next_bool = not row[f"within 1A of {atom} in {struct}"] in [0,"0"]
        near_bool = near_bool | next_bool
    if near_bool:
        return 0.5
    # atom bound 
    atom_bool = row[f"atom binding of {atom} in {structs[0]}"]>0
    for struct in structs[1:]:
        next_bool = row[f"atom binding of {atom} in {struct}"]>0
        atom_bool = atom_bool | next_bool
    if atom_bool:
        return 0.25
    return 0

def change_value_bfactor(row,u,col,name="solvent"):
    if row[name] == "HOH":
        index = u.select_atoms(f"{WAT_SELECTION} and resnum {row.residue_number}").indices
    elif row[name] == "MG":
        index = u.select_atoms(f"{MG_SELECTION} and resnum {row.residue_number}").indices
    else:
        index = u.select_atoms(f"{RNA_SELECTION} and resnum {int(row.residue_number)}").indices
    return index,row[col]

def get_density_interpolator(mrc_f):
    
    import mrcfile
    from scipy.interpolate import RegularGridInterpolator
    if '.dx'== mrc_f[-3:]:
        from gridData import Grid 
        grid = Grid(mrc_f)
        xs,ys,zs = grid.midpoints
        mapvalues = grid.grid
    else:
        mrc = mrcfile.open(mrc_f)
        xs = (np.arange(mrc.header.nxstart,mrc.header.mx)*mrc.voxel_size.x)+mrc.header.origin.x
        ys = (np.arange(mrc.header.nystart,mrc.header.my)*mrc.voxel_size.y)+mrc.header.origin.y
        zs = (np.arange(mrc.header.nzstart,mrc.header.mz)*mrc.voxel_size.z)+mrc.header.origin.z
        mapvalues = mrc.data
        # sometimes, but not always, its a 123 orientation instead off 321
        if mrc.header.mapc == 1:
            mapvalues = np.swapaxes(mapvalues,0,2)
    get_mapvalue = RegularGridInterpolator((xs,ys,zs),mapvalues, method='linear')
    return get_mapvalue

