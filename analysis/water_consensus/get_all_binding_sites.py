# import utils
import sys, os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0,parent_dir)
from utils import *

import MDAnalysis
from MDAnalysis.analysis import align

# input
output_bind = "all_solvent_bindingsite.csv"
output_cols = ["model","solvent","residue_number","close_binding_site","binding_site"]


# align 2.3 to 2.2 and save aligned pdb
u22 = MDAnalysis.Universe(PDB22_F)
u23 = MDAnalysis.Universe(PDB23_F)

rmsds = align.alignto(u23,  # mobile
                      u22,  # reference
                      select=CORE_SELECTION,match_atoms=True) 
u23.atoms.write(PDB23_ALIGNED)

# get all binding spots and coordinates
binders,close_binders = [], []
number, identity, models = [], [], []

for model in [PDB22_F,PDB23_F]+OTHER_STRUCTS:
    u22 = MDAnalysis.Universe(model)

    # get all waters and their binding sites
    waters = u22.select_atoms(WAT_SELECTION)
    for atom in waters:
        number.append(atom.resnum)
        identity.append(atom.resname)
        models.append(MODEL_RENAME[model]) 

        binding_site = u22.select_atoms(f"({RNA_SELECTION}) and around {WAT_MAX_DIST} (resname {atom.resname} and resnum {atom.resnum} and name {atom.name})")
        binders.append([str(atom.resnum)+"_"+atom.resname+"_"+atom.name if atom.name not in ["OP1","OP2"] else str(atom.resnum)+"_"+atom.resname+"_"+"OP" for atom in binding_site])
        binding_site = u22.select_atoms(f"({RNA_SELECTION}) and around {WAT_MED_DIST} (resname {atom.resname} and resnum {atom.resnum} and name {atom.name})")
        close_binders.append([str(atom.resnum)+"_"+atom.resname+"_"+atom.name if atom.name not in ["OP1","OP2"] else str(atom.resnum)+"_"+atom.resname+"_"+"OP" for atom in binding_site])
                
    # get all mg and their binding sites
    mgs = u22.select_atoms(MG_SELECTION)
    for atom in mgs:
        number.append(atom.resnum)
        identity.append(atom.resname)
        models.append(MODEL_RENAME[model]) 

        binding_site = u22.select_atoms(f"({RNA_SELECTION}) and around {MG_MAX_DIST} (resname {atom.resname} and resnum {atom.resnum} and name {atom.name})")
        binders.append([str(atom.resnum)+"_"+atom.resname+"_"+atom.name if atom.name not in ["OP1","OP2"] else str(atom.resnum)+"_"+atom.resname+"_"+"OP" for atom in binding_site])
        binding_site = u22.select_atoms(f"({RNA_SELECTION}) and around {MG_MED_DIST} (resname {atom.resname} and resnum {atom.resnum} and name {atom.name})")
        close_binders.append([str(atom.resnum)+"_"+atom.resname+"_"+atom.name if atom.name not in ["OP1","OP2"] else str(atom.resnum)+"_"+atom.resname+"_"+"OP" for atom in binding_site])
        
        
# save all information
binders = [join_RNA_atom_names(atoms) for atoms in binders]
close_binders = [join_RNA_atom_names(atoms) for atoms in close_binders]
df = pd.DataFrame(np.array([models,identity,number,close_binders,binders]).T,columns=output_cols)
df.to_csv(output_bind,index=False)
