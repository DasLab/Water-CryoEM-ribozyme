import os, sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0,parent_dir)
from utils import *
import MDAnalysis
from MDAnalysis.analysis import align
from sys import argv


# input
output_bind = "all_solvent_bindingsite.csv"
output_cols = ["model","solvent","residue_number","close_binding_site","binding_site"]
out_data = "all_solvent_consensus_status.csv"
#mnumX = int(argv[1])
#mnumY = int(argv[2])
#models = df.model.unique()
ref = argv[1] #models[mnumX]
model = argv[2]#models[mnumY]
#print(model,ref)
model_dir = '../../../models/SWIM/models/'

PDB22_F = f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb'
PDB23_F = f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb'

MODEL_RENAME ={f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.2_5_0.7',

               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens2_Q0.7_Qres0.6_cycle999.pdb':'2.2_2_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens3_Q0.7_Qres0.6_cycle999.pdb':'2.2_3_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens4_Q0.7_Qres0.6_cycle999.pdb':'2.2_4_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens6_Q0.7_Qres0.6_cycle999.pdb':'2.2_6_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens7_Q0.7_Qres0.6_cycle999.pdb':'2.2_7_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens10_Q0.7_Qres0.6_cycle999.pdb':'2.2_10_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens15_Q0.7_Qres0.6_cycle999.pdb':'2.2_15_0.7',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens20_Q0.7_Qres0.6_cycle999.pdb':'2.2_20_0.7',

               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.5_Qres0.6_cycle999.pdb':'2.2_5_0.5',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.6_Qres0.6_cycle999.pdb':'2.2_5_0.6',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.8_Qres0.6_cycle999.pdb':'2.2_5_0.8',
               f'{model_dir}2.2A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.9_Qres0.6_cycle999.pdb':'2.2_5_0.9',

               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.3_5_0.7',

               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens2_Q0.7_Qres0.6_cycle999.pdb':'2.3_2_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens3_Q0.7_Qres0.6_cycle999.pdb':'2.3_3_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens4_Q0.7_Qres0.6_cycle999.pdb':'2.3_4_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens6_Q0.7_Qres0.6_cycle999.pdb':'2.3_6_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens7_Q0.7_Qres0.6_cycle999.pdb':'2.3_7_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens10_Q0.7_Qres0.6_cycle999.pdb':'2.3_10_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens15_Q0.7_Qres0.6_cycle999.pdb':'2.3_15_0.7',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens20_Q0.7_Qres0.6_cycle999.pdb':'2.3_20_0.7',

               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.5_Qres0.6_cycle999.pdb':'2.3_5_0.5',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.6_Qres0.6_cycle999.pdb':'2.3_5_0.6',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.8_Qres0.6_cycle999.pdb':'2.3_5_0.8',
               f'{model_dir}2.3A-G414_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.9_Qres0.6_cycle999.pdb':'2.3_5_0.9',

               f'{model_dir}2.2A-G414-halfA_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.2_halfA',
               f'{model_dir}2.2A-G414-halfB_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.2_halfB',
               f'{model_dir}2.2A-G414-nohalf_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.2_nohalf',
               f'{model_dir}2.3A-G414-halfA_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.3_halfA',
               f'{model_dir}2.3A-G414-halfB_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.3_halfB',
               f'{model_dir}2.3A-G414-nohalf_cdist3.2_distI4.5_maxwatD3.4_seg3_dens5_Q0.7_Qres0.6_cycle999.pdb':'2.3_nohalf'}
OTHER_STRUCTS = [x for x in MODEL_RENAME.keys() if x not in [PDB22_F,PDB23_F]]
MODEL_RENAME_R = {y:x for x,y in MODEL_RENAME.items()}


def local_align_and_get_dist(refernece,compare,atom_name,res_num,to_compare_sele,local_alignment_dist=10,save_pdb=None,min_dist_overlap=1):
    ref = MDAnalysis.Universe(refernece)
    comp = MDAnalysis.Universe(compare)
    # get local residues for local alignment
    current_Mg = f'resname {atom_name} and resnum {res_num}'
    resnums_to_align_ = list(ref.select_atoms(f'{RNA_SELECTION} and around {local_alignment_dist} ({current_Mg})').residues.resnums)
    resnums_to_align = resnums_to_align_.copy() 

    # get rid of non overlapping RNA, eg mutants or unmodeled nt
    for resnum in resnums_to_align_:
        x=ref.select_atoms(f'{RNA_SELECTION} and resnum {resnum}')
        y=comp.select_atoms(f'{RNA_SELECTION} and resnum {resnum}') 
        xx=list(x.names)
        yy=list(y.names)
        xx.sort()
        yy.sort()
        if len(x) == 0 or len(y) == 0:
            resnums_to_align.remove(resnum)
        elif x[0].resname != y[0].resname or len(x) != len(y) or xx != yy:
            resnums_to_align.remove(resnum)
    if len(resnums_to_align) == 0:
        return 0

    # now make the selection, make sure ordered
    resnums_to_align = "".join([f'resnum {x} or ' for x in resnums_to_align])[:-4]
    atom_neam_ = list(set(ref.select_atoms(f'{RNA_SELECTION} and ({resnums_to_align})').names))
    atom_neam_ +=list(set(comp.select_atoms(f'{RNA_SELECTION} and ({resnums_to_align})').names))
    atom_neam_.sort()
    ref_atoms = sum([ref.select_atoms(f'{RNA_SELECTION} and ({resnums_to_align}) and name {x}') for x in atom_neam_])
    comp_atoms = sum([comp.select_atoms(f'{RNA_SELECTION} and ({resnums_to_align}) and name {x}') for x in atom_neam_])
    rmsds = align.alignto(comp_atoms,ref_atoms)

    # get current solvent potisions, and calculate distance
    current_cord = ref.select_atoms(current_Mg).positions
    compare_cord = comp.select_atoms(to_compare_sele).positions
    diff = compare_cord.astype(float)-current_cord.astype(float)
    dist = np.linalg.norm(diff,axis=1)

    # get number of overlapping solvent
    num_close = sum(dist<=min_dist_overlap)

    # align the ref to other and save just the current solvent positoin, but now aligned.
    if save_pdb is not None:
        ref = MDAnalysis.Universe(refernece)
        comp = MDAnalysis.Universe(compare)
        current_Mg = f'resname {atom_name} and resnum {res_num}'
        resnums_to_align = ref.select_atoms(f'{RNA_SELECTION} and around {local_alignment_dist} ({current_Mg})').residues.resnums
        resnums_to_align = "".join([f'resnum {x} or ' for x in resnums_to_align])[:-4]
        rmsds = align.alignto(ref,  # mobile
                              comp,  # reference
                              select=f'{RNA_SELECTION} and ({resnums_to_align})',
                              match_atoms=True)
        current_cord = ref.select_atoms(current_Mg).positions
        ref = PandasPdb().read_pdb(refernece)
        new = PandasPdb()
        df = pd.DataFrame(columns=ref.df['ATOM'].columns,index=[0])
        df.atom_number = 1
        df.record_name = "ATOM"
        df.atom_name = atom_name
        df.residues_name = "A"
        df.chain_id = "N"
        df.residue_number = int(res_num)
        df.line_idx = 1
        df.b_factor = 0
        df.occupancy = 1
        df.x_coord = current_cord[0][0]
        df.y_coord = current_cord[0][1]
        df.z_coord = current_cord[0][2]
        df.element_symbol = "C"
        df.charge = 0
        df = df.fillna('')
        new.df["ATOM"] = df
        new.to_pdb(f'{save_pdb}_{atom_name}{res_num}.pdb')
    return num_close


# from the bindingsite of 2.2 and 2.3, get whether there is overlap
# also check out overlap after a local alignment 
df = pd.read_csv(output_bind)
tqdm.pandas()
compare_wat = df[(df.model==model)&(df.solvent=="HOH")].copy()
compare_wat_bind = compare_wat["binding_site"].to_numpy()
compare_wat_bind_close = compare_wat["close_binding_site"].to_numpy()

compare_mg = df[(df.model==model)&(df.solvent=="MG")].copy()
compare_mg_bind = compare_mg["binding_site"].to_numpy()
compare_mg_bind_close = compare_mg["close_binding_site"].to_numpy()

# for all ions and water get whether it is within 1A of current solvent
#if MODEL_RENAME_R[ref] == PDB22_F and MODEL_RENAME_R[model] == PDB23_F:
#    save_pdb = "../models/aligned_solvent/22in23"
#elif MODEL_RENAME_R[ref] == PDB23_F and MODEL_RENAME_R[model] == PDB22_F:
#    save_pdb = "../models/aligned_solvent/23in22"
#else:
save_pdb = None
dftemp =  df[df.model==ref].copy()

dftemp[f"within 1A of mg in {model}"] = dftemp.progress_apply(lambda row: local_align_and_get_dist(MODEL_RENAME_R[ref],MODEL_RENAME_R[model],row.solvent,row.residue_number,"resname MG",save_pdb=save_pdb),axis=1)
dftemp[f"within 1A of wat in {model}"] = dftemp.progress_apply(lambda row: local_align_and_get_dist(MODEL_RENAME_R[ref],MODEL_RENAME_R[model],row.solvent,row.residue_number,"resname HOH",save_pdb=save_pdb),axis=1)

dftemp.to_csv(f'{model}_{ref}_intermediate_{out_data}',index=False)


# for all ions and water get whether its binding site overlaps with current solvent
# if close_binding_site empty, fill with closest atom
dftemp[f"exact binding spot of mg in {model}"] = dftemp.apply(lambda row: get_number_exact_overlap(compare_mg_bind,row["close_binding_site"],compare_mg_bind_close) if (row["close_binding_site"]!="" and row["close_binding_site"]==row["close_binding_site"]) else get_number_exact_overlap(compare_mg_bind,row["binding_site"],compare_mg_bind_close),axis=1)
dftemp[f"exact binding spot of wat in {model}"] = dftemp.apply(lambda row: get_number_exact_overlap(compare_wat_bind,row["close_binding_site"],compare_wat_bind_close) if (row["close_binding_site"]!="" and row["close_binding_site"]==row["close_binding_site"]) else get_number_exact_overlap(compare_wat_bind,row["binding_site"],compare_wat_bind_close),axis=1)
dftemp[f"atom binding of wat in {model}"] = dftemp.apply(lambda row: get_number_partial_overlap(compare_wat_bind,row["binding_site"]),axis=1)
dftemp[f"atom binding of mg in {model}"] = dftemp.apply(lambda row: get_number_partial_overlap(compare_mg_bind,row["binding_site"]),axis=1)

dftemp.to_csv(f'{model}_{ref}_{out_data}',index=False)
