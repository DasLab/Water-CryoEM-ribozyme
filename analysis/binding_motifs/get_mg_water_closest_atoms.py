# import utils
import sys, os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0,parent_dir)
from utils import *

# input
out22 = "22A_mg_water_closest_atoms.csv"
out23 = "23A_mg_water_closest_atoms.csv"

# get closest atom
def get_closest_atom(pdb_f):
    # read pdb and keep only heavy atoms
    pdb = PandasPdb().read_pdb(pdb_f)
    pdb.df["ATOM"] = pdb.df["ATOM"][pdb.df["ATOM"].element_symbol!="H"].reset_index()
    
    mindists = []
    minatoms = []
    mintypes = []

    for i in range(len(pdb.df["HETATM"])):
        ref = pdb.df["HETATM"].loc[i,['x_coord','y_coord','z_coord']].to_numpy()

        # get closest RNA
        atomdist = pdb.distance(xyz=ref, records=('ATOM',))
        atommin = atomdist.min()

        # get closest solvent
        hetdist = pdb.distance(xyz=ref, records=('HETATM',))
        # get rid of compare to itself
        hetdist[hetdist==0] = np.inf
        # get rid of partial occupancy of same type and near
        if float(pdb.df["HETATM"].loc[i,"occupancy"])<1:
            hetdist[(hetdist<3.5) & (pdb.df["HETATM"].occupancy.astype(float)<1) & (pdb.df["HETATM"].residue_name == pdb.df["HETATM"].loc[i,"residue_name"])] = np.inf
        hetmin = hetdist.min()

        # get total minimum and label it
        if atommin<hetmin:
            mindist = atommin
            minatom = pdb.df["ATOM"].iloc[atomdist.idxmin()].atom_name
            if minatom in ["OP1","OP2"]:
                minatom = "OP"
            minatomtype = f"{minatom[0]}ᴿᴺᴬ"
        else:
            mindist = hetmin
            minatom = pdb.df["HETATM"].iloc[hetdist.idxmin()].atom_name
            if minatom=="O":
                minatomtype = "H₂O"
            else:
                minatomtype = "Mg²⁺"
        mindists.append(mindist)
        minatoms.append(minatom)
        mintypes.append(minatomtype)

    dist = pdb.df["HETATM"].copy()
    dist["min_dist"] = mindists
    dist["closest_atom"] = minatoms
    dist["closest_atom_type"] = mintypes
    return dist

dist22 = get_closest_atom(PDB22_F)#
dist23 = get_closest_atom(PDB23_F)#
dist23["bond"] = dist23.atom_name.apply(lambda x: "Mg²⁺" if x=="MG" else "H₂O") + "-" + dist23.closest_atom_type
dist22["bond"] = dist22.atom_name.apply(lambda x: "Mg²⁺" if x=="MG" else "H₂O") + "-" + dist22.closest_atom_type
dist22.to_csv(out22,index=False)
dist23.to_csv(out23,index=False)
