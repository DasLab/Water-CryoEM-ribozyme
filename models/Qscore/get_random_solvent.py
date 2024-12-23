# script to get random positions in the solvent shell
# python get_random_solvent.py DOTS_PER_A PDB OUTPUT
# python get_random_solvent.py 0.6 PDB23_F ../models/random_solvent_23.pdb
# python get_random_solvent.py 0.6 PDB22_F ../models/random_solvent_22.pdb

from utils import *


# inputs
dots_per_A = float(argv[1]) 
pdb = argv[2]
if pdb == 'PDB23_F':
    pdb = PDB23_F
elif pdb == 'PDB22_F':
    pdb = PDB22_F
output = argv[3]


def get_random_sampled_solvent_locations(pdb_f,out):
    df22 = PandasPdb().read_pdb(pdb_f)

    solvent_sample = []
    solvent_bad = []
    all_coords = df22.df["ATOM"][df22.df["ATOM"].element_symbol!="H"][['x_coord','y_coord','z_coord']]
    for i,coords in tqdm(all_coords.iterrows()):

        # get solvent area around atom
        coords = coords.to_numpy()
        minxg = ceil((coords[0]-WAT_MAX_DIST)*dots_per_A)/dots_per_A
        maxxg = floor((coords[0]+WAT_MAX_DIST)*dots_per_A)/dots_per_A
        minyg = ceil((coords[1]-WAT_MAX_DIST)*dots_per_A)/dots_per_A
        maxyg = floor((coords[1]+WAT_MAX_DIST)*dots_per_A)/dots_per_A
        minzg = ceil((coords[2]-WAT_MAX_DIST)*dots_per_A)/dots_per_A
        maxzg = floor((coords[2]+WAT_MAX_DIST)*dots_per_A)/dots_per_A

        # for all possible x,y,z
        sample_dist = 1/dots_per_A
        xs = np.arange(minxg,maxxg+sample_dist,sample_dist)
        ys = np.arange(minyg,maxyg+sample_dist,sample_dist)
        zs = np.arange(minzg,maxzg+sample_dist,sample_dist)
        for x in xs:
            for y in ys:
                for z in zs:
                    coord = [x,y,z]
                    dist = np.linalg.norm(coord-coords)
                    
                    # if too close, add as bad grid point
                    if dist<SOLVENT_MIN_DIST:
                        solvent_bad.append(coord)
                        while coord in solvent_sample:
                            solvent_sample.remove(coord)
                    # otherwise add to good grid points
                    elif (dist<WAT_MAX_DIST) and (coord not in solvent_bad) and (coord not in solvent_sample):
                        solvent_sample.append(coord)

    print(len(solvent_sample),len(solvent_bad))
    
    # make pdb with C to represent all of these locations
    new = PandasPdb()
    df = pd.DataFrame(columns=df22.df['ATOM'].columns)
    df.atom_number = range(1,len(solvent_sample)+1)
    df.record_name = "ATOM"
    df.atom_name = "CA"
    df.residues_name = "A"
    df.chain_id = "N"
    df.residue_number = 1
    df.line_idx = range(1,len(solvent_sample)+1)
    df.b_factor = 0
    df.occupancy = 1
    df.x_coord = np.array(solvent_sample)[:,0]
    df.y_coord = np.array(solvent_sample)[:,1]
    df.z_coord = np.array(solvent_sample)[:,2]
    df.element_symbol = "C"
    df.charge = 0
    df = df.fillna('')
    new.df["ATOM"] = df
    new.to_pdb(out)

get_random_sampled_solvent_locations(pdb,output)
