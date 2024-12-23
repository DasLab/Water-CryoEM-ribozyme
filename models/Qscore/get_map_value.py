import sys, os
import pandas as pd # /home/rachael/chimerax-1.5/bin/python3.9 -m pip install --target=/home/rachael/chimerax-1.5/lib/python3.9/site-packages pandas
import numpy as np
from chimerax.core.commands import run
from chimerax.atomic import all_atoms

mapPath = sys.argv[1]
molPath = sys.argv[2]
name = sys.argv[3]
outMolPath = "outputs/"+os.path.basename(molPath) + "__"+name+"__" + os.path.basename(mapPath) + ".csv"
run(session,f"open {mapPath}")
run(session,f"open {molPath}")
run(session,"measure mapvalues #1 atoms #2")

atoms = all_atoms(session)
mapvalue = [a.mapvalue for a in atoms]
atom_name = atoms.names
x_coord,y_coord,z_coord = atoms.coords.T
residue_number = atoms.residues.numbers
residue_name = atoms.residues.names
df = pd.DataFrame(np.array([atom_name,residue_number,residue_name,x_coord,y_coord,z_coord,mapvalue]).T,columns=['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','mapvalue'])
df.to_csv(outMolPath,index=False)
run(session,"exit")
