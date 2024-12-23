# code to split up the random solvent pdb into many

# python split_rand_pdb_for_processing.py ../models/random_solvent_22.pdb ../models/random_solvent_22_part
# python split_rand_pdb_for_processing.py ../models/random_solvent_23.pdb ../models/random_solvent_23_part

from sys import argv

infile = argv[1]
outloc = argv[2]

atoms = open(infile,'r').readlines()
f = open(f'{outloc}0.pdb','w')
num_per = 500
for i,atom in enumerate(atoms):
    if len(atom) == 84:
        atom = atom[:4]+atom[5:]
    if len(atom) != 83:
        print("ERROR")
    if (i+1)%num_per==0:
        print(i)
        f.close()
        f = open(f'{outloc}{(i+1)//num_per}.pdb','w')
    atom_num = atom.split()[1]
    new_num = str((int(atom_num)%num_per)+1)
    new_num = (" "*(6-len(new_num))) + new_num
    atom = atom[:5] + new_num + atom[11:]
    atom = atom[:22] + new_num[2:] + atom[26:]
    f.write(atom) 
f.close()

