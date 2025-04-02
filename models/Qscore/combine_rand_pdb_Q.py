# combine all the Q results from the split randomization runs
# python combine_rand_pdb_Q.py ../models/random_solvent_23_part 2.3 ../models/random_solvent_22_part  2.2 random_sovlent_shell_Q.csv random_sovlent_shell_halfQ.csv

from sys import argv

infileA = argv[1]
mapA = argv[2]
infileB = argv[3]
mapB = argv[4]
output = argv[5]
output2 = argv[6]

random = []
for csv in glob(f'{infileA}*__Q__*sh.*csv'):
    col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','Qscore']
    temp = pd.read_csv(csv,names=col_names)
    temp['residue_number'] = temp.residue_number.str[:-2].astype(int)
    temp['map'] = mapA
    random.append(temp)
for csv in glob(f'{infileB}*__Q__*sh.*csv'):
    col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','Qscore']
    temp = pd.read_csv(csv,names=col_names)
    temp['residue_number'] = temp.residue_number.str[:-2].astype(int)
    temp['map'] = mapB
    random.append(temp)
random = pd.concat(random)
random.to_csv(output, index=False)


random = []
for csv in glob(f'{infileA}*__Q__*sh_half_A*.csv'):
    csv2 = csv[:-10]+"B.ccp4.csv"
    col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','QscoreA']
    temp = pd.read_csv(csv,names=col_names)
    temp['residue_number'] = temp.residue_number.str[:-2].astype(int)
    temp['map'] = mapA
    col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','QscoreB']
    temp2 = pd.read_csv(csv2,names=col_names)
    temp['QscoreB'] = temp2.QscoreB
    random.append(temp)
for csv in glob(f'{infileB}*__Q__*sh_half_A*.csv'):
    csv2 = csv[:-10]+"B.ccp4.csv"
    col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','QscoreA']
    temp = pd.read_csv(csv,names=col_names)
    temp['residue_number'] = temp.residue_number.str[:-2].astype(int)
    temp['map'] = mapB
    col_names = ['atom_name','residue_number','residue_name','x_coord','y_coord','z_coord','QscoreB']
    temp2 = pd.read_csv(csv2,names=col_names)
    temp['QscoreB'] = temp2.QscoreB
    random.append(temp)
random = pd.concat(random)
random.to_csv(output2,index=False)
