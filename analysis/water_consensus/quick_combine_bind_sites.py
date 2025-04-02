# after running bind_sites_parra.py for all pairs, use this to combine all results

import pandas as pd
from glob import glob


for i in range(29):
    dfs = []
    for j in range(29):
        df = pd.read_csv(f'{i}_{j}_all_solvent_consensus_status.csv')
        dfs.append(df)
    dfs = pd.concat(dfs)
    dfs.to_csv(f'{i}_all_solvent_consensus_status.csv',index=False)

df = pd.read_csv(f'0_all_solvent_consensus_status.csv')
for i in range(1,29):
    df2 = pd.read_csv(f'{i}_all_solvent_consensus_status.csv')
    df = df.merge(df2,on=['model','solvent','residue_number','close_binding_site','binding_site'],validate='1:1')
df.to_csv(f'all_solvent_consensus_status.csv',index=False)

with open('../models/22in23.pdb','w') as new:
    for f in glob('../models/aligned_solvent/22in23_*.pdb'):
        old = open(f)
        new.write('\n'.join(old.readlines()))

with open('../models/23in22.pdb','w') as new:
    for f in glob('../models/aligned_solvent/23in22_*.pdb'):
        old = open(f)
        new.write('\n'.join(old.readlines()))

