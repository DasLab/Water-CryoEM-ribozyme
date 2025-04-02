# import utils
import sys, os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0,parent_dir)
from utils import *

# inputs
# rmsd cut effects residence time so for now not doing this
mg_in = f"../water_consensus/Mg_reports_bind_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
na_in = f"../water_consensus/Na_reports_bind_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
wat_in = f"../water_consensus/wat_reports_bind_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"

mg_2out = f"Mg_vals_bind_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
na_2out = f"Na_vals_bind_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
wat_2out = f"wat_vals_bind_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"

mg_out = f"Mg_restimes_bind_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
na_out = f"Na_restimes_bind_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
wat_out = f"wat_restimes_bind_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"


def get_all_atom_cat(bind):
    atoms = [get_atom_category_md(atom) for atom in bind.split()]
    atoms.sort()
    return ','.join(atoms)

# get bind sites for every time point in all simulations
def get_bind_site_time(dfin,out2,out):
    mg_reports = pd.read_csv(dfin)
    sims = []
    solvs = []
    binds = []
    presents = []
    for j,row in mg_reports.iterrows():
        present = []
        bind = []
        for i in range(400):
            bind_site = row[str(i)]
            if bind_site == bind_site:
                if bind_site not in bind:
                    if present == []: present.append([])
                    elif bind != []:
                        present.append([np.nan]*len(present[-1]))
                    bind.append(bind_site)
                for b,p in zip(bind,present):
                    if b == bind_site:
                        p.append(1)
                    else:
                        p.append(np.nan)
            else: 
                if present == []: present.append([])
                for x in present:
                    x.append(np.nan)
        if len(bind) != 0:
            presents.extend(present)
            sims.extend([row.simulation]*len(bind))
            solvs.extend([row.other_atom]*len(bind))
            binds.extend(bind)
    for i,present in enumerate(presents):
        present.append(sims[i])
        present.append(solvs[i])
        present.append(binds[i])

    df = pd.DataFrame(presents).rename(columns={400:'simulation',401:'other_atom',402:'bind_site'})
    df2 = df.copy()
    df2['atom_types'] = df2.bind_site.apply(lambda bind: get_all_atom_cat(bind))
    df2['bind_name'] = df2.bind_site.apply(lambda atoms: join_RNA_atom_names(atoms.split(),add_21=True))
    df2.to_csv(out2,index=False)
    keep_col = ['other_atom', 'bind_site','simulation']
    cols = df.drop(columns=keep_col).columns
    df['res_time'] = df[cols].apply(lambda row: get_resident_times(row.to_list().copy(),allow_skip=ALLOW_SKIP_FOR_BOUND,times=[]), axis=1)
    keep_col.append('res_time')
    df = df[keep_col].explode('res_time')
    df = df[~df.res_time.isnull()]
    df['atom_types'] = df.bind_site.apply(lambda bind: get_all_atom_cat(bind))
    df['bind_name'] = df.bind_site.apply(lambda atoms: join_RNA_atom_names(atoms.split(),add_21=True))
    df.to_csv(out,index=False)
    return df

get_bind_site_time(mg_in,mg_2out,mg_out)
get_bind_site_time(wat_in,wat_2out,wat_out)
get_bind_site_time(na_in,na_2out,na_out)

