from utils import *


# input
in_data = "all_solvent_consensus_status.csv"
mg_in = f"Mg_reports_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
mg_bind_out = f"Mg_reports_bind_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
mg_out = f"Mg_binding_spots_and_cooords_agreement_md_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
mg_md_out = f"Mg_md_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
na_in = f"Na_reports_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
na_bind_out = f"Na_reports_bind_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
wat_in = f"wat_reports_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
wat_bind_out = f"wat_reports_bind_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
water_out = f"wat_binding_spots_and_cooords_agreement_md_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
water_md_out = f"wat_md_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"


# get bind sites for every time point in all simulations
'''
mg_reports = pd.read_csv(mg_in)
cols = mg_reports.drop(columns=['other_atom', 'RNA_atom',"simulation"]).columns
mg_reports[["RNA_atom"]+list(cols)] = mg_reports[["RNA_atom"]+list(cols)].apply(lambda row: replace_dist_with_rna_atom(row), axis=1)
mg_reports = mg_reports.drop(columns="RNA_atom").groupby(["simulation","other_atom"]).agg(lambda x: ' '.join(x.dropna())).reset_index()
mg_reports.to_csv(mg_bind_out,index=False)

na_reports = pd.read_csv(na_in)
cols = na_reports.drop(columns=['other_atom', 'RNA_atom',"simulation"]).columns
na_reports[["RNA_atom"]+list(cols)] = na_reports[["RNA_atom"]+list(cols)].apply(lambda row: replace_dist_with_rna_atom(row), axis=1)
na_reports = na_reports.drop(columns="RNA_atom").groupby(["simulation","other_atom"]).agg(lambda x: ' '.join(x.dropna())).reset_index()
na_reports.to_csv(na_bind_out,index=False)

wat_reports = pd.read_csv(wat_in)
cols = wat_reports.drop(columns=['other_atom', 'RNA_atom',"simulation"]).columns
wat_reports[["RNA_atom"]+list(cols)] = wat_reports[["RNA_atom"]+list(cols)].apply(lambda row: replace_dist_with_rna_atom(row), axis=1)
wat_reports = wat_reports.drop(columns="RNA_atom").groupby(["simulation","other_atom"]).agg(lambda x: ' '.join(x.dropna())).reset_index()
wat_reports.to_csv(wat_bind_out,index=False)
'''

def get_summary_bind_sites(solvent,in_file,label_name,out_file,out_md_file):
    # get water summary
    df22_water = pd.read_csv(in_data)
    df22_water = df22_water[(df22_water.model.isin(['2.2Å','2.3Å'])) & (df22_water.solvent==solvent)]
    df22_water_bind = df22_water[df22_water.model=='2.2Å']["binding_site"].to_numpy()
    df23_water_bind = df22_water[df22_water.model=='2.3Å']["binding_site"].to_numpy()
    df22_water_closebind = df22_water[df22_water.model=='2.2Å']["close_binding_site"].to_numpy()
    df23_water_closebind = df22_water[df22_water.model=='2.3Å']["close_binding_site"].to_numpy()
    reports_all = pd.read_csv(in_file)
    cols = reports_all.drop(columns=['other_atom', 'simulation']).columns

    # for each simulation get count of binding sites
    all_fracs = []
    counts_exact = pd.Series(dtype='float64')
    for sim in SIMULATIONS:
        print(f"Working on {sim}.")
        reports = reports_all[reports_all.simulation==sim].copy()
        counts = pd.Series(dtype='float64')
        for col in cols:
            reports[str(col)] = reports[str(col)].apply(lambda x: x if pd.isnull(x) else join_RNA_atom_names(x.split(),add_21=True))
            curr_counts = reports[str(col)].value_counts()
            curr_counts[curr_counts>1] = 1 #restrict 1 per simulation
            counts = counts.add(curr_counts,fill_value=0)
            counts_exact = counts_exact.add(curr_counts,fill_value=0)

        fracs = counts/len(cols)
        all_fracs.append(fracs)

    # get average and max times and save
    fracs_max = pd.DataFrame(all_fracs).max()
    total_frames = len(cols)*len(reports_all.simulation.unique())
    fracs = counts_exact/total_frames
    df22_water_mdsim = pd.DataFrame([fracs,fracs_max]).T.reset_index().rename(columns={"index":"RNA_atoms",0:"avg_time",1:"max_time"})

    # get overlap of simulation binding site to 2.2A and 2.3A sites
    df22_water_mdsim["2.2A exact match"] = df22_water_mdsim.RNA_atoms.apply(lambda x: get_number_exact_overlap(df22_water_bind,x,df22_water_closebind))
    df22_water_mdsim["2.2A partial match"] = df22_water_mdsim.RNA_atoms.apply(lambda x: get_number_partial_overlap(df22_water_bind,x)) # TODO did not change
    df22_water_mdsim["2.3A exact match"] = df22_water_mdsim.RNA_atoms.apply(lambda x: get_number_exact_overlap(df23_water_bind,x,df23_water_closebind))
    df22_water_mdsim["2.3A partial match"] = df22_water_mdsim.RNA_atoms.apply(lambda x: get_number_partial_overlap(df23_water_bind,x))

    # get overlap of 2.2A and 2.3A binding sites to simulation binding sites
    mdsim_bind = df22_water_mdsim["RNA_atoms"].to_numpy()
    df22_water["MD partial match"] = df22_water.close_binding_site.apply(lambda x: get_number_partial_overlap(mdsim_bind,x))
    df22_water["MD exact match"] = df22_water.apply(lambda x: get_number_exact_overlap2(mdsim_bind,x.close_binding_site,x.binding_site),axis=1)

    # for all 2.2A and 2.3A binding sites get time something there in simulation
    df22_water["MD % exact"] = df22_water.close_binding_site.apply(lambda x: frac_best(fracs,x))
    df22_water["MD max% exact"] = df22_water.close_binding_site.apply(lambda x: frac_best(fracs_max,x))

    # get category name for md waters
    df22_water_mdsim["2A category"] = df22_water_mdsim.apply(lambda row: get_category_both_one_none(row,"2.2A","2.3A",label_name), axis=1)
    df22_water_mdsim.to_csv(out_md_file,index=False)

    # get category name for cryoem waters
    df22_water["2.3A category"] = df22_water.apply(lambda row: get_category_bind_overlap(row,"2.3Å","2.3Å",label_name), axis=1)
    df22_water["2.2A category"] = df22_water.apply(lambda row: get_category_bind_overlap(row,"2.2Å","2.2Å",label_name), axis=1)
    df22_water.to_csv(out_file,index=False)

    
get_summary_bind_sites(solvent="HOH",in_file=wat_bind_out,label_name="Water",out_file=water_out,out_md_file=water_md_out)
get_summary_bind_sites(solvent="MG",in_file=mg_bind_out,label_name="Mg",out_file=mg_out,out_md_file=mg_md_out)
