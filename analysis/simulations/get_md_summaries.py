from utils import *
import MDAnalysis
from MDAnalysis.analysis.rms import RMSD


# inputs
mg_report_f = f"Mg_reports_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}.csv"
mg_report_cut_f = f"Mg_reports_dist{MG_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
na_report_f = f"Na_reports_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}.csv"
na_report_cut_f = f"Na_reports_dist{NA_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"
wat_report_f = f"wat_reports_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}.csv"
wat_report_cut_f = f"wat_reports_dist{WAT_MED_DIST}_time{TIME_FOR_BOUND}_skip{ALLOW_SKIP_FOR_BOUND}_rmsdcut{SIM_RMSD_CUTOFF}.csv"


mg_reports = []
na_reports = []
wat_reports = []

for sim in SIMULATIONS:
    print(f"Working on {sim}.")
    mg_report = pd.read_csv(f"../simulations/simulations_distances/{sim}_Mg.csv")
    na_report = pd.read_csv(f"../simulations/simulations_distances/{sim}_Na.csv")
    wat_report = pd.read_csv(f"../simulations/simulations_distances/{sim}_wat.csv")
    cols = mg_report.drop(columns=['other_atom', 'RNA_atom']).columns

    # remove any dists less than max_dist
    mg_report[cols] = mg_report[cols][mg_report[cols] <= MG_MED_DIST]
    na_report[cols] = na_report[cols][na_report[cols] <= NA_MED_DIST]
    wat_report[cols] = wat_report[cols][wat_report[cols] <= WAT_MED_DIST]

    # remove any rows that have total < time_limit hits
    mg_report = mg_report[mg_report.isnull().sum(axis=1) <= len(cols)-TIME_FOR_BOUND]
    na_report = na_report[na_report.isnull().sum(axis=1) <= len(cols)-TIME_FOR_BOUND]
    wat_report = wat_report[wat_report.isnull().sum(axis=1) <= len(cols)-TIME_FOR_BOUND]

    # remove any dists that do not fit time limit with allow_skip if applicable
    mg_report[cols] = mg_report[cols].apply(lambda row: pd.Series(correct_for_time_limit_allow_skip(row.to_list(),TIME_FOR_BOUND,ALLOW_SKIP_FOR_BOUND)),axis=1)
    na_report[cols] = na_report[cols].apply(lambda row: pd.Series(correct_for_time_limit_allow_skip(row.to_list(),TIME_FOR_BOUND,ALLOW_SKIP_FOR_BOUND)),axis=1)
    wat_report[cols] = wat_report[cols].apply(lambda row: pd.Series(correct_for_time_limit_allow_skip(row.to_list(),TIME_FOR_BOUND,ALLOW_SKIP_FOR_BOUND)),axis=1)

    # remove any rows that have total < time_limit hits
    mg_report = mg_report[mg_report.isnull().sum(axis=1) <= len(cols)-TIME_FOR_BOUND]
    na_report = na_report[na_report.isnull().sum(axis=1) <= len(cols)-TIME_FOR_BOUND]
    wat_report = wat_report[wat_report.isnull().sum(axis=1) <= len(cols)-TIME_FOR_BOUND]

    mg_report["simulation"] = sim
    na_report["simulation"] = sim
    wat_report["simulation"] = sim
    
    mg_reports.append(mg_report)
    na_reports.append(na_report)
    wat_reports.append(wat_report)
    
mg_reports = pd.concat(mg_reports)
na_reports = pd.concat(na_reports)
wat_reports = pd.concat(wat_reports)

mg_reports.to_csv(mg_report_f,index=False)
na_reports.to_csv(na_report_f,index=False)
wat_reports.to_csv(wat_report_f,index=False)

# get rmsd to 2.2 and 2.3 
u22 = MDAnalysis.Universe("../input_files/simulations/thr22.prmtop","../input_files/simulations/thr22.restrt")
u23 = MDAnalysis.Universe("../input_files/simulations/thr23.prmtop","../input_files/simulations/thr23.restrt")

dataframes = []
for i,sim in enumerate(SIMULATIONS):
    print(sim)
    u = MDAnalysis.Universe(f"../simulations/simulations/{sim}_nowater.thr.prmtop", f"../simulations/simulations/{sim}_md-nowater.nc")

    R22 = RMSD(u.select_atoms(CORE_SELECTION_RENUMBERED), u22.select_atoms(CORE_SELECTION_RENUMBERED),
                                      superposition=CORE_SELECTION_RENUMBERED)
    R23 = RMSD(u.select_atoms(CORE_SELECTION_RENUMBERED), u23.select_atoms(CORE_SELECTION_RENUMBERED),
                                      superposition=CORE_SELECTION_RENUMBERED)
    R22.run()
    R23.run()
    rmsd22 = R22.results.rmsd.T 
    rmsd23 = R23.results.rmsd.T 
    time = rmsd22[1]/1000
    rmsd22 = rmsd22[2]
    rmsd23 = rmsd23[2]
    data = pd.DataFrame(np.array([time,rmsd22,rmsd23]).T,columns=["Time (ns)","RMSD (Å)","rmsd23"])
    data['rmsd22_20avg'] = data["RMSD (Å)"].rolling(20,center=True).mean()
    data['rmsd23_20avg'] = data.rmsd23.rolling(20,center=True).mean()
    data['simulation'] = sim
    data["frame"] = data.index
    dataframes.append(data)
dataframes = pd.concat(dataframes)
dataframes["sim_frame"] = dataframes.simulation+dataframes.frame.astype(str)
dataframes.to_csv("22_23_core_rmsds.csv")

# filter frames by core rmsd cutoff and save
frames_to_cutdf = dataframes[(dataframes['RMSD (Å)']>SIM_RMSD_CUTOFF) & (dataframes['rmsd23']>SIM_RMSD_CUTOFF)]
frames_to_cut = {}
for sim in SIMULATIONS:
    frames_to_cut[sim] = frames_to_cutdf[frames_to_cutdf.simulation==sim].frame.to_list()
    
mg_reports = pd.read_csv(mg_report_f)
for sim in SIMULATIONS:
    mg_reports.loc[mg_reports.simulation==sim,[str(x) for x in frames_to_cut[sim]]] = np.nan
mg_reports.to_csv(mg_report_cut_f,index=False)

na_reports = pd.read_csv(na_report_f)
for sim in SIMULATIONS:
    na_reports.loc[na_reports.simulation==sim,[str(x) for x in frames_to_cut[sim]]] = np.nan
na_reports.to_csv(na_report_cut_f,index=False)

wat_reports = pd.read_csv(wat_report_f)
for sim in SIMULATIONS:
    wat_reports.loc[wat_reports.simulation==sim,[str(x) for x in frames_to_cut[sim]]] = np.nan
wat_reports.to_csv(wat_report_cut_f,index=False)
