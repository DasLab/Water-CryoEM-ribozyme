# Comparison of SWIM criteria

In response to Review #3.1 and #3.2 we analyzed the changes in number of waters identified and consensus between the two SWIM models from the two cryo-EM maps. We varied the minimum Q-score and the minimum density at the peak where the water is placed. We additionally investigated the effects of running the analysis without the half-map criteria or half-maps alone.

## Obtain models

The models were obtained using SWIM, the scripts can be found in [models/SWIM/run_swim_metrics.sh](../../../models/SWIM/run_swim_metrics.sh) and [models/SWIM/run_swim_half.sh](../../../models/SWIM/run_swim_half.sh). The models produced are found in [models/SWIM/models/](../../../models/SWIM/models/).

## Comparing the models

The models were compared by running the code below, analysis was then done in [anlayse_consensus.ipynb](anlayse_consensus.ipynb).

```
python get_all_binding_sites_SWIM_comparison.py

python bind_sites_parra.py 2.2_5_0.7 2.3_5_0.7

python bind_sites_parra.py 2.2_3_0.7 2.3_3_0.7
python bind_sites_parra.py 2.2_4_0.7 2.3_4_0.7
python bind_sites_parra.py 2.2_6_0.7 2.3_6_0.7
python bind_sites_parra.py 2.2_7_0.7 2.3_7_0.7

python bind_sites_parra.py 2.2_5_0.5 2.3_5_0.5
python bind_sites_parra.py 2.2_5_0.6 2.3_5_0.6
python bind_sites_parra.py 2.2_5_0.8 2.3_5_0.8
python bind_sites_parra.py 2.2_5_0.9 2.3_5_0.9

python bind_sites_parra.py 2.2_2_0.7 2.3_2_0.7 
python bind_sites_parra.py 2.2_10_0.7 2.3_10_0.7
python bind_sites_parra.py 2.2_15_0.7 2.3_15_0.7
python bind_sites_parra.py 2.2_20_0.7 2.3_20_0.7

python bind_sites_parra.py 2.2_halfA 2.2_halfB
python bind_sites_parra.py 2.3_halfA 2.3_halfB
python bind_sites_parra.py 2.2_nohalf 2.3_nohalf

python get_consensus_bind_sites.py
```

