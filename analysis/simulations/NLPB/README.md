# Use APBS to identify expected local ion concentrations

## Prepare

First, did minor edit to make the pdb file ready, namely to the 5' and 3' residues: changed RG5 and RG3, added HO5' and removed OP2 OP3 and P. Then was ready to make pqr using the following command:

```
pdb2pqr30 --userff=AMBER_DEShaw2018.DAT --usernames=AMBER_DEShaw2018.names --apbs-input=apbs.in --assign-only thr_for_pdb2pqr.pdb thr.pqr > thr_pdb2pqr.out
```
## APBS edits
In order to differentiate the ions, had to add a loop to the APBS code which would save the density for each individual ion the file edits around line 1099 (the VDT_NDENS loop) can be found here and is located in the APBS code at [src/mg/vpmg.c](vpmg.c). For each ion, this file should be edited for the loop to only add in one ion contributions,  APBS compiled as usual, and then apbs run.

## Run APBS
Run APBS as usual, for me I have to load in some python libraries, remember if wanting to save individual atom densities, need to implement the above APBS edits.

```
export LD_LIBRARY_PATH=${CONDA_ENV}/lib:$LD_LBRARY_PATH
${APBS_LOC}/build/UILD_MALLOC\=ON/src/apbs apbs_all.in
```
## Results
Due to size of files, the raw output are not found here, but the results, after converting from .dx to .mrc, can be found here, [ion_density_mg.mrc](ion_density_mg.mrc) and [ion_density_na.mrc](ion_density_na.mrc).