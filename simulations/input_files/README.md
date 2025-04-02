## Input coordinates
Starting PDB can be found here, [con2_v2_20Mg.pdb](con2_v2_20Mg.pdb).

For nMg replace "MG   MG " with " nMg nMg".
For mMg replace "MG   MG " with " mMg mMg".

## Running MD - pmemd.cuda

Below is the full running script, the required input files are all listed below.
```
tleap -s -f tleap.in > tleap.out

pmemd.cuda -O -i mini_1.in -o mini_1.out -p thr.prmtop -c thr.rst7 -r mini_1.rst -ref thr.rst7
pmemd.cuda -O -i mini_2.in -o mini_2.out -p thr.prmtop -c mini_1.rst -r mini_2.rst -ref mini_1.rst
pmemd.cuda -O -i mini_3.in -o mini_3.out -p thr.prmtop -c mini_2.rst -r mini_3.rst -ref mini_2.rst
 
pmemd.cuda -i heat.in -o heat.out -p thr.prmtop -c mini_3.rst -r heat.rst -x heat.mdcrd -ref mini_3.rst -e heat.mden

pmemd.cuda -i eq_1.in -o eq_1.out -p thr.prmtop -c heat.rst -r eq_1.rst -x eq_1.mdcrd -ref heat.rst -e eq_1.mden
pmemd.cuda -i eq_2.in -o eq_2.out -p thr.prmtop -c eq_1.rst -r eq_2.rst -x eq_2.mdcrd -ref eq_1.rst -e eq_2.mden
pmemd.cuda -i eq_3.in -o eq_3.out -p thr.prmtop -c eq_2.rst -r eq_3.rst -x eq_3.mdcrd -ref eq_2.rst -e eq_3.mden
pmemd.cuda -i eq_4.in -o eq_4.out -p thr.prmtop -c eq_3.rst -r eq_4.rst -x eq_4.mdcrd -ref eq_3.rst -e eq_4.mden
pmemd.cuda -i eq_5.in -o eq_5.out -p thr.prmtop -c eq_4.rst -r eq_5.rst -x eq_5.mdcrd -ref eq_4.rst -e eq_5.mden
pmemd.cuda -i eq_6.in -o eq_6.out -p thr.prmtop -c eq_5.rst -r eq_6.rst -x eq_6.mdcrd -ref eq_5.rst -e eq_6.mden
pmemd.cuda -i eq_7.in -o eq_7.out -p thr.prmtop -c eq_6.rst -r eq_7.rst -x eq_7.mdcrd -ref eq_6.rst -e eq_7.mden
pmemd.cuda -i eq_8.in -o eq_8.out -p thr.prmtop -c eq_7.rst -r eq_8.rst -x eq_8.mdcrd -ref eq_7.rst -e eq_8.mden
pmemd.cuda -i eq_9.in -o eq_9.out -p thr.prmtop -c eq_8.rst -r eq_9.rst -x eq_9.mdcrd -ref eq_8.rst -e eq_9.mden
pmemd.cuda -i eq_10.in -o eq_10.out -p thr.prmtop -c eq_9.rst -r eq_10.rst -x eq_10.mdcrd -ref eq_9.rst -e eq_10.mden
pmemd.cuda -i eq_11.in -o eq_11.out -p thr.prmtop -c eq_10.rst -r eq_11.rst -x eq_11.mdcrd -ref eq_10.rst -e eq_11.mden
pmemd.cuda -i eq_12.in -o eq_12.out -p thr.prmtop -c eq_11.rst -r eq_12.rst -x eq_12.mdcrd -ref eq_11.rst -e eq_12.mden
pmemd.cuda -i eq_13.in -o eq_13.out -p thr.prmtop -c eq_12.rst -r eq_13.rst -x eq_13.mdcrd -ref eq_12.rst -e eq_13.mden
pmemd.cuda -i eq_14.in -o eq_14.out -p thr.prmtop -c eq_13.rst -r eq_14.rst -x eq_14.mdcrd -ref eq_13.rst -e eq_14.mden
pmemd.cuda -i eq_15.in -o eq_15.out -p thr.prmtop -c eq_14.rst -r eq_15.rst -x eq_15.mdcrd -ref eq_14.rst -e eq_15.mden
pmemd.cuda -i eq_16.in -o eq_16.out -p thr.prmtop -c eq_15.rst -r eq_16.rst -x eq_16.mdcrd -ref eq_15.rst -e eq_16.mden
pmemd.cuda -i eq_17.in -o eq_17.out -p thr.prmtop -c eq_16.rst -r eq_17.rst -x eq_17.mdcrd -ref eq_16.rst -e eq_17.mden
pmemd.cuda -i eq_18.in -o eq_18.out -p thr.prmtop -c eq_17.rst -r eq_18.rst -x eq_18.mdcrd -ref eq_17.rst -e eq_18.mden
pmemd.cuda -i eq_19.in -o eq_19.out -p thr.prmtop -c eq_18.rst -r eq_19.rst -x eq_19.mdcrd -ref eq_18.rst -e eq_19.mden
pmemd.cuda -i eq_20.in -o eq_20.out -p thr.prmtop -c eq_19.rst -r eq_20.rst -x eq_20.mdcrd -ref eq_19.rst -e eq_20.mden

pmemd.cuda -i md.in -o md.out -p thr.prmtop -c eq_20.rst -r md.rst7 -x md.nc -inf mdinfo -e md.mden
```

### Obtaining starting coordinates and addiving solvent and ions (tleap.in)

#### nMg and mMg amber parameter file (mMg_nMg.dat)
```
mMg and nMg parameters for use with TIP4p-d (https://github.com/bio-phys/optimizedMgFFs/blob/main/ff_Mg_tip4pd.itp 21-12-03) from paper Grotz, Cruz-León, and Schwierz JCTC 2021 RNA parameters from from parm14jpq.dat in amber20
MASS
mMg 24.305
nMg 24.305

NONBON
  mMg 0.538978 148.5430382
  nMg 0.544591 162.5262591

LJEDIT
  C  mMg              1.908000    0.086000  0.714416 48.484447
  CA mMg              1.908000    0.086000  0.714416 48.484447
  CB mMg              1.908000    0.086000  0.714416 48.484447
  CI mMg              1.908000    0.109400  0.714416 48.484447
  CP mMg              1.908000    0.086000  0.714416 48.484447
  CS mMg              1.908000    0.086000  0.714416 48.484447
  CQ mMg              1.908000    0.086000  0.714416 48.484447
  CT mMg              1.908000    0.109400  0.714416 48.484447
  C5 mMg              1.908000    0.086000  0.714416 48.484447
  C4 mMg              1.908000    0.086000  0.714416 48.484447
  H  mMg              0.600000    0.015700  0.714416 48.484447
  H1 mMg              1.387000    0.015700  0.714416 48.484447
  H2 mMg              1.287000    0.015700  0.714416 48.484447
  HA mMg              1.459000    0.015000  0.714416 48.484447
  H4 mMg              1.409000    0.015000  0.714416 48.484447
  H5 mMg              1.359000    0.015000  0.714416 48.484447
  HO mMg              0.000000    0.000000  0.714416 48.484447
  NA mMg              1.824000    0.170000  0.714416 48.484447
  NB mMg              1.824000    0.170000  0.714416 48.484447
  NC mMg              1.824000    0.170000  0.714416 48.484447
  N2 mMg              1.824000    0.170000  0.714416 48.484447
  N* mMg              1.824000    0.170000  0.714416 48.484447
  O  mMg              1.661200    0.210000  0.714416 48.484447
  O2 mMg              1.661200    0.210000  0.714416 48.484447
  OH mMg              1.721000    0.210400  0.714416 48.484447
  OS mMg              1.683700    0.170000  0.714416 48.484447
  P  mMg              2.100000    0.200000  0.714416 48.484447
  C  nMg              1.908000    0.086000  0.611630 46.0599418
  CA nMg              1.908000    0.086000  0.611630 46.0599418
  CB nMg              1.908000    0.086000  0.611630 46.0599418
  CI nMg              1.908000    0.109400  0.611630 46.0599418
  CP nMg              1.908000    0.086000  0.611630 46.0599418
  CS nMg              1.908000    0.086000  0.611630 46.0599418
  CQ nMg              1.908000    0.086000  0.611630 46.0599418
  CT nMg              1.908000    0.109400  0.611630 46.0599418
  C5 nMg              1.908000    0.086000  0.611630 46.0599418
  C4 nMg              1.908000    0.086000  0.611630 46.0599418
  H  nMg              0.600000    0.015700  0.611630 46.0599418
  H1 nMg              1.387000    0.015700  0.611630 46.0599418
  H2 nMg              1.287000    0.015700  0.611630 46.0599418
  HA nMg              1.459000    0.015000  0.611630 46.0599418
  H4 nMg              1.409000    0.015000  0.611630 46.0599418
  H5 nMg              1.359000    0.015000  0.611630 46.0599418
  HO nMg              0.000000    0.000000  0.611630 46.0599418
  NA nMg              1.824000    0.170000  0.611630 46.0599418
  NB nMg              1.824000    0.170000  0.611630 46.0599418
  NC nMg              1.824000    0.170000  0.611630 46.0599418
  N2 nMg              1.824000    0.170000  0.611630 46.0599418
  N* nMg              1.824000    0.170000  0.611630 46.0599418
  O  nMg              1.661200    0.210000  0.611630 46.0599418
  O2 nMg              1.661200    0.210000  0.611630 46.0599418
  OH nMg              1.721000    0.210400  0.611630 46.0599418
  OS nMg              1.683700    0.170000  0.611630 46.0599418
  P  nMg              2.100000    0.200000  0.611630 46.0599418

END
```

#### DESRES - random ion (x6)
```
source leaprc.water.tip4pew
source leaprc.RNA.Shaw

thr = loadpdb con2_v2_20Mg.pdb
check thr 

solvateoct thr TIP4PDBOX 9.0
addionsrand thr MG 75 
addionsrand thr Na+ 196

check thr
saveamberparm thr thr.prmtop thr.rst7
savePdb thr thr_MD.pdb
quit
```

#### DESRES - electrostatically-driven ion placement (x4)
```
source leaprc.water.tip4pew
source leaprc.RNA.Shaw

thr = loadpdb con2_v2_20Mg.pdb
#check thr 

solvateoct thr TIP4PDBOX 9.0
addions thr MG 75 
addions thr Na+ 196

# remove problem child Na+ 
# that is outside solvent sphere
remove thr thr.54688 
# add it randomly back in
addionsrand thr Na+ 1

check thr
saveamberparm thr thr.prmtop thr.rst7
savePdb thr thr_MD.pdb
quit
```

#### mMg - random ion (x6)
```
source leaprc.RNA.OL3
source leaprc.water.tip4pew
source leaprc.water.tip4pd
mnMG = loadAmberParams mMg_nMg.dat
i = createAtom mMg mMg 2.0 
set i element "Mg" 
r = createResidue mMg 
add r i 
mMg = createUnit mMg
add mMg r
thr = loadpdb con2_v2_20mMg.pdb
check thr 
solvateoct thr TIP4PDBOX 9.0

addionsrand thr mMg 75 
addionsrand thr Na+ 196

check thr
saveamberparm thr thr.prmtop thr.rst7
savePdb thr thr_MD.pdb
quit
```

#### mMg - electrostatically-driven ion placement (x4)
```
source leaprc.RNA.OL3
source leaprc.water.tip4pew
source leaprc.water.tip4pd
mnMG = loadAmberParams mMg_nMg.dat
i = createAtom mMg mMg 2.0 
set i element "Mg" 
r = createResidue mMg 
add r i 
mMg = createUnit mMg
add mMg r
thr = loadpdb con2_v2_20mMg.pdb
check thr 

addions thr mMg 75 
addions thr Na+ 196
solvateoct thr TIP4PDBOX 9.0

check thr
saveamberparm thr thr.prmtop thr.rst7
savePdb thr thr_MD.pdb
quit
```

#### nMg - random ion (x6)
```
source leaprc.RNA.OL3
source leaprc.water.tip4pew
source leaprc.water.tip4pd
mnMG = loadAmberParams mMg_nMg.dat
i = createAtom nMg nMg 2.0 
set i element "Mg" 
r = createResidue nMg 
add r i 
nMg = createUnit nMg
add nMg r

thr = loadpdb con2_v2_20nMg.pdb
check thr 

solvateoct thr TIP4PDBOX 9.0

addionsrand thr nMg 75 
addionsrand thr Na+ 196

check thr
saveamberparm thr thr.prmtop thr.rst7
savePdb thr thr_MD.pdb
quit
```

#### nMg - electrostatically-driven ion placement (x4)
```
source leaprc.RNA.OL3
source leaprc.water.tip4pew
source leaprc.water.tip4pd
mnMG = loadAmberParams mMg_nMg.dat
i = createAtom nMg nMg 2.0 
set i element "Mg" 
r = createResidue nMg 
add r i 
nMg = createUnit nMg
add nMg r

thr = loadpdb con2_v2_20nMg.pdb
check thr 

solvateoct thr TIP4PDBOX 9.0

addions thr nMg 75 
addions thr Na+ 196

check thr
saveamberparm thr thr.prmtop thr.rst7
savePdb thr thr_MD.pdb
quit
```


### minimization (mini_*.in)
```
&cntrl
imin=1,
ntx=1,
irest=0,
ntpr=50,
ntf=1,
ntb=1,
cut=9.0,
nsnb=10,
ntr=1,
maxcyc=1000,
ncyc=500,
ntmin=1,
restraintmask=':1-386',
restraint_wt=XXX,
&end
&ewald
ew_type = 0, skinnb = 1.0,
&end
```
XXX = (1) 25.0, (2) 20.0, (3) 15.0


### heat (heat.in)
```
&cntrl
imin=0, ntx=1, ntpr=500,
ntwr=500, ntwx=500, ntwe=500,
nscm=5000,
ntf=2, ntc=2,
ntb=1, ntp=0,
nstlim=100000, t=0.0, dt=0.002,
cut=9.0,
tempi=100.0, ntt=1,
ntr=1,nmropt=1,
restraintmask=':1-386',
restraint_wt=15.0,
&end
&wt type='TEMP0',istep1=0,
istep2=6250, value1=0.0,value2=100.0
&end
&wt type='TEMP0',istep1=6251,
istep2=68750, value1=100.0,value2=300.0,
&end
&wt type='TEMP0',istep1=68751,
istep2=100000, value1=300.0,value2=300.0,
&end
&wt type='END', &end
```

### equilibrate (eq_*.in)
```
&cntrl
imin=0, ntx=5, ntpr=500,
ntwr=500, ntwx=500,
ntwe=500,
nscm=5000,
ntf=2, ntc=2,
ntb=2, ntp=1, tautp=0.2,
taup=0.2,
nstlim=YYY,
t=0.0, dt=0.002,
cut=9.0,
ntt=1,
ntr=1,
irest=1,
restraintmask=':1-386',
restraint_wt=XXX,
&end
&ewald
ew_type = 0, skinnb = 1.0,
&end`
```

XXX = (1) 15, (2) 13, (3) 11, (4) 9, (5) 7, (6) 5, (7) 4, (8) 3, (9) 2, (10) 1, (11) 0.9, (12) 0.8, (13) 0.7, (14) 0.6, (15) 0.5, (16) 0.4, (17) 0.3, (18) 0.2, (19) 0.1, (20) 0.0
YYY = (1-10) 1000000 (11- ) 500000


### production (md.in)
```
&cntrl
imin=0, ntx=5,
ntpr=100000,ntwr=100000,
ntwx=100000,ntwe=100000,
nscm=100000,
ntf=2, ntc=2,
ntb=2, ntp=1, tautp=5.0, taup=5.0,
nstlim=200000000, t=0.0, dt=0.002,
cut=9.0,
ntt=1,
irest=1,
iwrap=1,
ioutfm=1,
&end
&ewald
ew_type = 0, skinnb = 1.0,
&end
```

## Obtaining trajectories - cpptraj
```
parm thr.prmtop
trajin md.nc 1 last 5 # every 1ns
autoimage
rms fit (:1-386)&!@H=

trajout YYY_md-full.nc netcdf # save trajectory in NetCDF format
run

closest 20000 (:1-386)&!@H= oxygen closestout YYY_closest_water.txt outprefix YYY_somewater
autoimage
rms fit2 (:1-386)&!@H=
trajout YYY_md-somewater.nc netcdf
run

strip :WAT outprefix YYY_nowater # get rid of waters
autoimage 
rms fit3 (:1-386)&!@H=
trajout YYY_md-nowater.nc netcdf # save trajectory
```
