
SEG=3 # sigma threshold to segmenet at
DENS=5 # sigma threshold to filter at
Q=0.7 
QRES=0.6
CDIST_MG=3.0
CDIST_WAT=3.2
DIST_I_I=4.5
CYCLES=999 # should try only 1 cycle?
MAXDISTW=3.4

NAME=2.3A-G414
MAP=../../maps/Con1_2.3A_sh.mrc
HALFA=../../maps/Con1_2.3A_sh_half_A.mrc
HALFB=../../maps/Con1_2.3A_sh_half_B.mrc 
MODEL=inputs/9cbw_noH_noSol_414.pdb

OUT=models/${NAME}_cdist${CDIST_WAT}_distI${DIST_I_I}_maxwatD${MAXDISTW}_seg${SEG}_dens${DENS}_Q${Q}_Qres${QRES}_cycle${CYCLES}.pdb

chimera --nogui --script "SWIM_script.py \
    --dmapPath $MAP \
    --half1Path $HALFA \
    --half2Path $HALFB \
    --molPath $MODEL \
    --outMolPath $OUT \
    --segSigma $SEG \
    --thrSigma $DENS \
    --minQ $Q \
    --minQRes $QRES \
    --sigQ 0.6 \
    --minDistI 1.8 \
    --maxDistI 2.5 \
    --minDistW 2.5 \
    --maxDistW $MAXDISTW \
    --minCW $CDIST_WAT \
    --minCI $CDIST_MG \
    --minDist_ion_ion $DIST_I_I \
    --cycles $CYCLES"  > ${OUT}.out

NAME=2.2A-G414
MAP=../../maps/Con2-2.2A_sh.mrc
HALFA=../../maps/Con2-2.2A_sh_half_A.mrc
HALFB=../../maps/Con2-2.2A_sh_half_B.mrc
MODEL=inputs/9cbu_noH_noSol_414.pdb

OUT=models/${NAME}_cdist${CDIST_WAT}_distI${DIST_I_I}_maxwatD${MAXDISTW}_seg${SEG}_dens${DENS}_Q${Q}_Qres${QRES}_cycle${CYCLES}.pdb

chimera --nogui --script "SWIM_script.py \
    --dmapPath $MAP \
    --half1Path $HALFA \
    --half2Path $HALFB \
    --molPath $MODEL \
    --outMolPath $OUT \
    --segSigma $SEG \
    --thrSigma $DENS \
    --minQ $Q \
    --minQRes $QRES \
    --sigQ 0.6 \
    --minDistI 1.8 \
    --maxDistI 2.5 \
    --minDistW 2.5 \
    --maxDistW $MAXDISTW \
    --minCW $CDIST_WAT \
    --minCI $CDIST_MG \
    --minDist_ion_ion $DIST_I_I \
    --cycles $CYCLES"  > ${OUT}.out
