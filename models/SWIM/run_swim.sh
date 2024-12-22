
SEG=3 # sigma threshold to segmenet at
DENS=5 # sigma threshold to filter at
Q=0.7 
QRES=0.6
CDIST_MG=3.0
CDIST_WAT=3.2
CYCLES=999 # should try only 1 cycle?

NAME=2.3A-G414
MAP=../../maps/Con1_2.3A_sh.mrc
HALFA=../../maps/Con1_2.3A_sh_half_A.mrc
HALFB=../../maps/Con1_2.3A_sh_half_B.mrc 
MODEL=inputs/9cbw_noH_noSol_414.pdb

OUT=models/${NAME}_cdist${CDIST_WAT}_seg${SEG}_dens${DENS}_Q${Q}_Qres${QRES}_cycle${CYCLES}.pdb

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
    --maxDistW 3.5 \
    --minCW $CDIST_WAT \
    --minCI $CDIST_MG \
    --cycles $CYCLES"  > ${OUT}.out

NAME=2.2A-G414
MAP=../../maps/Con2-2.2A_sh.mrc
HALFA=../../maps/Con2-2.2A_sh_half_A.mrc
HALFB=../../maps/Con2-2.2A_sh_half_B.mrc
MODEL=inputs/9cbu_noH_noSol_414.pdb

#MAP=../Greg_map_check/maps/Con2-2.2A_sh.mrc
#HALFA=../Greg_map_check/maps/Con2-2.2A_sh_half_A.ccp4
#HALFB=../Greg_map_check/maps/Con2-2.2A_sh_half_B.ccp4

OUT=models/${NAME}_cdist${CDIST_WAT}_seg${SEG}_dens${DENS}_Q${Q}_Qres${QRES}_cycle${CYCLES}.pdb

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
    --maxDistW 3.5 \
    --minCW $CDIST_WAT \
    --minCI $CDIST_MG \
    --cycles $CYCLES"  > ${OUT}.out
