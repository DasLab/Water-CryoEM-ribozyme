# edited from segger/download/Segger/SWIM_script.py
# downloaded from git on 2024 Dec 17

import argparse
from tqdm import tqdm

# RCK added a save with Q-based Bfactors
def get_B_from_Q(Q):
    return 150 * (1-Q)

def main():

    ###########################################################################
    # read in the arguments

    parser = argparse.ArgumentParser(description="SWIM - Segment-guided Water and Ion Modeling, [placing water/ions using density maps.")

    parser.add_argument("--dmapPath", type=str, required=True, help="Path to the density map file.")
    parser.add_argument("--half1Path", type=str, default=None, help="Path to the first half-map.")
    parser.add_argument("--half2Path", type=str, default=None, help="Path to the second half-map.")

    parser.add_argument("--molPath", type=str, required=True, help="Path to the input molecular structure, without water or ions.")
    parser.add_argument("--outMolPath", type=str, required=True, help="Path to the output molecular structure (PDB format).")

    parser.add_argument("--segSigma", type=float, default=3.0, help="Contour level at which density is segmented.") # RCK added as was not present previously
    parser.add_argument("--thrSigma", type=float, default=3.0, help="Contour level above which water/ions are placed.")
    parser.add_argument("--useQ", type=bool, default=True, help="Use Q-score in placing water/ions.")
    parser.add_argument("--minQ", type=float, default=0.7, help="Minimum Q-score for placing water/ions.")
    parser.add_argument("--minQRes", type=float, default=0.6, help="Minimum Q-score for nearest residue.")  # RCK added as was not present previously
    parser.add_argument("--sigQ", type=float, default=0.6, help="Sigma value for Q-score calculation, use 0.6 for 1.5A or lower maps, 0.4 for 1.0 to 1.5A.")
    parser.add_argument("--toChain", type=str, default='', help="Chain to use for placement (leave empty for auto-selection).")
    parser.add_argument("--ionType", type=str, default="MG", help="Ion type to use (e.g., MG).")
    parser.add_argument("--minDistI", type=float, default=1.8, help="Minimum distance for ion placement in Angstroms.")
    parser.add_argument("--maxDistI", type=float, default=2.5, help="Maximum distance for ion placement in Angstroms.")
    parser.add_argument("--minDistW", type=float, default=2.5, help="Minimum distance for water placement in Angstroms.")
    parser.add_argument("--maxDistW", type=float, default=3.5, help="Maximum distance for water placement in Angstroms.")
    parser.add_argument("--minCW", type=float, default=3.0, help="Minimum distance for C to be declared clash to water in Angstroms.")
    parser.add_argument("--minCI", type=float, default=3.0, help="Minimum distance for C to be declared clash to ion in Angstroms.")

    parser.add_argument("--cycles", type=int, default=999, help="Number of cycles to continue, will automatically stop after convergence.") # RCK added as not present

    args = parser.parse_args()

    dmapPath = args.dmapPath
    half1Path = args.half1Path
    half2Path = args.half2Path
    molPath = args.molPath
    outMolPath = args.outMolPath
    segSigma = args.segSigma # RCK added
    thrSigma = args.thrSigma
    useQ = args.useQ
    minQ = args.minQ
    minQRes = args.minQRes # RCK added
    sigQ = args.sigQ
    toChain = args.toChain
    ionType = args.ionType
    minDistI, maxDistI = args.minDistI, args.maxDistI
    minDistW, maxDistW = args.minDistW, args.maxDistW
    min_C_dist_water = args.minCW # RCK added
    min_C_dist_ion = args.minCI # RCK added
    cycles = args.cycles # RCK added

    ###########################################################################
    # open files

    from VolumeViewer import open_volume_file

    def openMap ( mapPath ) :
        openedMap = None
        try :
            openedMap = open_volume_file ( mapPath, 'ccp4')[0]
            print " - opened map: %s" % mapPath
        except :
            print(" - could not open map: {}".format(mapPath))
            exit(1)
        return openedMap

    dmap = openMap ( dmapPath )

    hMapA, hMapB = None, None
    if half1Path :
        hMapA = openMap ( half1Path )
    if half2Path :
        hMapB = openMap ( half2Path )

    from chimera import openModels
    mol = openModels.open ( molPath, type='PDB' )[0]
    print " - opened model: %s" % molPath

    ###########################################################################
    # segement maps

    M = dmap.data.full_matrix()
    from numpy import std, average
    sdev = std(M)
    avg = average(M)
    mapThreshold = avg + segSigma * sdev # RCK edited

    print " - for segement sigma of %.2f, threshold is %.4f" % (segSigma, mapThreshold)

    from Segger import regions
    smod = regions.Segmentation(dmap.name, dmap)
    smod.calculate_watershed_regions ( dmap, mapThreshold )

    print " - got %d regions" % len ( smod.regions )

    ###########################################################################
    # calculate Qs
    # RCK the minQres in SWIM just seems to be ignore if Q has not already been calculated
    # RCK so calculate it first


    from Segger import qscores
    from Segger import gridm
    from Segger.SWIM import get_N_mg_binders
    import _multiscale
    g1 = gridm.Grid ()
    g1.FromAtomsLocal ( mol.atoms, 3.0 )#nearAtoms, 3.0 ) # RCK TEMP test
    minD, maxD = qscores.MinMaxD ( dmap )

    nearAtoms = [at for at in mol.atoms if not at.element.name == "H"]
    points = _multiscale.get_atom_coordinates ( nearAtoms, transformed = False )
    ptGrid = gridm.Grid()
    ptGrid.FromPoints ( points, 3.0 )

    #minDA, maxDA = qscores.MinMaxD ( hMapA )
    #minDB, maxDB = qscores.MinMaxD ( hMapB )
    #for at in tqdm ( nearAtoms ) :
    # RCK actually it needs residue values
    numNoQ = 0
    for r in tqdm(mol.residues) :
        qs = []
        for at in r.atoms :
            # TEST if it actually was H-averaged Q..., was not
            # but still for sanity of speed of run use models without H
            if at.element.name != "H" :
            #at.Q = qscores.QscoreG ( [at], dmap, sigQ, agrid=g1, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0)  # TODO test Qscore vs QscoreG
            #at.Q = qscores.Qscore ( [at], dmap, sigQ, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0)  # RCK I have no ideea why but this gives very different Qs
                #at.Q = qscores.Qscore_ ( [at], dmap, sigQ, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=minD, maxD=maxD, fitg=0)  # RCK also try G but without agrid?? or that distance different??
                # had dRAD as 1 to stay with SWIM, but those scores seem wrong, change to default 0.5 did not change anython
                # scores with Qscore_ seem much too high, and thhose with Qscore too low
                #at.Q = qscores.Qscore ( [at], dmap, sigQ, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=minD, maxD=maxD, fitg=0)

                atPt = at.coord()
                atPt = [atPt.x, atPt.y, atPt.z]
                xfI = dmap.openState.xform
                at.Q = qscores.QscorePt3 ( atPt, xfI, dmap, sigQ, ptGrid=ptGrid, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                # this copies what seems to be the latest (if 1) statement in the mapq code?

                at.bfactor = get_B_from_Q(at.Q)
                qs.append(at.Q)
        r.Q = sum(qs)/len(qs) # RCK as Greg has done before just do even weighting in averaging, NOT mass weighted
        #QA = qscores.Qscore_ ( r.atoms, hMapA, sigQ, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minDA, maxD=maxDA, fitg=0)
        #QB = qscores.Qscore_ ( r.atoms, hMapB, sigQ, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minDB, maxD=maxDB, fitg=0)
        print(r.id,r.Q,qs)#,r.Q,sum(QA)/len(QA),sum(QB)/len(QB))

    ###########################################################################
    # run SWIM

    from Segger import SWIM

    nearAtoms = [at for at in mol.atoms if not at.element.name == "H"]


    # RCK added as SWIM method describe iteratively adding water/ions
    # RCK note that nearAtoms stays same so we only look at RNA first solvation shell as we claim to do in the paper
    for i in range(cycles):
        # RCK edited as goSWIM params changed
        addW, addI = SWIM.goSWIM ( dmap, smod, mol, nearAtoms, 
            toChain = toChain, 
            minDistI = minDistI, maxDistI = maxDistI, 
            minDistW = minDistW, maxDistW = maxDistW, 
            hMapA = hMapA, hMapB = hMapB, 
            minQ = minQ, minQRes = minQRes, # RCK added minQRes
            sigQ = sigQ, ionType = ionType,
            # RCK added below as not previously present
            minWaterZ = thrSigma, minIonZ = thrSigma,
            min_C_dist_water = min_C_dist_water,
            min_C_dist_ion = min_C_dist_ion)

        # RCK added a save with Q-based Bfactors
        for at in addW+addI:
            #at.Q = qscores.QscoreG ( [at], dmap, sigQ, agrid=g1, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
            at.residue.Q = at.Q
            at.bfactor = get_B_from_Q(at.Q)
            if half1Path:
                print(at.residue.id,at.residue.type,at.Q,at.QhalfA,at.QhalfB)
            else:
                print(at.residue.id,at.residue.type,at.Q)

        print "Added %d waters, %d ions" % ( len(addW), len(addI) )
        if len(addW) + len(addI) == 0:
            print "Cycle %d no new water or ions added, stopping" % (i+1)
            break

    # when no new ion or waters or max cycles, stop and save
    print " - saving to: %s" % outMolPath
    import chimera
    chimera.PDBio().writePDBfile ( [mol], outMolPath )

    print ""


if __name__ == "__main__":
    main()
