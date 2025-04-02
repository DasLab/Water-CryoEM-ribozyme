import argparse
from tqdm import tqdm
import pandas as pd
import numpy as np
from numpy import std, average, array, ones
from VolumeViewer import open_volume_file
from Segger import regions
from Matrix import chimera_xform
from FitMap import locate_maximum, overlap_and_correlation


'''
chimera --nogui --script "get_density_peaks.py \
    --dmapPath nt_wat.mrc \
    --segSigma 1 \
    --minNumPt 2 \
    --minSigma 2 \
    --outPath nt_wat_peaks.csv"  

segSigma 1 minNumPt 2 minSigma 2
    this gives 1097 peaks
    41/44% cryoEM within 1A
    64% cryoEM within 2A
    changing minNumPt to 1 did absolutely nothing ?? is this already assumed in teh segmentation
segSigma 0.5 minSigma 1
    this gives 2382 peaks
    50/52 1A
    84/85 2A
segSigma 0.1 minSigma 0.2
    this gives 2624 peaks
    50/53 1A
    85/85 2A
    --> essentially no change -- aka limited by the amount of segments?
segSigma -1 minSigma 0   (this is dens of 1, essentially everything...)
    this gives 2637 peaks
    50/53 1A
    85/85 2A

    > ${OUT}.out

to install pandas simple, chimera python not setup ssl
/home/rachael/.local/UCSF-Chimera64-1.14rc/bin/python -m pip install ~/Downloads/six-1.17.0-py2.py3-none-any.whl 
/home/rachael/.local/UCSF-Chimera64-1.14rc/bin/python -m pip install ~/Downloads/python_dateutil-2.9.0.post0-py2.py3-none-any.whl
/home/rachael/.local/UCSF-Chimera64-1.14rc/bin/python -m pip install ~/Downloads/numpy-1.7.2-cp27-cp27mu-manylinux1_x86_64.whl
/home/rachael/.local/UCSF-Chimera64-1.14rc/bin/python -m pip install ~/Downloads/pytz-2024.2-py2.py3-none-any.whl
/home/rachael/.local/UCSF-Chimera64-1.14rc/bin/python -m pip install ~/Downloads/pandas-0.18.1-cp27-cp27mu-manylinux1_x86_64.whl 
'''

# taken from segger, just cleaned up
def PtsToMax ( pts, dmap ) :

    apts = array ( pts, dtype=np.float32 )

    wts = ones ( len(pts), np.float32 )

    darray = dmap.data.matrix()

    xyz_to_ijk_tf = dmap.data.xyz_to_ijk_transform

    move_tf, stats = locate_maximum(apts, wts,
                                    darray, xyz_to_ijk_tf,
                                    max_steps = 1000,
                                    ijk_step_size_min = 0.01,
                                    ijk_step_size_max = 0.5,
                                    optimize_translation = True,
                                    optimize_rotation = True,
                                    metric = 'sum product',
                                    request_stop_cb = None)

    xf = chimera_xform ( move_tf )

    mpts = [None] * len(pts)
    for i, pt in enumerate(pts) :
        pt = chimera.Point ( apts[i][0], apts[i][1], apts[i][2] )
        xpt = xf.apply (pt)
        mpts[i] = xpt

    return mpts, stats['average map value']


# taken from segger, just cleaned up
def RegPt ( reg, segMap, mode="tomax" , doLargeMoveCheck=True) :

    P, ctr, val = None, None, None
    # what point in the region to use...
    if mode == "ctr" :

        # use the center of all points in the region
        # - may not be close to highest value or peak
        ctr = reg.center_of_points()
        P = chimera.Point(ctr[0],ctr[1],ctr[2])

    elif mode == "hval" :
        # use the highest value in the region
        # - may not be close to the center, but it is the peak

        ctr, maxD = None, -1e9
        rpts = reg.map_points()
        map_values = segMap.interpolated_values ( rpts, segMap.openState.xform )

        maxValPt = None
        for pt, val in zip(rpts, map_values) :
            if val > maxD :
                maxValPt, maxD = pt, val

        P = chimera.Point(maxValPt[0], maxValPt[1], maxValPt[2])
        ctr = [maxValPt[0], maxValPt[1], maxValPt[2]]
        val = maxD

    elif mode == "tomax" :
        # go to interpolated maximum...
        # - interestingly this can be different than the voxel with
        # - the highest value, due to interpolation used

        ctr, maxD = None, -1e9
        rpts = reg.map_points()
        map_values = segMap.interpolated_values ( rpts, segMap.openState.xform )

        maxValPt = None
        for pt, val in zip(rpts, map_values) :
            if val > maxD :
                maxValPt, maxD = pt, val

        ctr = maxValPt
        ctrP = chimera.Point(ctr[0], ctr[1], ctr[2])

        pts, avgMapV = PtsToMax ( [ctr], segMap )
        maxPt = pts[0]
        maxValP = chimera.Point ( maxPt[0], maxPt[1], maxPt[2] )

        # if movement is too large to maximum, likely this
        # is not a well-separated blob, so ignore it
        V = maxValP - ctrP
        if doLargeMoveCheck and V.length > 2*segMap.data.step[0]:
            P = None
            val = V.length
        else:
        #if V.length <= 2*segMap.data.step[0] :
            P = maxValP
            val = avgMapV
        #else :
        #    P = None
        #    val = V.length

    elif mode == "com" :
        # use the center of mass
        rpts = reg.map_points()

        map_values = segMap.interpolated_values ( rpts, segMap.openState.xform )

        ctr, sum = array ( [0,0,0] ), 0.0
        for pt, val in zip(rpts, map_values) :
            ctr += pt * val
            sum += val
        ctr = ctr / sum
        P = chimera.Point(ctr[0],ctr[1],ctr[2])

    return [P, val]

# taken from segger
def openMap ( mapPath ) :
    openedMap = None
    try :
        openedMap = open_volume_file ( mapPath, 'ccp4')[0]
        print " - opened map: %s" % mapPath
    except :
        print(" - could not open map: {}".format(mapPath))
        exit(1)
    return openedMap

def main():

    ###########################################################################
    # read in the arguments

    parser = argparse.ArgumentParser(description="get desnity peaks using segger.")
    parser.add_argument("--dmapPath", type=str, required=True, help="Path to the density map file.")
    parser.add_argument("--outPath", type=str, required=True, help="Path to the output csv of the peaks.")
    parser.add_argument("--segSigma", type=float, default=3.0, help="Contour level at which density is segmented.") 
    parser.add_argument("--minNumPt", type=int, default=2, help="Minimum numbers of voxels in a region to find the peak density.")      # 2 is what is in SWIM_script 
    parser.add_argument("--minSigma", type=float, default=5, help="Minimum density to count a peak.") 

    args = parser.parse_args()

    dmapPath = args.dmapPath
    segSigma = args.segSigma 
    minNumPt = args.minNumPt
    minSigma = args.minSigma
    outPath = args.outPath
   
    ###########################################################################
    # open files

    dmap = openMap ( dmapPath )

    ###########################################################################
    # segement maps

    M = dmap.data.full_matrix()
    
    sdev = std(M[M != 0]) # ignore perfect 0s when calculating stdev, important for the cutoff MD densities
    avg = average(M[M != 0])
    print(sdev,avg)
    mapThreshold = avg + segSigma * sdev 

    print " - for segement sigma of %.2f, threshold is %.4f" % (segSigma, mapThreshold)
    
    smod = regions.Segmentation(dmap.name, dmap)
    smod.calculate_watershed_regions ( dmap, mapThreshold )
    # this is all defaults, not sure what defualts are, seems like there are no other options?

    print " - got %d regions" % len ( smod.regions )

    ###########################################################################
    # get peaks

    regs = list(smod.regions)
    print " - processing %d segments..." % len(regs)

    n_regs = []
    data = []
    for ri, reg in enumerate ( regs ) :

        modi = 500 
        if ri % modi == 0 :
            print "Processing region %d/%d" % (ri+1, len(regs))

        npts = len(reg.points())
        #if reg.surface_piece:
        #    reg.hide_surface()
        if npts >= minNumPt : # this seems arbitary, so instead just process all?, was > 3, so >=2
            P, val = RegPt ( reg, dmap, mode="tomax",doLargeMoveCheck=False )
            sig = (val - avg) / sdev
            # P are the x,y,z coords
            if P != None and sig>minSigma:
                n_regs.append ( [val, sig, P, reg] )
                data.append([P[0],P[1],P[2],val,sig,npts])

    # sort regions by value
    n_regs.sort ( reverse=True, key=lambda x: x[0] )
    pd.DataFrame(data,columns=['x','y','z','dens','sig','npts']).to_csv(outPath, index=False)
    print ""
    print " - %d filtered segments" % len(n_regs)
    print " - saved ", outPath


if __name__ == "__main__":
    main()
