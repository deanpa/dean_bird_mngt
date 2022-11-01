#!/usr/bin/env python

import sys
import numpy as np
from rios import applier

infiles = applier.FilenameAssociations()
infiles.lcdb = sys.argv[1]

outfiles = applier.FilenameAssociations()
outfiles.newArray = sys.argv[2]

def indxLCDB(info, inputs, outputs):
    # empty array to population
    newClassArray = np.zeros_like(inputs.lcdb).astype(np.uint8)
    # specify classes as list
    nonHab = [1, 5, 6, 10, 12, 14, 16, 20, 21, 22, 45, 46, 70]
    grassScrub = [2, 15, 30, 33, 40, 41, 43, 44]
    shrub = [33, 47, 50, 51, 52, 55, 56, 58, 80, 81]
    exoticForest = [64, 71]
    broadLeaf = [54, 68]
    beech = [69]
    # put into a tuple and loop thru to reclass lcdb
    habClass = (nonHab, grassScrub, shrub, exoticForest, broadLeaf, beech)
    for i in range(1, len(habClass)):
        # find values in lcdb == habClass_i
        mask_i = np.in1d(inputs.lcdb, habClass[i]).reshape(np.shape(inputs.lcdb))
        newClassArray[mask_i] = i + 1
    outputs.newArray = newClassArray


applier.apply(indxLCDB, infiles, outfiles)

##  call on command line:
## ./reclassLCDB.py /home/dean/AA_Mega_Local/LCR_CurrentResearch/Walker_WarmForests/Modelling/Kea_Broadscale/GIS/NewRasters/lcdb_extent.tif /home/dean/workfolder/projects/dean_kiwi/kea/keaData/lcdb_kea_indx.img


##  ./reclassLCDB.py kiwi_data/lcdbBase_Fiord.tif kiwi_data/lcdbINDX.img
