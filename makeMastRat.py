#!/usr/bin/env python

import sys
import os
import numpy as np
from osgeo import gdal
from rios import rat
from rios import applier

def riosCopyRaster(info, inputs, outputs):
    newRast = np.where(np.isnan(inputs.lcdb[0]), 0, inputs.lcdb[0] + 1)
    newRast = newRast.astype(np.uint32)
    outputs.newlcdb = np.expand_dims(newRast, 0)

def mastRAT(infile, outfile):
    inputs = applier.FilenameAssociations()
    inputs.lcdb = infile

    outputs = applier.FilenameAssociations()
    outputs.newlcdb = outfile

    controls = applier.ApplierControls()
    controls.setThematic(True)
    controls.setOutputDriverName('KEA')

    applier.apply(riosCopyRaster, inputs, outputs, controls=controls)

#    newRast = gdal.Open(outputs.newlcdb, gdal.GA_Update)

    classNames = ["noData", "nonForest", "otherForest", "podocarpBroad", "pureBeech", 
        "mixedBeech", "shrub"] 
    rat.writeColumn(outputs.newlcdb, "Class_Names", classNames, colUsage=gdal.GFU_Name)
    rat.writeColumn(outputs.newlcdb, "RodentHab", [0, 0, 1, 1, 1, 1, 1])
    rat.writeColumn(outputs.newlcdb, "Masts", [0, 0, 0, 1, 1, 1, 0])


# ./makeMastRat.py eco5Hab.kea eco5_RAT.kea


########            Main function
#######
def main():
    baseDir = os.getenv('BIRDPROJDIR', default = '.')
    # SET PATH TO DATA - SHOULD EXIST ALREADY 
    inputDataPath = os.path.join(baseDir, 'SpeciesProjects', 'Kea', 'Data')

    inEco5 = os.path.join(inputDataPath,'eco5Hab.kea')
    outRaster = os.path.join(inputDataPath,'eco5_RAT.kea')

    mastRAT(inEco5, outRaster)
    

if __name__ == '__main__':
    main()


