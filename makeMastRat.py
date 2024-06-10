#!/usr/bin/env python

import sys
import numpy as np
from osgeo import gdal
from rios import rat
from rios import applier

def riosCopyRaster(info, inputs, outputs):
    newRast = np.where(np.isnan(inputs.lcdb[0]), 0, inputs.lcdb[0] + 1)
    newRast = newRast.astype(np.uint32)
    outputs.newlcdb = np.expand_dims(newRast, 0)

inputs = applier.FilenameAssociations()
inputs.lcdb = sys.argv[1]

outputs = applier.FilenameAssociations()
outputs.newlcdb = sys.argv[2]

controls = applier.ApplierControls()
controls.setThematic(True)
controls.setOutputDriverName('KEA')

applier.apply(riosCopyRaster, inputs, outputs, controls=controls)

newRast = gdal.Open(outputs.newlcdb, gdal.GA_Update)

classNames = ["noData", "nonForest", "otherForest", "podocarpBroad", "pureBeech", "mixedBeech", "shrub"] 
rat.writeColumn(newRast, "Class_Names", classNames, colUsage=gdal.GFU_Name)
rat.writeColumn(newRast, "Masts", [0, 0, 0, 1, 1, 1, 0])


# ./makeMastRat.py eco5Hab.kea eco5_RAT.kea
