#!/usr/bin/env python

import sys
import numpy
from osgeo import gdal
from rios import rat
from rios import applier

def addClass(info, inputs, outputs):
###    lowbeech = (inputs.lcdb[0] == 6) & (inputs.dem[0] < 200)
###    newlcdb = numpy.where(lowbeech, 7, inputs.lcdb[0]).astype(numpy.uint8)
#    # make Resolution and Secretary Islands non habitat
#    islandMask = (inputs.islands[0] == 1) & (newlcdb > 0)
#    newlcdb = numpy.where(islandMask, 1, newlcdb)
###    outputs.newlcdb = numpy.expand_dims(newlcdb, 0)

    newlcdb = inputs.lcdb[0].astype(numpy.uint8)
    outputs.newlcdb = numpy.expand_dims(newlcdb, 0)

inputs = applier.FilenameAssociations()
inputs.lcdb = sys.argv[1]
#inputs.dem = sys.argv[2]
#inputs.islands = sys.argv[3]

outputs = applier.FilenameAssociations()
outputs.newlcdb = sys.argv[2]


controls = applier.ApplierControls()
controls.setThematic(True)

applier.apply(addClass, inputs, outputs, controls=controls)


classNames = ["", "NonHabitat", "GrassScrub", "Shrub", "ProdForest", "Broadleaf", "Native"]
rat.writeColumn(outputs.newlcdb, "Class_Names", classNames, colUsage=gdal.GFU_Name)


rat.writeColumn(outputs.newlcdb, "Rodent_CC", [0.0, 0.0, 200.0, 150.0, 100.0, 300.0, 300])
rat.writeColumn(outputs.newlcdb, "Rodent_MastCC", [0.0, 0.0, 200.0, 150.0, 100.0, 300.0, 2387])
rat.writeColumn(outputs.newlcdb, "Rodent_CrashCC", [0.0, 0.0,  200.0, 150.0, 100.0, 300.0, 2187])
rat.writeColumn(outputs.newlcdb, "Masts", [0, 0, 0, 0, 0, 0, 1])
#   ./makeSeedDensity.py SpeciesProjects/Kea/Data/lcdb200_kea_classed.img SpeciesProjects/Kea/Data/dem200_kea.img SpeciesProjects/Kea/Data/seed_KeaTemp.img
#   ./makeSeedDensity.py SpeciesProjects/Kea/Data/lcdb200_kea_classed.img  SpeciesProjects/Kea/Data/seed_KeaTemp.img





#rat.writeColumn(outputs.newlcdb, "Rodent_CC", [0.0, 0.0,  75.0, 150.0, 300.0, 450.0, 450.0])
#rat.writeColumn(outputs.newlcdb, "Rodent_MastCC", [0.0, 0.0, 75.0, 150.0, 300.0, 5000.0, 5000.0])
#rat.writeColumn(outputs.newlcdb, "Rodent_CrashCC", [0.0, 0.0,  75.0, 150.0, 300.0, 50., 50.0])
#rat.writeColumn(outputs.newlcdb, "Masts", [0, 0, 0, 0, 0, 1, 1])


#   command line call:

#   ./makeSeedDensity.py SpeciesProjects/Kea/Data/lcdb_kea_extent.img SpeciesProjects/Kea/Data/dem200_kea.img SpeciesProjects/Kea/Data/seed_Kea.img

#   ./makeSeedDensity.py $POFPROJDIR/kiwi_data/lcdbINDX.img $POFPROJDIR/kiwi_data/dem.tif $POFPROJDIR/kiwi_data/islandSecRes.tif $POFPROJDIR/kiwi_data/rodentK.img


# ./addlowbeechclass.py $POFPROJDIR/kiwi_data/lcdbINDX.img $POFPROJDIR/kiwi_data/dem.tif $POFPROJDIR/kiwi_data/islandSecRes.tif $POFPROJDIR/kiwi_data/preBurninK.img $POFPROJDIR/kiwi_data/postBurninK.img

#   ./addlowbeechclass.py kiwi_data/lcdbINDX.img kiwi_data/dem.tif kiwi_data/islandSecRes.tif kiwi_data/lcdb_lowbeech.img
