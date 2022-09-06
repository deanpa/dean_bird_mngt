
import os
import csv
import numpy as np
import pickle
from osgeo import osr
from osgeo import ogr
ogr.UseExceptions()
from osgeo import gdal
gdal.UseExceptions()
from rios import rat
from numba import jit

COMPRESSED_HFA = ['COMPRESSED=YES']

sr = osr.SpatialReference()
sr.ImportFromEPSG(2193) # always nztm?
NZTM_WKT = sr.ExportToWkt()

class FormatError(Exception):
    "Was unable to read the format"

class KiwiData(object):
    """
    Holds the data used for the model run. This is calculated from 
    the parameters. Use the pickleSelf() method to save this to a pickle.
    """
    def __init__(self, params):
        self.params = params
        
        self.seasAdjRes = np.genfromtxt(self.params.seasAdjResFile, delimiter=",",names=True)

        self.rodentGeoTrans, self.rodentNcols, self.rodentNrows, originalExtent = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[0]))
        self.rodentExtentMask = originalExtent > 0

        # do the other extentMasks (stoats and kiwis)
        self.stoatGeoTrans, self.stoatNcols, self.stoatNrows, stoatMask = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[1]))

        # DEM at rodent resolution resample with 'Average'
        self.DEM = self.resampleRasterRodent(self.params.DEM, self.params.resolutions[0],
                    GDTMethod = gdal.GDT_Float32, GRAMethod = gdal.GRA_Average)

        self.kClasses = self.resampleRasterRodent(self.params.kClasses, self.params.resolutions[0],
                    GDTMethod = gdal.GDT_UInt16, GRAMethod = gdal.GRA_NearestNeighbour)

        self.kiwiGeoTrans, self.kiwiNcols, self.kiwiNrows, kiwiMask = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[2]))


        print('kiwinrows', self.kiwiNcols, self.kiwiNrows, type(self.kiwiNrows))


        # add the self.rodentMaxAltitude to the self.rodentExtentMask
        self.rodentExtentMask[self.DEM > self.params.rodentMaxAltitude] = 0

        ### AREAS TRAPPED IN RECENT TIMES -- at rodent resol
        self.islands = self.resampleRasterRodent(self.params.islands, 
                    self.params.resolutions[0],
                    GDTMethod = gdal.GDT_UInt16, GRAMethod = gdal.GRA_NearestNeighbour)

        self.islands[~self.rodentExtentMask] = 0

####        self.rodentIslandMask = (self.islands == 1)

        # Use to get stoat areas per zone -- still at rodent resolution 
        self.rodentExtentForStoats = self.rodentExtentMask.copy()

        # remove non-habitat from rodentExtentMask -- at rodent resolution
        self.rodentExtentMask[self.kClasses == 0] = 0


####        # remove rodents from Resolution and Secretary islands
####        self.rodentExtentMask[self.rodentIslandMask] = 0

        # get kiwiCorrectionK to scale pixels near water or high elevation
        (self.kiwiCorrectionK, self.stoatExtentMask, self.stoatPercentArea) = (
            scaleKiwiStoatMask(self.params.resolutions[0], 
            self.params.resolutions[2], self.kiwiNcols, self.kiwiNrows,
            self.DEM, originalExtent, self.params.kiwiMaxAltitude, 
            self.params.stoatMaxAltitude, self.rodentExtentForStoats))


        self.kiwiExtentMask = self.kiwiCorrectionK > 0.0
        self.stoatExtentMask = self.stoatExtentMask > 0

        self.rodentPercentArea = np.where(self.rodentExtentMask, 1.0, 0)
        self.kiwiKDummy = np.where(self.kiwiExtentMask, 1.0, 0)

        ## KIWI K MAP FOR SENSITIVITY TEST
        # self.kiwiKMap = getKiwiKMap(self.kiwiNcols, self.kiwiNrows,  
        #     self.kiwiCorrectionK, self.kiwiExtentMask, self.params.kiwiSurv, 
        #     self.params.kiwiSurvDecay, self.params.kiwiRecDecay,
        #     self.params.kiwiProd)
        self.kiwiKMap = getKiwiKMap(self.kiwiNcols, self.kiwiNrows,  
            self.kiwiCorrectionK, self.kiwiExtentMask, self.params.kiwiSurv, 
            self.params.kiwiSurvDDcoef, self.params.kiwiRecDDcoef,
            self.params.kiwiProd, self.params.kiwiTheta)

        kMapDS = gdal.Open(self.params.kClasses)
        ### MASK FOR CELLS THAT CAN MAST????
        self.mastingLU = rat.readColumn(kMapDS, "Masts") > 0
        ### 'CC' is the carrying capacity at per ha level.
        self.paramRodentCCLU = rat.readColumn(kMapDS, "Rodent_CC")
        self.paramRodentMastCCLU = rat.readColumn(kMapDS, "Rodent_MastCC")
        self.paramRodentCrashCCLU = rat.readColumn(kMapDS, "Rodent_CrashCC")

        # for assessing tracking tunnel rates in calculation.py
        self.beechMask = self.mastingLU[self.kClasses] & self.rodentExtentForStoats
    


        """
        driver = gdal.GetDriverByName('HFA')
        ds = driver.Create('islandTemp.img', self.rodentNcols, self.rodentNrows, 
                1, gdal.GDT_Byte)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.rodentGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.islands)
        del ds

        
        driver = gdal.GetDriverByName('HFA')
        ds = driver.Create('KiwiCorrectTemp.img', self.kiwiNcols, self.kiwiNrows, 
                1, gdal.GDT_Float32)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.kiwiGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.kiwiCorrectionK)
        del ds

        driver = gdal.GetDriverByName('HFA')
        ds = driver.Create('kiwiExtentMaskTemp.img', self.stoatNcols, self.stoatNrows, 
                1, gdal.GDT_Byte)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.stoatGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.kiwiExtentMask)
        del ds
     
        driver = gdal.GetDriverByName('GTiff')
        ds = driver.Create('demTemp.tif', self.rodentNcols, self.rodentNrows, 
                1, gdal.GDT_Float32)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.rodentGeoTrans)
        band = ds.GetRasterBand(1)
        bantuid.WriteArray(self.DEM)
        del ds


        driver = gdal.GetDriverByName('HFA')
        ds = driver.Create('rodentExtTemp.img', self.rodentNcols, self.rodentNrows, 
                1, gdal.GDT_Byte)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.rodentGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.rodentExtentMask)
        del ds

        
        ds = driver.Create('kclassTemp.img', self.rodentNcols, self.rodentNrows, 
                1, gdal.GDT_Byte)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.rodentGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.kClasses)
        del ds

        """
                

        (self.rodentControlList, self.rodentAreaDictByMgmt, self.pixelsInControlAndBeechMask,
                self.controlAndBeechMask) = self.readAndResampleControlForRodents()


#        print('pixelsInControlAndBeechMask', self.pixelsInControlAndBeechMask)


        (self.kiwiSpatialDictByMgmt, self.kiwiAreaDictByMgmt, 
            self.stoatAreaDictByMgmt) = self.readAndResampleControlForKiwis()


#        print('in dict', 'ALL' in self.rodentSpatialDictByMgmt.keys())
#        print('rodentSpatialDictByMgmt', self.rodentSpatialDictByMgmt.keys())
#        print('self.rodentControlDictByYear', self.rodentControlDictByYear.keys())
    


    def rasterizeShape(self, inshape, resolution, extent=None):
        """
        Rasterize a given shapefile mask at the chosen resolution. 
        Returns a numpy array with the data with zeros and ones.
        """
        ogrDataset = ogr.Open(inshape)
        if ogrDataset is None:
            msg = "Unable to read %s as a vector" % inshape
            raise FormatError(msg)
        ogrLayer = ogrDataset.GetLayer()

        if extent is None:
            x0, x1, y0, y1 = ogrLayer.GetExtent()
        else:
            x0, x1, y0, y1 = extent

        # modify extent to allow for nesting within stoat raster
        xmin = np.floor(x0 / self.params.resolutions[1]) * self.params.resolutions[1]
        xmax = np.ceil(x1 / self.params.resolutions[1]) * self.params.resolutions[1]
        ymin = np.floor(y0 / self.params.resolutions[1]) * self.params.resolutions[1]
        ymax = np.ceil(y1 / self.params.resolutions[1]) * self.params.resolutions[1]

        ncols = int((xmax - xmin) / resolution)
        nrows = int((ymax - ymin) / resolution)

        # Create temp file as a raster
        driver = gdal.GetDriverByName('HFA')
        gdalDataset = driver.Create('extentraster.img', ncols, nrows, 1,
                gdal.GDT_Byte, COMPRESSED_HFA)

        geoTrans = [xmin, resolution, 0, ymax, 0, -resolution]
        gdalDataset.SetGeoTransform(geoTrans)

        gdalDataset.SetProjection(NZTM_WKT)

        gdal.RasterizeLayer(gdalDataset, [1], ogrLayer, burn_values=[1])
        gdalDataset.FlushCache()

        gdalBand = gdalDataset.GetRasterBand(1)
        data = gdalBand.ReadAsArray()

        return geoTrans, ncols, nrows, data


    def resampleRasterRodent(self, infile, resolution, GDTMethod, GRAMethod):
        """
        Resample the input file to the given resolution and bounds
        Note: Assumes rodent res etc
        """
        driver = gdal.GetDriverByName('HFA')
        gdalDataset = driver.Create('resraster.img', self.rodentNcols, self.rodentNrows, 1,
                GDTMethod, COMPRESSED_HFA)
        gdalDataset.SetGeoTransform(self.rodentGeoTrans)

        gdalDataset.SetProjection(NZTM_WKT)

        inDS = gdal.Open(infile)
        src_wkt = inDS.GetProjection()

        gdal.ReprojectImage(inDS, gdalDataset, src_wkt, NZTM_WKT, 
                GRAMethod)

        outBand = gdalDataset.GetRasterBand(1)
        data = outBand.ReadAsArray()

        return data

    def resampleRasterKiwi(self, infile, resolution, GDTMethod, GRAMethod):
        """
        Resample the input file to the given resolution and bounds
        Note: Assumes kiwi res etc
        """
        driver = gdal.GetDriverByName('HFA')
        gdalDataset = driver.Create('resraster.img', self.kiwiNcols, self.kiwiNrows, 1,
                GDTMethod, COMPRESSED_HFA)
        gdalDataset.SetGeoTransform(self.kiwiGeoTrans)

        gdalDataset.SetProjection(NZTM_WKT)

        inDS = gdal.Open(infile)
        src_wkt = inDS.GetProjection()

        gdal.ReprojectImage(inDS, gdalDataset, src_wkt, NZTM_WKT, 
                GRAMethod)

        outBand = gdalDataset.GetRasterBand(1)
        data = outBand.ReadAsArray()

        return data

    def readAndResampleControlForRodents(self):
        """
        Go through the control file and rasterize to the rodent
        resolution. Returns a list with (mask, startYear, Revisit, inFilePath).

        Also returns a dictionary
        keyed on the original shape file name so that analysis
        can be performed on each management area.
        """
        controlList = [] # what we return - values are
        # (mask, startYear, Revisit, inFilePath)
        spatialDict = {} # keyed on path - has mask for each file so 
            # we don't need to resample again if we already have the path

        ### AREA DICTIONARY FOR CALCULATING DENSITY
        rodentAreaDict = {}

        ### Dict for the number of pixels in beech in each mgmt zone
        pixelsInControlAndBeechMask = {}
        controlAndBeechMask = {}
        # work out the extent to use
        x0, y1 = gdal.ApplyGeoTransform(self.rodentGeoTrans, 0, 0)
        x1, y0 = gdal.ApplyGeoTransform(self.rodentGeoTrans, self.rodentNcols, self.rodentNrows)
        extent = [x0, x1, y0, y1]

        firstRow = True
        with open(self.params.controlFile, newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if firstRow:
                    firstRow = False
                    continue

                if len(row) == 4:
                    shpFile, startYear, ctrlMth, revisit = row
                    startYear = int(startYear)
                    ctrlMth = int(ctrlMth)
                    revisit = int(revisit)

                    if self.params.controlPathPrefix is not None:
                        shpFile = os.path.join(self.params.controlPathPrefix, shpFile)

                    # check we already have this one resampled
                    # if so - just grab it.
                    if shpFile in spatialDict:
                        mask = spatialDict[shpFile]
                    else:
                        if not os.path.exists(shpFile):
                            raise IOError("Cannot find file %s" % shpFile)

                        try:
                            geoTrans, ncols, nrows, data = self.rasterizeShape(shpFile, 
                                    self.params.resolutions[0], extent=extent)
                        except FormatError:
                            # must be a raster input
                            data = self.resampleRasterRodent(shpFile, self.params.resolutions[0],
                                        GDTMethod = gdal.GDT_UInt16,
                                        GRAMethod = gdal.GRA_NearestNeighbour)
                        mask = (data == 1)

                        # save it for next time
                        spatialDict[shpFile] = mask

                    ### GET NUMBER OF RODENT PIXELS PER MGMT AREA
                    rodentAreaDict[shpFile] = np.sum(mask & self.rodentExtentMask)

                    ### USE THIS FOR ESTIMATING DENSITY IN BEECH
                    
                    controlAndBeechMask[shpFile] = self.beechMask & mask
                    pixelsInControlAndBeechMask[shpFile] = np.sum(controlAndBeechMask[shpFile])

                    # save data
                    controlList.append((mask, startYear, ctrlMth, revisit, shpFile))

        rodentAreaDict['ALL'] = np.sum(self.rodentExtentMask)
        pixelsInControlAndBeechMask['ALL'] = 0
        controlList.append((self.rodentExtentMask.copy(), 100, -1, 300, 'ALL'))
#        print('rodentAreaDict', rodentAreaDict)
        return(controlList, rodentAreaDict, pixelsInControlAndBeechMask, controlAndBeechMask)

    def readAndResampleControlForKiwis(self):
        """
        Similar to readAndResampleControlForRodents() above, but 
        resamples to Kiwi resolution and returns a dictionary
        keyed on the original shape file name so that analysis
        can be performed on each management area.
        """
        kiwiSpatialDict = {}  # keyed on file name
        kiwiAreaDict = {}
        stoatAreaDict = {}
        # work out the extent to use
        x0, y1 = gdal.ApplyGeoTransform(self.kiwiGeoTrans, 0, 0)
        x1, y0 = gdal.ApplyGeoTransform(self.kiwiGeoTrans, self.kiwiNcols, self.kiwiNrows)
        extent = [x0, x1, y0, y1]

        firstRow = True
        with open(self.params.controlFile, newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if firstRow:
                    firstRow = False
                    continue

                if len(row) == 4:
                    shpFile, startYear, ctrlMth, revisit = row
        
                    if self.params.controlPathPrefix is not None:
                        shpFile = os.path.join(self.params.controlPathPrefix, shpFile)

                    if shpFile not in kiwiSpatialDict:
                        # we haven't come accross this one before
                        if not os.path.exists(shpFile):
                            raise IOError("Cannot find file %s" % shpFile)

                        try:
                            geoTrans, ncols, nrows, data = self.rasterizeShape(shpFile, 
                                    self.params.resolutions[2], extent=extent)
                        except FormatError:
                            # must be a raster input
                            data = self.resampleRasterKiwi(shpFile, self.params.resolutions[2],
                                        GDTMethod = gdal.GDT_UInt16,
                                        GRAMethod = gdal.GRA_NearestNeighbour)
                        mask = (data == 1)
                        ### POPULATE kiwi and stoatAreaDict with km-sq in each mgmt zone
                        kiwiAreaDict[shpFile] = np.sum(data * self.kiwiCorrectionK)
                        stoatAreaDict[shpFile] = np.sum(data * self.stoatPercentArea)
                        # store it
                        kiwiSpatialDict[shpFile] = mask.copy()

        # put in a special key = 'ALL' that contains the extent mask
        # to make it easier when doing the stats
        kiwiSpatialDict['ALL'] = self.kiwiExtentMask
        kiwiAreaDict['ALL'] = np.sum(self.kiwiExtentMask * self.kiwiCorrectionK)
        stoatAreaDict['ALL'] = np.sum(self.stoatExtentMask * self.stoatPercentArea)

        return(kiwiSpatialDict, kiwiAreaDict, stoatAreaDict)



    def pickleSelf(self, fname):
        fileobj = open(fname, 'wb')
        pickle.dump(self, fileobj, protocol=4) # so we get large file support
        fileobj.close()

    @staticmethod
    def unpickleFromFile(fname):
        fileobj = open(fname, 'rb')
        data = pickle.load(fileobj)
        fileobj.close()
        return data

@jit
def scaleKiwiStoatMask(rodentResol, kiwiResol, kiwiNcols, kiwiNrows,
        DEM, originalExtent, kiwiMaxElev, stoatMaxElev, rodentExtentForStoats):
    """
    Calc proportion of pixels at rodent resolution that are suitable for kiwi
    """
    # number of rodent cells in one row within a kiwi cell
    oldPixPerNewPix = int(kiwiResol / rodentResol)
    # total number of rodent cells in one kiwi cell
    ncells = (oldPixPerNewPix)**2.0
    # new array to populate
    kiwiCorrectionK = np.zeros((kiwiNrows, kiwiNcols))
    stoatPercentArea = np.zeros((kiwiNrows, kiwiNcols))            # use to calc density in mgmt
    stoatExtentMask = np.zeros((kiwiNrows, kiwiNcols), np.uint8)
    # loop thru kiwi raster to populate
    for kiwiY in range(kiwiNrows):
        for kiwiX in range(kiwiNcols):
            kiwiTotal = 0.0
            stoatTotal = 0.0
            oldx = kiwiX * oldPixPerNewPix
            oldy = kiwiY * oldPixPerNewPix
            for x in range(oldPixPerNewPix):
                for y in range(oldPixPerNewPix):
                    addY = oldy + y
                    addX = oldx + x
                    ## get scale for kiwi
                    if DEM[addY, addX] <= kiwiMaxElev:
                        kiwiTotal += originalExtent[addY, addX]
                    if DEM[addY, addX] <= stoatMaxElev:
                        stoatTotal += originalExtent[addY, addX]
                    ## get count of rodent cells for stoats
#                    if stoatExtentMask[kiwiY, kiwiX] == 1:
#                        continue
                    # if rodent present in stoat pixel then indicate
                    if rodentExtentForStoats[addY, addX]:
                        stoatExtentMask[kiwiY, kiwiX] = 1
            stoatPercentArea[kiwiY, kiwiX] = stoatTotal / ncells
            kiwiCorrectionK[kiwiY, kiwiX] = kiwiTotal / ncells
    return(kiwiCorrectionK, stoatExtentMask, stoatPercentArea)



@jit(nopython=True)
# def getKiwiKMap(kiwiNcols, kiwiNrows, kiwiCorrectionK,
#         kiwiExtentMask, kiwiSurv, kiwiSurvDecay, kiwiRecDecay, kiwiProd):
def getKiwiKMap(kiwiNcols, kiwiNrows, kiwiCorrectionK,
        kiwiExtentMask, kiwiSurv, kiwiSurvDDcoef, kiwiRecDDcoef, kiwiProd, kiwiTheta):
    """
    ## MAKE KIWI EQUILIBRIUM POP DENSITY BY PIXEL FOR SENSITIVITY TEST
    """
    kiwi_KMap = np.zeros((kiwiNrows, kiwiNcols))
    n0 = 10.0
    # loop thru kiwi raster to populate
    for row in range(kiwiNrows):
        for col in range(kiwiNcols):
            if ~kiwiExtentMask[row, col]:
                continue
            N = n0
            prp = kiwiCorrectionK[row, col]
            for i in range(15):
                surv_i= kiwiSurv * np.exp(-((N/(kiwiSurvDDcoef*prp))**kiwiTheta))
                NStar = N * surv_i
                recRate = kiwiProd *np.exp(-(N/(kiwiRecDDcoef*prp))**kiwiTheta)
                N = (1 + recRate) * NStar
            kiwi_KMap[row, col] = N
    return(kiwi_KMap)
