
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

class KeaData(object):
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

        # do the other extentMasks (stoats and keas)
        self.stoatGeoTrans, self.stoatNcols, self.stoatNrows, stoatMask = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[1]))

        # DEM at rodent resolution resample with 'Average'
        self.DEM = self.resampleRasterRodent(self.params.DEM, self.params.resolutions[0],
                    GDTMethod = gdal.GDT_Float32, GRAMethod = gdal.GRA_Average)

        self.kClasses = self.resampleRasterRodent(self.params.kClasses, self.params.resolutions[0],
                    GDTMethod = gdal.GDT_UInt16, GRAMethod = gdal.GRA_NearestNeighbour)

        self.keaGeoTrans, self.keaNcols, self.keaNrows, keaMask = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[2]))


        print('keanrows', self.keaNcols, self.keaNrows, type(self.keaNrows))


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

        # get keaCorrectionK to scale pixels near water or high elevation
        self.keaCorrectionK = scaleKeaMask(self.params.resolutions[0], 
            self.params.resolutions[2], self.keaNcols, self.keaNrows,
            self.DEM, originalExtent, self.params.keaMaxAltitude)
        self.keaExtentMask = self.keaCorrectionK > 0.0

        (self.stoatExtentMask, self.stoatPercentArea) = (
            scaleStoatMask(self.params.resolutions[0], 
            self.params.resolutions[1], self.stoatNcols, self.stoatNrows,
            self.DEM, originalExtent, self.params.stoatMaxAltitude, 
            self.rodentExtentForStoats))
        self.stoatExtentMask = self.stoatExtentMask > 0


        self.rodentPercentArea = np.where(self.rodentExtentMask, 1.0, 0)
        self.keaKDummy = np.where(self.keaExtentMask, 1.0, 0)

        ## KEA K MAP FOR SENSITIVITY TEST
        # self.keaKMap = getKeaKMap(self.keaNcols, self.keaNrows,  
        #     self.keaCorrectionK, self.keaExtentMask, self.params.keaSurv, 
        #     self.params.keaSurvDecay, self.params.keaRecDecay,
        #     self.params.keaProd)
        self.keaKMap = getKeaKMap(self.keaNcols, self.keaNrows,  
            self.keaCorrectionK, self.keaExtentMask, self.params.keaSurv, 
            self.params.keaSurvDDcoef, self.params.keaRecDDcoef,
            self.params.keaProd, self.params.keaTheta)

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
        ds = driver.Create('KeaCorrectTemp.img', self.keaNcols, self.keaNrows, 
                1, gdal.GDT_Float32)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.keaGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.keaCorrectionK)
        del ds

        driver = gdal.GetDriverByName('HFA')
        ds = driver.Create('keaExtentMaskTemp.img', self.stoatNcols, self.stoatNrows, 
                1, gdal.GDT_Byte)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.stoatGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.keaExtentMask)
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


        (self.keaSpatialDictByMgmt, self.keaAreaDictByMgmt) = self.readAndResampleControlForKeas()

        (self.stoatSpatialDictByMgmt, self.stoatAreaDictByMgmt) = self.readAndResampleControlForStoats()

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
        # # changed to kea since they have the larger resolution???? but then no data for edges
        # xmin = np.floor(x0 / self.params.resolutions[2]) * self.params.resolutions[2]
        # xmax = np.ceil(x1 / self.params.resolutions[2]) * self.params.resolutions[2]
        # ymin = np.floor(y0 / self.params.resolutions[2]) * self.params.resolutions[2]
        # ymax = np.ceil(y1 / self.params.resolutions[2]) * self.params.resolutions[2]

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

    def resampleRaster(geoTrans, nCol, nRow, infile, resolution, GDTMethod, GRAMethod):
        """
        Resample the input file to the given resolution and bounds
        """
        driver = gdal.GetDriverByName('HFA')
        gdalDataset = driver.Create('resraster.img', nCol, nRow, 1,
                GDTMethod, COMPRESSED_HFA)
        gdalDataset.SetGeoTransform(geoTrans)

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

    def readAndResampleControlForKeas(self):
        """
        Similar to readAndResampleControlForRodents() above, but 
        resamples to Kea resolution and returns a dictionary
        keyed on the original shape file name so that analysis
        can be performed on each management area.
        """
        keaSpatialDict = {}  # keyed on file name
        keaAreaDict = {}
        # work out the extent to use
        x0, y1 = gdal.ApplyGeoTransform(self.keaGeoTrans, 0, 0)
        x1, y0 = gdal.ApplyGeoTransform(self.keaGeoTrans, self.keaNcols, self.keaNrows)
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

                    if shpFile not in keaSpatialDict:
                        # we haven't come accross this one before
                        if not os.path.exists(shpFile):
                            raise IOError("Cannot find file %s" % shpFile)

                        try:
                            geoTrans, ncols, nrows, data = self.rasterizeShape(shpFile, 
                                    self.params.resolutions[2], extent=extent)
                        except FormatError:
                            # must be a raster input
                            data = self.resampleRaster(self.keaGeoTrans, self.keaNcols,
                                        self.keaNrows, shpFile, self.params.resolutions[2],
                                        GDTMethod = gdal.GDT_UInt16,
                                        GRAMethod = gdal.GRA_NearestNeighbour)
                        mask = (data == 1)
                        ### POPULATE kea and stoatAreaDict with km-sq in each mgmt zone
                        keaAreaDict[shpFile] = np.sum(data * self.keaCorrectionK)
                        # store it
                        keaSpatialDict[shpFile] = mask.copy()

        # put in a special key = 'ALL' that contains the extent mask
        # to make it easier when doing the stats
        keaSpatialDict['ALL'] = self.keaExtentMask
        keaAreaDict['ALL'] = np.sum(self.keaExtentMask * self.keaCorrectionK)

        return(keaSpatialDict, keaAreaDict)

    def readAndResampleControlForStoats(self):
        """
        Similar to readAndResampleControlForRodents() above, but 
        resamples to Stoat resolution and returns a dictionary
        keyed on the original shape file name so that analysis
        can be performed on each management area.
        """
        stoatSpatialDict = {}  # keyed on file name
        stoatAreaDict = {}
        # work out the extent to use
        x0, y1 = gdal.ApplyGeoTransform(self.stoatGeoTrans, 0, 0)
        x1, y0 = gdal.ApplyGeoTransform(self.stoatGeoTrans, self.stoatNcols, self.stoatNrows)
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

                    if shpFile not in stoatSpatialDict:
                        # we haven't come accross this one before
                        if not os.path.exists(shpFile):
                            raise IOError("Cannot find file %s" % shpFile)

                        try:
                            geoTrans, ncols, nrows, data = self.rasterizeShape(shpFile, 
                                    self.params.resolutions[1], extent=extent)
                        except FormatError:
                            # must be a raster input
                            data = self.resampleRaster(self.stoatGeoTrans, self.stoatNcols,
                                        self.stoatNrows, shpFile, self.params.resolutions[1],
                                        GDTMethod = gdal.GDT_UInt16,
                                        GRAMethod = gdal.GRA_NearestNeighbour)
                        mask = (data == 1)
                        ### POPULATE stoatAreaDict with km-sq in each mgmt zone
                        stoatAreaDict[shpFile] = np.sum(data * self.stoatPercentArea)
                        # store it
                        stoatSpatialDict[shpFile] = mask.copy()

        # put in a special key = 'ALL' that contains the extent mask
        # to make it easier when doing the stats
        stoatSpatialDict['ALL'] = self.stoatExtentMask
        stoatAreaDict['ALL'] = np.sum(self.stoatExtentMask * self.stoatPercentArea)

        return(stoatSpatialDict, stoatAreaDict)


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
def scaleKeaMask(rodentResol, keaResol, keaNcols, keaNrows,
        DEM, originalExtent, keaMaxElev):
    """
    Calc proportion of pixels at rodent resolution that are suitable for kea
    """
    # number of rodent cells in one row within a kea cell
    oldPixPerNewPix = int(keaResol / rodentResol)
    # total number of rodent cells in one kea cell
    ncells = (oldPixPerNewPix)**2.0
    # new array to populate
    keaCorrectionK = np.zeros((keaNrows, keaNcols))
    # loop thru kea raster to populate
    for keaY in range(keaNrows):
        for keaX in range(keaNcols):
            keaTotal = 0.0
            oldx = keaX * oldPixPerNewPix
            oldy = keaY * oldPixPerNewPix
            for x in range(oldPixPerNewPix):
                for y in range(oldPixPerNewPix):
                    addY = oldy + y
                    addX = oldx + x
                    ## get scale for kea
                    if DEM[addY, addX] <= keaMaxElev:
                        keaTotal += originalExtent[addY, addX]
            keaCorrectionK[keaY, keaX] = keaTotal / ncells
    return(keaCorrectionK)

def scaleStoatMask(rodentResol, stoatResol, stoatNcols, stoatNrows,
        DEM, originalExtent, stoatMaxElev, rodentExtentForStoats):
    """
    Calc proportion of pixels at rodent resolution that are suitable for stoats
    """
    # number of rodent cells in one row within a kea cell
    oldPixPerNewPix = int(stoatResol / rodentResol)
    # total number of rodent cells in one kea cell
    ncells = (oldPixPerNewPix)**2.0
    # new array to populate
    stoatPercentArea = np.zeros((stoatNrows, stoatNcols))            # use to calc density in mgmt
    stoatExtentMask = np.zeros((stoatNrows, stoatNcols), np.uint8)
    # loop thru kea raster to populate
    for stoatY in range(stoatNrows):
        for stoatX in range(stoatNcols):
            stoatTotal = 0.0
            oldx = stoatX * oldPixPerNewPix
            oldy = stoatY * oldPixPerNewPix
            for x in range(oldPixPerNewPix):
                for y in range(oldPixPerNewPix):
                    addY = oldy + y
                    addX = oldx + x
                    ## get scale for stoats
                    if DEM[addY, addX] <= stoatMaxElev:
                        stoatTotal += originalExtent[addY, addX]
                    ## get count of rodent cells for stoats
#                    if stoatExtentMask[keaY, keaX] == 1:
#                        continue
                    # if rodent present in stoat pixel then indicate
                    if rodentExtentForStoats[addY, addX]:
                        stoatExtentMask[stoatY, stoatX] = 1
            stoatPercentArea[stoatY, stoatX] = stoatTotal / ncells
    return(stoatExtentMask, stoatPercentArea)


@jit(nopython=True)
# def getKeaKMap(keaNcols, keaNrows, keaCorrectionK,
#         keaExtentMask, keaSurv, keaSurvDecay, keaRecDecay, keaProd):
def getKeaKMap(keaNcols, keaNrows, keaCorrectionK,
        keaExtentMask, keaSurv, keaSurvDDcoef, keaRecDDcoef, keaProd, keaTheta):
    """
    ## MAKE KEA EQUILIBRIUM POP DENSITY BY PIXEL FOR SENSITIVITY TEST
    """
    kea_KMap = np.zeros((keaNrows, keaNcols))
    n0 = 10.0
    # loop thru kea raster to populate
    for row in range(keaNrows):
        for col in range(keaNcols):
            if ~keaExtentMask[row, col]:
                continue
            N = n0
            prp = keaCorrectionK[row, col]
            for i in range(15):
                surv_i= keaSurv[3] * np.exp(-((N/(keaSurvDDcoef*prp))**keaTheta))
                NStar = N * surv_i
                recRate = keaProd *np.exp(-(N/(keaRecDDcoef*prp))**keaTheta)
                N = (1 + recRate) * NStar
            kea_KMap[row, col] = N
    return(kea_KMap)
