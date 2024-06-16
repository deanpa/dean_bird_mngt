
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
from numba import njit

COMPRESSED_HFA = ['COMPRESSED=YES']

sr = osr.SpatialReference()
sr.ImportFromEPSG(2193) # always nztm?
NZTM_WKT = sr.ExportToWkt()

class FormatError(Exception):
    "Was unable to read the format"

class PreyData(object):
    """
    Holds the data used for the model run. This is calculated from 
    the parameters. Use the pickleSelf() method to save this to a pickle.
    """
    def __init__(self, params):
        self.params = params
        
        #self.seasAdjRes = np.genfromtxt(self.params.seasAdjResFile, delimiter=",",names=True)
        (self.seasAdj, self.mastSeasAdj, self.crashSeasAdj) = self.readInSeasResourceArrays()
          

        self.rodentGeoTrans, self.rodentNcols, self.rodentNrows, originalExtent = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[0]))
        self.rodentExtentMask = originalExtent > 0

        # do the other extentMasks (stoats and preys)
        self.stoatGeoTrans, self.stoatNcols, self.stoatNrows, stoatMask = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[1]))

        # DEM at rodent resolution resample with 'Average'
        self.DEM = self.resampleRasterRodent(self.params.DEM, self.params.resolutions[0],
                    GDTMethod = gdal.GDT_Float32, GRAMethod = gdal.GRA_Average)

        self.kClasses = self.resampleRasterRodent(self.params.kClasses, 
            self.params.resolutions[0], GDTMethod = gdal.GDT_UInt16, 
                    GRAMethod = gdal.GRA_NearestNeighbour)

        ## GET EXTENT TO READ IN PREY HABITAT MAP
        ogrDatasetTmp = ogr.Open(self.params.extentShp)
        ogrLayerTmp = ogrDatasetTmp.GetLayer()
        x0, x1, y0, y1 = ogrLayerTmp.GetExtent()
        fullExtTmp = [x0, x1, y0, y1]

        print('fullExtTmp', fullExtTmp)
        del ogrDatasetTmp, ogrLayerTmp
        ## RASTERISE PREY HABITAT AT PREY RESOLUTION AND FULL EXTENT
        self.preyGeoTrans, self.preyNcols, self.preyNrows, preyMask = (
            self.rasterizeShape(self.params.preyHabitatShp, self.params.resolutions[2],
                extent = fullExtTmp))
        # extentMasks for prey
        preyGeoTrans, preyNcols, preyNrows, tmpPreyMask = (
            self.rasterizeShape(self.params.extentShp, self.params.resolutions[2]))
        ## CLIP PREY MASK TO FULL EXTENT MASK

        print('sum preyMask', np.sum(preyMask), 'sum tmpMask', np.sum(tmpPreyMask))


        preyMask = preyMask & tmpPreyMask

        print('AFTER sum preyMask', np.sum(preyMask), 'sum tmpMask', np.sum(tmpPreyMask))

        del (preyGeoTrans, preyNcols, preyNrows, tmpPreyMask)

        ## FINEST RESOLUTION AMONG ALL SPECIES
        self.finestResol = np.min(self.params.resolutions)

        ## RASTERISE PREY HABITAT AT FINEST RESOLUTION FOR K CORRECTION
        (self.preyKCorrGeoTrans, self.preyKCorrNcols, self.preyKCorrNrows, 
            preyKCorrMask) = self.rasterizeShape(self.params.preyHabitatShp, 
            self.finestResol, extent = fullExtTmp)


        tmpGeoTrans, tmpNcols, tmpNrows, tmpExtent = (
            self.rasterizeShape(self.params.extentShp, self.finestResol))

        ## CLIP PREY K CORRECTION RASTER TO FULL EXTENT AT FINE RESOL
        preyKCorrMask = preyKCorrMask & tmpExtent
        del (tmpGeoTrans, tmpNcols, tmpNrows, tmpExtent)

        # add the self.rodentMaxAltitude to the self.rodentExtentMask and kClasses
        self.rodentExtentMask[self.DEM > self.params.rodentMaxAltitude] = 0
        self.kClasses[self.DEM > self.params.rodentMaxAltitude] = 0

        ### AREAS TRAPPED IN RECENT TIMES -- at rodent resol
        self.islands = self.resampleRasterRodent(self.params.islands, 
                    self.params.resolutions[0],
                    GDTMethod = gdal.GDT_UInt16, GRAMethod = gdal.GRA_NearestNeighbour)

        self.islands[~self.rodentExtentMask] = 0

        # Use to get stoat areas per zone -- still at rodent resolution 
        self.rodentExtentForStoats = self.rodentExtentMask.copy()


        kMapDS = gdal.Open(self.params.kClasses)
        ### MASK FOR CELLS THAT CAN MAST????
        rodentHabLU = rat.readColumn(kMapDS, "RodentHab") == 1
        rodentHabMask = rodentHabLU[self.kClasses]


        # remove non-habitat from rodentExtentMask -- at rodent resolution
        self.rodentExtentMask[~rodentHabMask] = 0

        # get preyCorrectionK to scale pixels near water or high elevation
        self.preyCorrectionK = scalePreyMask(self.finestResol, self.params.resolutions[2], 
            self.preyNcols, self.preyNrows, self.DEM, preyKCorrMask, self.params.preyMaxAltitude)
        ## MAKE PREY EXTENT MASK
        self.preyExtentMask = self.preyCorrectionK > 0.0


## DELETE BELOW IF ABOVE WORKS

#        # get preyCorrectionK to scale pixels near water or high elevation
#        self.preyCorrectionK = scalePreyMask(self.params.resolutions[0], 
#            self.params.resolutions[2], self.preyNcols, self.preyNrows,
#            self.DEM, originalExtent, self.params.preyMaxAltitude)
#########################
        
        ## GET REACTIVE CONTROL MONTHS FROM PARAMETERS
        self.getReactCtrlMth()

        (self.stoatExtentMask, self.stoatPercentArea) = (
            scaleStoatMask(self.params.resolutions[0], 
            self.params.resolutions[1], self.stoatNcols, self.stoatNrows,
            self.DEM, originalExtent, self.params.stoatMaxAltitude, 
            self.rodentExtentForStoats))
        self.stoatExtentMask = self.stoatExtentMask > 0


        self.rodentPercentArea = np.where(self.rodentExtentMask, 1.0, 0)

        self.preyKDummy = np.where(self.preyExtentMask > 0, 1.0, 0)

        maxPreyFec = self.params.preyPropBreedpa[1]*self.params.preyFec[3,1]
        ## PREY K MAP FOR SENSITIVITY TEST
        self.preyKMap = getPreyKMap(self.preyNcols, self.preyNrows,  
            self.preyCorrectionK, self.preyExtentMask, self.params.preySurv, 
            self.params.preySurvDDcoef, self.params.preyRecDDcoef,
            maxPreyFec, self.params.preyTheta)

#        kMapDS = gdal.Open(self.params.kClasses)
        ### MASK FOR CELLS THAT CAN MAST????
        self.mastingLU = rat.readColumn(kMapDS, "Masts") > 0
        ### 'CC' is the carrying capacity at per ha level.
        # self.paramRodentCCLU = rat.readColumn(kMapDS, "Rodent_CC")
        # self.paramRodentMastCCLU = rat.readColumn(kMapDS, "Rodent_MastCC")
        # self.paramRodentCrashCCLU = rat.readColumn(kMapDS, "Rodent_CrashCC")


#        print('self.mastingLU', self.mastingLU)
        #       self.paramRodentCCLU,'rodentMastCCLU', self.paramRodentMastCCLU,
        #       'rodentCrashCCLU', self.paramRodentCrashCCLU)
        # print('self.kClasses',self.kClasses)

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
        ds = driver.Create('PreyCorrectTemp.img', self.preyNcols, self.preyNrows, 
                1, gdal.GDT_Float32)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.preyGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.preyCorrectionK)
        del ds

        driver = gdal.GetDriverByName('GTiff')
        ds = driver.Create('demTemp.tif', self.rodentNcols, self.rodentNrows, 
                1, gdal.GDT_Float32)
        ds.SetProjection(NZTM_WKT)
        ds.SetGeoTransform(self.rodentGeoTrans)
        band = ds.GetRasterBand(1)
        band.WriteArray(self.DEM)
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

        (self.preySpatialDictByMgmt, self.preyAreaDictByMgmt) = (
            self.readAndResampleControlForPrey())

        (self.stoatSpatialDictByMgmt, self.stoatAreaDictByMgmt) = (
            self.readAndResampleControlForStoats())

        ## READ IN LEAD POINT DATA FROM DIRECTORY IF NOT NONE, MAKE 2D ARRAY OF PROBABILITIES
        self.readInLeadPointData()
        print('Do Lead Poisoning: ', self.doLeadPoisoning)
        if self.doLeadPoisoning:
            fullPreyExt = np.array(fullExtTmp)

            makeProbLeadDeath(self.pLeadDeath3D, self.preyExtentMask, self.preAdultMonthPLead, 
                self.adultMonthPLead, self.params.preySigma, self.nLeadPts, self.leadPoints, 
                fullPreyExt, self.params.resolutions[2], self.preyNcols, self.preyNrows)
        
            driver = gdal.GetDriverByName('KEA')
            ds = driver.Create('pLeadDiePreAdult.kea', self.preyNcols, self.preyNrows, 
                    1, gdal.GDT_Float32)
            ds.SetProjection(NZTM_WKT)
            ds.SetGeoTransform(self.preyGeoTrans)
            band = ds.GetRasterBand(1)
            band.WriteArray(self.pLeadDeath3D[0])
            del ds
            
            driver = gdal.GetDriverByName('KEA')
            ds = driver.Create('pLeadDieAdult.kea', self.preyNcols, self.preyNrows, 
                    1, gdal.GDT_Float32)
            ds.SetProjection(NZTM_WKT)
            ds.SetGeoTransform(self.preyGeoTrans)
            band = ds.GetRasterBand(1)
            band.WriteArray(self.pLeadDeath3D[1])
            del ds
        
            driver = gdal.GetDriverByName('KEA')
            ds = driver.Create('annualPLeadDiePreAdult.kea', self.preyNcols, 
                self.preyNrows, 1, gdal.GDT_Float32)
            ds.SetProjection(NZTM_WKT)
            ds.SetGeoTransform(self.preyGeoTrans)
            band = ds.GetRasterBand(1)
            band.WriteArray(self.pLeadDeath3D[2])
            del ds





    def readInLeadPointData(self):
        """
        ## READ IN LEAD POINT DATA FROM DIRECTORY IF NOT NONE
        """
        if self.params.leadPointData is not None:
            self.doLeadPoisoning = True        
            leadPtsTmp = np.genfromtxt(self.params.leadPointData, delimiter=",", names=True,
                dtype=['f8', 'f8'])
#                dtype=['S10', 'f8', 'f8'])
            xx = leadPtsTmp['xcoord']
            yy = leadPtsTmp['ycoord']
            self.nLeadPts = len(xx)
            self.leadPoints = np.hstack([np.expand_dims(xx, 1),np.expand_dims(yy, 1)])
            self.preAdultMonthPLead = 1.0 - np.exp(np.log(1.0 - self.params.pLeadMax['preAdult']) / 12.0)
            self.adultMonthPLead = 1.0 - np.exp(np.log(1.0 - self.params.pLeadMax['adult']) / 12.0)
#            self.monthlyLeadProb = {'preAdult' : preAdultMonthPLead, 'adult' : adultMonthPLead}
            shpPreyKMap = np.shape(self.preyKMap)
            self.pLeadDeath3D = np.zeros((3, shpPreyKMap[0], shpPreyKMap[1]), dtype = np.float32)
        else:
            self.doLeadPoisoning = False
            self.pLeadDeath3D = None


    def getReactCtrlMth(self):
        """
        ## GET REACTIVE CONTROL MONTHS
        """
        sumReactMth = self.params.reactiveAssessMth + self.params.reactiveCtrlDelay
        ## IF GREATER THAN MONTH 11 (12 MONTHS IN YEAR)
        if sumReactMth > 11:
            self.reactiveCtrlMth = sumReactMth - 12
            self.jumpYearCtrl = True
        else:
            self.reactiveCtrlMth = sumReactMth
            self.jumpYearCtrl = False
        #self.prescrptCtrlMth = self.reactiveCtrlMth  now set in params and allowed to be different

    def readInSeasResourceArrays(self):
        """
        Reads in resource file returning arrays for each of non mast, mast, crash
        The arrays are vegtype(0-5) * mths (0-11) shape
         
        """
        mcd={}
        with open(self.params.seasAdjResFile, 'r', newline='') as csvf:
            content = csv.reader(csvf, delimiter=",")
            next(content) #skip header
            for row in content:
                if not row[0] in mcd:
                    mcd[row[0]]=[]
                for value in row[2:]:
                    mcd[row[0]].append(value);  
        #convert values to float          
        res = {}
        for key, value in mcd.items():
            res[key] = [float(item) for item in value]
        
        #extract values to array: rows=vegtype, cols=mth
        seasAdj = np.array(res['nonmast'])
        seasAdj = np.reshape(seasAdj,(-1,12))    
        mastSeasAdj = np.array(res['mast'])
        mastSeasAdj = np.reshape(mastSeasAdj,(-1,12))    
        crashSeasAdj = np.array(res['postmast'])
        crashSeasAdj = np.reshape(crashSeasAdj,(-1,12))    
        
        print('preProcessing seasAdj', seasAdj[:,0])
        
        return(seasAdj, mastSeasAdj, crashSeasAdj)

    
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
        # # changed to prey since they have the larger resolution???? but then no data for edges
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

        self.controlAreaBool = []
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
                    shpFile, startYear, revisit, controlIndicator, = row
                    startYear = int(startYear)
                    revisit = int(revisit)
                    controlIndicator = eval(controlIndicator)
                    ## BOOLEAN ARRAY INDICATING WHETHER CONTROL SHOULD OCCUR IN AREA
                    self.controlAreaBool.append(controlIndicator)
                    
                    shpFile = os.path.join(self.params.inputDataPath, shpFile)

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
                    controlList.append((mask, startYear, revisit, controlIndicator, shpFile))
        ## IF SUMMARISE RESULTS OVER FULL EXTENT ADD TO 
        if self.params.summariseFullExtent:
            # put in a special key = 'ALL' that contains the extent mask
            # to make it easier when doing the stats
            rodentAreaDict['ALL'] = np.sum(self.rodentExtentMask)
            pixelsInControlAndBeechMask['ALL'] = 0
            controlList.append((self.rodentExtentMask.copy(), 100, -1, 300, 'ALL'))
#        print('rodentAreaDict', rodentAreaDict)
        return(controlList, rodentAreaDict, pixelsInControlAndBeechMask, controlAndBeechMask)

    def readAndResampleControlForPrey(self):
        """
        Similar to readAndResampleControlForRodents() above, but 
        resamples to Prey resolution and returns a dictionary
        keyed on the original shape file name so that analysis
        can be performed on each management area.
        """
        preySpatialDict = {}  # keyed on file name
        preyAreaDict = {}
        # work out the extent to use
        x0, y1 = gdal.ApplyGeoTransform(self.preyGeoTrans, 0, 0)
        x1, y0 = gdal.ApplyGeoTransform(self.preyGeoTrans, self.preyNcols, self.preyNrows)
        extent = [x0, x1, y0, y1]

        firstRow = True
        with open(self.params.controlFile, newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if firstRow:
                    firstRow = False
                    continue

                if len(row) == 4:
                    shpFile, startYear, revisit, controlIndicator = row
        
                    shpFile = os.path.join(self.params.inputDataPath, shpFile)

                    if shpFile not in preySpatialDict:
                        # we haven't come accross this one before
                        if not os.path.exists(shpFile):
                            raise IOError("Cannot find file %s" % shpFile)

                        try:
                            geoTrans, ncols, nrows, data = self.rasterizeShape(shpFile, 
                                    self.params.resolutions[2], extent=extent)
                        except FormatError:
                            # must be a raster input
                            data = self.resampleRaster(self.preyGeoTrans, self.preyNcols,
                                        self.preyNrows, shpFile, self.params.resolutions[2],
                                        GDTMethod = gdal.GDT_UInt16,
                                        GRAMethod = gdal.GRA_NearestNeighbour)
                        mask = (data == 1)
                        ### POPULATE prey and stoatAreaDict with km-sq in each mgmt zone
                        preyAreaDict[shpFile] = np.sum(data * self.preyCorrectionK)
                        # store it
                        preySpatialDict[shpFile] = mask.copy()

        ## IF SUMMARISE RESULTS OVER FULL EXTENT ADD TO 
        if self.params.summariseFullExtent:
            # put in a special key = 'ALL' that contains the extent mask
            # to make it easier when doing the stats
            preySpatialDict['ALL'] = self.preyExtentMask
            preyAreaDict['ALL'] = np.sum(self.preyExtentMask * self.preyCorrectionK)
        return(preySpatialDict, preyAreaDict)

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
                    shpFile, startYear, revisit, controlIndicator = row
        
                    shpFile = os.path.join(self.params.inputDataPath, shpFile)

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
        ## IF SUMMARISE RESULTS OVER FULL EXTENT ADD TO 
        if self.params.summariseFullExtent:
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

@njit
def scalePreyMask(finestResol, preyResol, preyNcols, preyNrows,
        DEM, preyKCorrMask, preyMaxElev):
    """
    Calc proportion of pixels at finest resolution that are suitable for prey
    """
    # number of rodent cells in one row within a prey cell
    oldPixPerNewPix = int(preyResol / finestResol)
    # total number of rodent cells in one prey cell
    ncells = (oldPixPerNewPix)**2.0
    # new array to populate
    preyCorrectionK = np.zeros((preyNrows, preyNcols))
    # loop thru prey raster to populate
    for preyY in range(preyNrows):
        for preyX in range(preyNcols):
            preyTotal = 0.0
            oldx = preyX * oldPixPerNewPix
            oldy = preyY * oldPixPerNewPix
            for x in range(oldPixPerNewPix):
                for y in range(oldPixPerNewPix):
                    addY = oldy + y
                    addX = oldx + x
                    ## get scale for prey
                    if DEM[addY, addX] <= preyMaxElev:
                        preyTotal += preyKCorrMask[addY, addX]
            preyCorrectionK[preyY, preyX] = preyTotal / ncells
    return(preyCorrectionK)


def scaleStoatMask(rodentResol, stoatResol, stoatNcols, stoatNrows,
        DEM, originalExtent, stoatMaxElev, rodentExtentForStoats):
    """
    Calc proportion of pixels at rodent resolution that are suitable for stoats
    """
    # number of rodent cells in one row within a prey cell
    oldPixPerNewPix = int(stoatResol / rodentResol)
    # total number of rodent cells in one prey cell
    ncells = (oldPixPerNewPix)**2.0
    # new array to populate
    stoatPercentArea = np.zeros((stoatNrows, stoatNcols))            # use to calc density in mgmt
    stoatExtentMask = np.zeros((stoatNrows, stoatNcols), np.uint8)
    # loop thru prey raster to populate
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
#                    if stoatExtentMask[preyY, preyX] == 1:
#                        continue
                    # if rodent present in stoat pixel then indicate
                    if rodentExtentForStoats[addY, addX]:
                        stoatExtentMask[stoatY, stoatX] = 1
            stoatPercentArea[stoatY, stoatX] = stoatTotal / ncells
    return(stoatExtentMask, stoatPercentArea)


@njit
def getPreyKMap(preyNcols, preyNrows, preyCorrectionK,
        preyExtentMask, preySurv, preySurvDDcoef, preyRecDDcoef, preyMaxFec, preyTheta):
    """
    ## MAKE PREY EQUILIBRIUM POP DENSITY BY PIXEL FOR SENSITIVITY TEST
    """
    prey_KMap = np.zeros((preyNrows, preyNcols))
    n0 = 10.0
    # loop thru prey raster to populate
    for row in range(preyNrows):
        for col in range(preyNcols):
            if ~preyExtentMask[row, col]:
                continue
            N = n0
            prp = preyCorrectionK[row, col]
            for i in range(15):
                surv_i= preySurv[3,1] * np.exp(-((N/(preySurvDDcoef*prp))**preyTheta))
                NStar = N * surv_i                
                recRate = preyMaxFec * np.exp(-(N/(preyRecDDcoef*prp))**preyTheta)
                N = (1 + recRate) * NStar
            prey_KMap[row, col] = N
    return(prey_KMap)


@njit
def makeProbLeadDeath(pLeadDeath3D, preyExtentMask, preAdultMonthPLead, adultMonthPLead, 
        preySigma, nLeadPts, leadPoints, fullExtTmp, resol, preyNcols, preyNrows):
    """
    ## MAKE 3D ARRAY AT PREY RESOL WITH PROB OF ENC, INT, DEATH
    ## FIRST LEVEL IS FOR PRE-ADULTS AND SECOND IS FOR ADULTS
    """
    ulx = fullExtTmp[0]
    uly = fullExtTmp[3]
    # loop thru prey raster to populate
    for row in range(preyNrows):
        for col in range(preyNcols):
            if ~preyExtentMask[row, col]:
                continue
            pNotDiePreAdult = 1.0
            pNotDieAdult = 1.0
            xx = ulx + col*resol
            yy = uly - row*resol
            ## LOOP THRU LEAD POINTS TO GET DISTANCE
            for pt in range(nLeadPts):
                dist = np.sqrt((xx - leadPoints[pt, 0])**2 + (yy - leadPoints[pt, 1])**2) 
                pDiePreAdult = preAdultMonthPLead * np.exp(-(dist**2) / 2.0 / (preySigma**2))
                pNotDiePreAdult = pNotDiePreAdult * (1.0 - pDiePreAdult)

                pDieAdult = adultMonthPLead * np.exp(-(dist**2) / 2.0 / (preySigma**2))
                pNotDieAdult = pNotDieAdult * (1.0 - pDieAdult)

            pLeadDeath3D[0, row, col] = 1.0 - pNotDiePreAdult
            pLeadDeath3D[1, row, col] = 1.0 - pNotDieAdult
            pLeadDeath3D[2, row, col] = 1.0 - pNotDiePreAdult**12


