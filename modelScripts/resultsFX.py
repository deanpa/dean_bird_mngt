#!/usr/bin/env python

import os
import shutil
import tempfile
import subprocess
import csv
import numpy as np
from osgeo import gdal
from osgeo import gdalconst
from osgeo import osr
from osgeo import ogr
from scipy.stats.mstats import mquantiles
import pylab
#from wand.image import Image
#from wand.drawing import Drawing
from modelScripts import calcresults
#from modelScripts import preProcessing


COMPRESSED_HFA = ['COMPRESSED=YES']
sr = osr.SpatialReference()
sr.ImportFromEPSG(2193) # always nztm?
NZTM_WKT = sr.ExportToWkt()

#########################################
#
#### ON NESI, HAVE TO: module load FFmpeg/3.4.5-GCCcore-7.4.0 
#
#########################################

x =1

# resize to this percent for each rodent image
# stoat and prey resizes are calculated from this
RODENT_RESIZE_PERCENT = 175  #30

# these may need tweaking
RODENT_DENSITY_RANGE = (0.0, 80.0)
STOAT_DENSITY_RANGE = (0.0, 8.0)
PREY_DENSITY_RANGE = (0.0, 18.0)

# used for creating colour images of the densities
COLOUR_TABLE = 'colourtable.npy'

# labels
TEXT_LOCATION = (25, 10)

# output 
#OUTPUT_MOVIE = 'results.wmv'

FRAMERATE = 0.5

COLOUR_RAMP_FRACTION = 0.1 # of the image width

## NUMBER OF YEARS OVER WHICH CALC PREY ANN GROWTH RATE
annGrowthYears = 5 #5      ## CHANGE TO 5


def writeTmpArray(params, data, results):
    # 'MastT', 'preyDensity', 'rodentDensity'
    driver = gdal.GetDriverByName('HFA')
    ds = driver.Create('mastTmp.img', data.rodentNcols, data.rodentNrows, 
         1, gdal.GDT_Byte)
    ds.SetProjection(NZTM_WKT)
    ds.SetGeoTransform(data.rodentGeoTrans)
    band = ds.GetRasterBand(1)
    band.WriteArray(results[0].popAllYears_3D['MastT'][-1])
    del ds



def processResults(params, data, results):
    ### Output data path to write graphs, imagees and movies
#    outputDataPath = os.path.join(os.getenv('PREYPROJDIR', default='.'), 
#            'modelResults', 'modC_Results')
    ## movie file path name
    movieFName = os.path.join(params.outputDataPath, 'results.wmv')

#    resultsDataPath = os.path.join(outputDataPath, 'results.pkl')
#    results = calcresults.PreyResults.unpickleFromFile(resultsDataPath)



    # if results isn't already a list it came from testcalc.py (rather than testcalmulti.py)
    # so turn it into a list
    if not isinstance(results, list):
        results = [results]

    # TODO: we only need the preProcessing at present to 
    # get the names of the control areas. Maybe they should be
    # in the results so we don't have to??
###    preProcessDataPath = os.path.join(outputDataPath, 'preProcData.pkl')
###    data = preProcessing.PreyData.unpickleFromFile(preProcessDataPath)

    # first, do the plots
    doTimeSeriesPlots(results, sorted(data.preySpatialDictByMgmt.keys()), params.outputDataPath)
    doMthlyTimeSeriesPlots(results, sorted(data.preySpatialDictByMgmt.keys()), params.outputDataPath)


###    # first, do the plots
###    doPreyPlots(results, sorted(data.preyControlDictByMgmt.keys()), params.outputDataPath)

    ## MAKE TABLE OF CONTROL COUNTS; MEANS AND 95% CI
    controlCountTable(results, params.outputDataPath)




    # then the movie
    makeMovie(results, movieFName, params.outputDataPath)


    #########
    ###############################################################
    # OPTIONAL: MAKE 3-D TIFF OF DENSITIES OVER TIME
###    tempTifName = os.path.join(params.outputDataPath, 'stoatDensity.tif')
###    gdt_type = gdalconst.GDT_Float32

#    rasterS = results[0].popAllYears_3D['stoatDensity']

#    print('raster', np.shape(rasterS))
###    writeTif(results, tempTifName, gdt_type, preProcessing.NZTM_WKT, 
###        data.stoatGeoTrans)
    ##############################################################
    #########






def makeMovie(results, movieFName, outputDataPath):
    # create temp dir to work in
    # TODO: do we need to be able to specify the dir this 
    # happens in?

    tempDir = tempfile.mkdtemp()

    # get the iteration for the movie info
    params = results[0].params
    popMovie = results[0].popAllYears_3D

    mastingPNG = os.path.join(tempDir, 'mastingMask.png')
    controlPNG = os.path.join(tempDir, 'controlMask.png')
    rodentPNG = os.path.join(tempDir, 'rodentDensity.png')
    stoatPNG = os.path.join(tempDir, 'stoatDensity.png')
    preyPNG = os.path.join(tempDir, 'preyDensity.png')

    for yearn in range(len(params.years)):
        yearName = params.years[yearn]

        # Note: Current directory
        thisFramePNG = 'frame_%02d.png' % yearn

        print('thisFramePNG', thisFramePNG)

        mastingMask = popMovie['MastT'][yearn]
        makeMaskPNG(tempDir, mastingMask, mastingPNG, 
                'Masting Year %d' % yearName, 
                RODENT_RESIZE_PERCENT)

        controlMask = popMovie['ControlT'][yearn]
        makeMaskPNG(tempDir, controlMask, controlPNG, 'Control',
                RODENT_RESIZE_PERCENT)


        rodentDensity = popMovie['rodentDensity'][yearn]
        makeColourMapPNG(tempDir, rodentDensity, rodentPNG, 'Rats', 
                RODENT_RESIZE_PERCENT, RODENT_DENSITY_RANGE)

        stoatDensity = popMovie['stoatDensity'][yearn]
        # fudge
        stoatResizePercent = ((rodentDensity.shape[0] / stoatDensity.shape[0]) 
                    * RODENT_RESIZE_PERCENT)
        makeColourMapPNG(tempDir, stoatDensity, stoatPNG, 'Stoats', 
            stoatResizePercent, STOAT_DENSITY_RANGE)

        preyDensity = popMovie['preyDensity'][yearn]

        # fudge
        preyResizePercent = ((rodentDensity.shape[0] / preyDensity.shape[0]) 
                    * RODENT_RESIZE_PERCENT)
        makeColourMapPNG(tempDir, preyDensity, preyPNG, 'Prey', 
                preyResizePercent, PREY_DENSITY_RANGE)

        # make the frame with all the inputs
        frameDataPath = os.path.join(outputDataPath, thisFramePNG)
        # subprocess.check_call(['montage', mastingPNG, controlPNG, 
        #           rodentPNG, stoatPNG, preyPNG,'-geometry', 
        #           '+2+2', frameDataPath])
        filesList =  [mastingPNG, controlPNG, rodentPNG, stoatPNG, preyPNG]
        newfig = Image()
        for fname in filesList:        
            with Image(filename=fname) as img:
                newfig.sequence.append(img)
                newfig.smush(False, 0)
            newfig.save(filename=frameDataPath)        

        
        #now make the movie
#    subprocess.check_call(['ffmpeg', '-framerate', str(FRAMERATE), '-i', 
#        os.path.join(params.outputDataPath, 'frame_%02d.png'),
#        '-vcodec', 'mpeg4', '-q:v', '1', '-y', 
#        '-loglevel', 'error', movieFName])

    # tidy up
    shutil.rmtree(tempDir)

def doPreyPlots(results, controlKeys, outputDataPath):

    nAreas, nYears = results[0].preyPropKMap_2D.shape
    nIterations = len(results)
    print('nIterations', nIterations)
    # go through the management zones and plot
    for i, key in enumerate(controlKeys):

        # this ended up easier than stacking 
        meansAllYears = []
        quantsAllYears = []
        for year in range(nYears):
            dataThisYear = []
            for iter in range(nIterations):
                val = results[iter].preyPropKMap_2D[i, year]
                dataThisYear.append(val)

            mean = np.mean(dataThisYear)
            meansAllYears.append(mean)

            quants = mquantiles(dataThisYear, prob=[0.025, 0.975])
            quantsAllYears.append(quants)

        doPlot(i, key, meansAllYears, quantsAllYears, outputDataPath)


def doTimeSeriesPlots(results, controlKeys, outputDataPath):
    burnin = results[0].params.burnin    
    rodentResol = results[0].params.resolutions[0]
    # number of hectares in a rodent pixel
    nHectInRodent = (rodentResol / 100.0)**2
    nAreas, nYears = results[0].preyDensity_2D.shape
    mngtYears = nYears - burnin
    annGrowYearsStop = burnin + annGrowthYears - 1
    annGrowYears = np.min([annGrowthYears, mngtYears])
    nIterations = len(results)
    print('nIterations', nIterations, nYears, mngtYears)
    popChangeArray = np.empty((nAreas, ), dtype=[('Area', 'U32'), ('TotalMean', float), 
            ('TotalLow_CI', float), ('TotalHigh_CI', float),
            ('AnnualMean', float), ('AnnualLow_CI', float), 
            ('AnnualHigh_CI', float), 
            ('ProbIncrease', float)])
    # go through the management zones and plot
    for i, key in enumerate(controlKeys):
        # this ended up easier than stacking 
        preyMeansAllYears = []
        preyQuantsAllYears = []
        stoatMeansAllYears = []
        stoatQuantsAllYears = []
        rodentMeansAllYears = []
        rodentQuantsAllYears = []
        for year in range(nYears):
            preyDataThisYear = []
            stoatDataThisYear = []
            rodentDataThisYear = []
            for iter in range(nIterations):
                ## POPULATE 1-D LISTS
                val = results[iter].preyDensity_2D[i, year]
                preyDataThisYear.append(val)
                val = results[iter].stoatDensity_2D[i, year]
                stoatDataThisYear.append(val)
                val = results[iter].rodentDensity_2D[i, year] / nHectInRodent
                rodentDataThisYear.append(val)
            preyMeansAllYears.append(np.mean(preyDataThisYear))
            preyQuantsAllYears.append(mquantiles(preyDataThisYear, prob=[0.025, 0.975]))
            stoatMeansAllYears.append(np.mean(stoatDataThisYear))
            stoatQuantsAllYears.append(mquantiles(stoatDataThisYear, prob=[0.025, 0.975]))
            rodentMeansAllYears.append(np.mean(rodentDataThisYear))
            rodentQuantsAllYears.append(mquantiles(rodentDataThisYear, prob=[0.025, 0.975]))
            ## FOR CALCULATION OF POPULATION CHANGE FOR TABLE
            if year == burnin:
                preyChange = np.array(preyDataThisYear).copy()

            if year == annGrowYearsStop:
                ## ANNUAL PROPORTION GROWTH            
                annualPercGrow = np.log(preyDataThisYear / preyChange) / annGrowYears

            if year == (nYears - 1):
                deltaPop = preyDataThisYear - preyChange
                probIncrease = np.sum(deltaPop > 0.0) / nIterations
                preyChange = deltaPop / preyChange
    
        popChangeArray[i][0] = os.path.basename(key)
        popChangeArray[i][1] = np.mean(preyChange)
        quants = mquantiles(preyChange, prob=[0.025, 0.975])
        popChangeArray[i][2] = quants[0]
        popChangeArray[i][3] = quants[1]


        ## ADD IN ANNUAL PERCENT GROWTH    
        popChangeArray[i][4] = np.mean(annualPercGrow)
        quants = mquantiles(annualPercGrow, prob=[0.025, 0.975])
        popChangeArray[i][5] = quants[0]
        popChangeArray[i][6] = quants[1]

        popChangeArray[i][7] = probIncrease
        doAreaPlot(i, key, preyMeansAllYears, preyQuantsAllYears, 
                stoatMeansAllYears, stoatQuantsAllYears, 
                rodentMeansAllYears, rodentQuantsAllYears, outputDataPath, burnin)
#    print('popChange', popChangeArray)
    tableFilePathName = os.path.join(outputDataPath, 'percentPreyChange.csv')
    np.savetxt(tableFilePathName, popChangeArray, 
        fmt=['%s', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'],
        comments = '', delimiter=',', 
        header='Area, TotalMean, TotalLow_CI, TotalHigh_CI, AnnualMean, AnnualLow_CI, AnnualHigh_CI, ProbIncrease')

def doMthlyTimeSeriesPlots(results, controlKeys, outputDataPath):   
    rodentResol = results[0].params.resolutions[0]
    # number of hectares in a rodent pixel
    nHectInRodent = (rodentResol / 100.0)**2
    nYrs = results[0].preyDensity_2D.shape[1]
    nAreas, nMths = results[0].preyDensity_2D_mth.shape
    nIterations = len(results)
    decYrs=np.linspace(start=0, stop=(nMths/12),endpoint=False,num=nMths)
    area =[]
    decimYear = np.empty((0), dtype=float)
    MeanRodents = np.empty((0),dtype=float)
    RodentsLow_CI = np.empty((0),dtype=float)
    RodentsHigh_CI = np.empty((0),dtype=float)
    MeanStoats = np.empty((0),dtype=float)
    StoatsLow_CI =np.empty((0),dtype=float)
    StoatsHigh_CI = np.empty((0),dtype=float)
    MeanPrey = np.empty((0),dtype=float)
    PreyLow_CI = np.empty((0),dtype=float)
    PreyHigh_CI = np.empty((0),dtype=float)

  
    # go through the management zones and plot
    for i, key in enumerate(controlKeys):
        # this ended up easier than stacking 
        preyMeansAllMths = []
        preyQuantsAllMths = []
        stoatMeansAllMths = []
        stoatQuantsAllMths = []
        rodentMeansAllMths = []
        rodentQuantsAllMths = []
        for mth in range(nMths):
            preyDataThisMth = []
            stoatDataThisMth = []
            rodentDataThisMth = []
            for iter in range(nIterations):
                ## POPULATE 1-D LISTS
                val = results[iter].preyDensity_2D_mth[i, mth]
                preyDataThisMth.append(val)
                val = results[iter].stoatDensity_2D_mth[i, mth]
                stoatDataThisMth.append(val)
                val = results[iter].rodentDensity_2D_mth[i, mth] / nHectInRodent
                rodentDataThisMth.append(val)
            preyMeansAllMths.append(np.mean(preyDataThisMth))
            preyQuantsAllMths.append(mquantiles(preyDataThisMth, prob=[0.025, 0.975]))
            stoatMeansAllMths.append(np.mean(stoatDataThisMth))
            stoatQuantsAllMths.append(mquantiles(stoatDataThisMth, prob=[0.025, 0.975]))
            rodentMeansAllMths.append(np.mean(rodentDataThisMth))
            rodentQuantsAllMths.append(mquantiles(rodentDataThisMth, prob=[0.025, 0.975]))
    
        doMthlyAreaPlot(i, key, preyMeansAllMths, preyQuantsAllMths, 
                stoatMeansAllMths, stoatQuantsAllMths, 
                rodentMeansAllMths, rodentQuantsAllMths, outputDataPath)
        
        area[(i*nMths):(i*nMths+nMths-1)]=np.repeat(os.path.basename(key), nMths)
        decimYear=np.append(decimYear,decYrs,axis=0)
        MeanRodents=np.append(MeanRodents,rodentMeansAllMths,axis=0)
        RodentsLow_CI=np.append(RodentsLow_CI,np.array(rodentQuantsAllMths)[...,0],axis=0)
        RodentsHigh_CI=np.append(RodentsHigh_CI,np.array(rodentQuantsAllMths)[...,1],axis=0)
        MeanStoats=np.append(MeanStoats,stoatMeansAllMths,axis=0)
        StoatsLow_CI=np.append(StoatsLow_CI,np.array(stoatQuantsAllMths)[...,0],axis=0)
        StoatsHigh_CI=np.append(StoatsHigh_CI,np.array(stoatQuantsAllMths)[...,1],axis=0)
        MeanPrey=np.append(MeanPrey,preyMeansAllMths,axis=0)
        PreyLow_CI=np.append(PreyLow_CI,np.array(preyQuantsAllMths)[...,0],axis=0)
        PreyHigh_CI=np.append(PreyHigh_CI,np.array(preyQuantsAllMths)[...,1],axis=0)
    tableFilePathName = os.path.join(outputDataPath, 'monthlyDensities.csv')
    with open (tableFilePathName, 'w', newline='') as f:
        headers = ['Area', 'Year','Rodents_Mean','Rodents_LowCI', 'Rodents_HighCI',
                   'Stoats_Mean','Stoats_LowCI', 'Stoats_HighCI','Prey_Mean',
                   'Prey_LowCI', 'Prey_HighCI']
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        writer.writerows(zip(area, decimYear, MeanRodents, RodentsLow_CI, 
                             RodentsHigh_CI, MeanStoats, StoatsLow_CI,
                             StoatsHigh_CI, MeanPrey, PreyLow_CI, PreyHigh_CI))


def controlCountTable(results, outputDataPath):
    countStructured = np.empty((1,), dtype=[('Mean', float), ('Low_CI', float), 
                    ('High_CI', float)])
    nIterations = len(results)
    countsThisIter = []
    for iter in range(nIterations):
        val = results[iter].controlCount
        countsThisIter.append(val)
    meanCount = np.mean(countsThisIter)
    quantsCount = mquantiles(countsThisIter, prob=[0.025, 0.975])
    countStructured['Mean'] = meanCount
    countStructured['Low_CI'] = quantsCount[0]
    countStructured['High_CI'] = quantsCount[1]
    tableFilePathName = os.path.join(outputDataPath, 'controlCount.csv')
    np.savetxt(tableFilePathName, countStructured, fmt=['%.2f', '%.2f', '%.2f'],
                    comments = '', delimiter=',', header='Mean, Low_CI, High_CI')


def makeMaskPNG(tempDir, mask, fname, title, resizePercent):
    """
    Create a rgb PNG for the mask
    """
    nrows, ncols = mask.shape

    # write to a .img file as we can't directly 
    # create a .png file
    maskIMG = os.path.join(tempDir, 'mask.img')

    driver = gdal.GetDriverByName('HFA')
    ds = driver.Create(maskIMG, ncols, nrows, 3, gdal.GDT_Byte, 
        ['COMPRESSED=YES'])

    for n in range(3):
        band = ds.GetRasterBand(n+1)
        band.WriteArray(np.where(mask, 255, 0))

    del ds

    # now create .png
    subprocess.check_call(['gdal_translate', '-of', 'PNG', '-outsize', 
        '%d%%' % resizePercent, '%d%%' % resizePercent, maskIMG, fname])
    os.remove(maskIMG)

    # write text
    # subprocess.check_call(['mogrify', '-fill', 'yellow', '-pointsize', '10', 
    #     '-draw', 'text %d, %d "%s"' % (TEXT_LOCATION[0], TEXT_LOCATION[1], title), 
    #     fname])
    with Image(filename=fname) as img:
        with Drawing() as ctx:
            ctx.fill_color = 'yellow'   
            ctx.font_size = 10
            ctx.text(TEXT_LOCATION[0], TEXT_LOCATION[1], title)
            ctx(img)
        img.save(filename=fname)

def makeColourMapPNG(tempDir, density, fname, title, resizePercent, densityRange):
    """
    Using the colour map create a rgb PNG for the density array
    """
    nrows, ncols = density.shape

    # read the map
    colourTable = np.load(COLOUR_TABLE)

    # write to a .img file as we can't directly 
    # create a .png file with GDAL
    densityIMG = os.path.join(tempDir, 'density.img')
    densityPNG = os.path.join(tempDir, 'density.png')
    rampIMG = os.path.join(tempDir, 'ramp.img')
    rampPNG = os.path.join(tempDir, 'ramp.png')

    driver = gdal.GetDriverByName('HFA')
    ds = driver.Create(densityIMG, ncols, nrows, 3, gdal.GDT_Byte, 
        ['COMPRESSED=YES'])

    # rescale density so the range is 0-75
    density = (density / (densityRange[1] - densityRange[0])) * 75
    density = np.clip(density, 0, 75).astype(np.uint8)
    #print('density', fname, density.min(), density.max())

    for n in range(3):
        band = ds.GetRasterBand(n+1)
        data = colourTable[n][density]
        band.WriteArray(data)

    del ds

    # now create .png
    subprocess.check_call(['gdal_translate', '-of', 'PNG', '-outsize', 
        '%d%%' % resizePercent, '%d%%' % resizePercent, densityIMG, densityPNG])

    os.remove(densityIMG)

    # now create the colour ramp. Get size of output
    ds = gdal.Open(densityPNG)
    nrows, ncols = ds.RasterYSize, ds.RasterXSize
    del ds

    rampWidth = int(COLOUR_RAMP_FRACTION * ncols)
    ramp = np.linspace(0, 75, nrows).astype(np.uint8)
    # duplicate to required width and rotate
    ramp = np.vstack([ramp] * rampWidth)
    ramp = np.rot90(ramp)

    # write to .img
    ds = driver.Create(rampIMG, rampWidth, nrows, 3, gdal.GDT_Byte, 
        ['COMPRESSED=YES'])

    for n in range(3):
        band = ds.GetRasterBand(n+1)
        data = colourTable[n][ramp]
        band.WriteArray(data)

    del ds

    # convert ramp to .png
    subprocess.check_call(['gdal_translate', '-of', 'PNG', rampIMG, rampPNG])
    os.remove(rampIMG)

    # # write labels on ramp first as montage writes text on all input images
    # subprocess.check_call(['mogrify', '-fill', 'black', '-pointsize', '8',
    #     '-draw', 'text %d, %d "%.1f"' % (0, 8, densityRange[1]),
    #     '-draw', 'text %d, %d "%.1f"' % (0, int(nrows/2), (densityRange[1] - densityRange[0]) / 2),
    #     '-draw', 'text %d, %d "%.1f"' % (0, nrows-2, densityRange[0]),
    #     rampPNG])
    with Image(filename=rampPNG) as img:
        with Drawing() as ctx:
            ctx.fill_color = 'black'   
            ctx.font_size = 8
            ctx.text(0, 8, f'{densityRange[1]}')
            ctx.text(0, int(nrows/2), f'{(densityRange[1] - densityRange[0]) / 2}')
            ctx.text(0, nrows-2, f'{densityRange[0]}')
            ctx(img)
        img.save(filename=rampPNG)  

    # do mosaic and write text - thankfully the title is far enough over that 
    # it doesn't appear in the ramp image also
    # subprocess.check_call(['montage', rampPNG, densityPNG, '-geometry', 
    #     '+0+0', '-fill', 'yellow', '-pointsize', '10',
    #     '-draw', 'text %d, %d "%s"' % (TEXT_LOCATION[0] + rampWidth, TEXT_LOCATION[1], title), fname])
    fg = Image(filename = densityPNG)
    with Drawing() as ctx:
        ctx.fill_color = 'yellow'   
        ctx.font_size = 10
        ctx.text(TEXT_LOCATION[0] + rampWidth, TEXT_LOCATION[1], title)
        ctx(fg)
    fg.composite(image=Image(filename = rampPNG), left = 0, top = 0)    
    fg.save(filename = fname)
    os.remove(densityPNG)
    os.remove(rampPNG)

def doAreaPlot(i, key, preyMeansAllYears, preyQuantsAllYears, 
                stoatMeansAllYears, stoatQuantsAllYears, 
                rodentMeansAllYears, rodentQuantsAllYears, outputDataPath, burnin):
    preyQuantsAllYears = np.array(preyQuantsAllYears)
    stoatQuantsAllYears = np.array(stoatQuantsAllYears)
    rodentQuantsAllYears = np.array(rodentQuantsAllYears)

    pylab.figure(figsize = (11,9))
    pylab.subplot(3,1,3)
    pylab.plot(preyMeansAllYears, linewidth=4, color = 'k')
    pylab.plot(preyQuantsAllYears[..., 0], linewidth=1, color = 'k')
    pylab.plot(preyQuantsAllYears[..., 1], linewidth=1, color = 'k')
    YMAX = np.max(preyQuantsAllYears[..., 1])
    pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.xlabel('Years', fontsize = 12)
    pylab.ylabel('Prey density \n $(ind. * km^{-2})$', fontsize = 12,
        rotation = 'horizontal', ha = 'right')

    pylab.subplot(3,1,2)
    pylab.plot(stoatMeansAllYears, linewidth=4, color = 'k')
    pylab.plot(stoatQuantsAllYears[..., 0], linewidth=1, color = 'k')
    pylab.plot(stoatQuantsAllYears[..., 1], linewidth=1, color = 'k')
    YMAX = np.max(stoatQuantsAllYears[..., 1])
    pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.ylabel('Stoat density \n $(ind. * km^{-2})$', fontsize = 12,
        rotation = 'horizontal', ha = 'right')

    pylab.subplot(3,1,1)
    pylab.plot(rodentMeansAllYears, linewidth=4, color = 'k')
    pylab.plot(rodentQuantsAllYears[..., 0], linewidth=1, color = 'k')
    pylab.plot(rodentQuantsAllYears[..., 1], linewidth=1, color = 'k')
    YMAX = np.max(rodentQuantsAllYears[..., 1])
    pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.ylabel('Rat density \n $(ind. * ha^{-1})$', fontsize = 12,
        rotation = 'horizontal', ha = 'right')

    title = os.path.basename(key)
    pylab.title(title)
    pylab.tight_layout()
    
    outname = 'density_plot_%s.png' % i

    filePathName = os.path.join(outputDataPath, outname)

    pylab.savefig(filePathName)
    pylab.cla()

def doMthlyAreaPlot(i, key, preyMeansAllMths, preyQuantsAllMths, 
                stoatMeansAllMths, stoatQuantsAllMths, 
                rodentMeansAllMths, rodentQuantsAllMths, outputDataPath):
    preyQuantsAllMths = np.array(preyQuantsAllMths)
    stoatQuantsAllMths = np.array(stoatQuantsAllMths)
    rodentQuantsAllMths = np.array(rodentQuantsAllMths)

    pylab.figure(figsize = (11,9))
    pylab.subplot(3,1,3)
    pylab.plot(preyMeansAllMths, linewidth=4, color = 'k')
    pylab.plot(preyQuantsAllMths[..., 0], linewidth=1, color = 'k')
    pylab.plot(preyQuantsAllMths[..., 1], linewidth=1, color = 'k')
    #YMAX = np.max(preyQuantsAllMths[..., 1])
    #pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.xlabel('Months', fontsize = 12)
    pylab.ylabel('Prey density \n $(ind. * km^{-2})$', fontsize = 12,
        rotation = 'horizontal', ha = 'right')

    pylab.subplot(3,1,2)
    pylab.plot(stoatMeansAllMths, linewidth=4, color = 'k')
    pylab.plot(stoatQuantsAllMths[..., 0], linewidth=1, color = 'k')
    pylab.plot(stoatQuantsAllMths[..., 1], linewidth=1, color = 'k')
    #YMAX = np.max(stoatQuantsAllMths[..., 1])
    #pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.ylabel('Stoat density \n $(ind. * km^{-2})$', fontsize = 12,
        rotation = 'horizontal', ha = 'right')

    pylab.subplot(3,1,1)
    pylab.plot(rodentMeansAllMths, linewidth=4, color = 'k')
    pylab.plot(rodentQuantsAllMths[..., 0], linewidth=1, color = 'k')
    pylab.plot(rodentQuantsAllMths[..., 1], linewidth=1, color = 'k')
    #YMAX = np.max(rodentQuantsAllMths[..., 1])
    #pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.ylabel('Rat density \n $(ind. * ha^{-1})$', fontsize = 12,
        rotation = 'horizontal', ha = 'right')

    title = os.path.basename(key)
    pylab.title(title)
    pylab.tight_layout()
    
    outname = 'density_plot_mthly_%s.png' % i

    filePathName = os.path.join(outputDataPath, outname)

    pylab.savefig(filePathName)
    pylab.cla()

def doPlot(i, key, meansAllYears, quantsAllYears, outputDataPath):
    quantsAllYears = np.array(quantsAllYears)

    #print('quantsAllYears', quantsAllYears)
    pylab.plot(meansAllYears, linewidth=3)
    pylab.plot(quantsAllYears[..., 0], linewidth=1)
    pylab.plot(quantsAllYears[..., 1], linewidth=1)
    pylab.xlabel('Years')
    pylab.ylabel('Prey proportion of prey_KMap')
    title = os.path.basename(key)
    pylab.title(title)

    
    outname = 'prop_plot_%s.png' % i

    filePathName = os.path.join(outputDataPath, outname)

    pylab.savefig(filePathName)
    pylab.cla()
    

def writeMultiTif(results, data, params, arrNames = None, yearIndx = None):
    """
    write single or multi-band tifs to directory
    """
    ## DEFAULT RASTER NAMES TO WRITE TO TIF, IF NOT SPECIFIED
    if arrNames is None:
        arrNames = ['MastT', 'ControlT', 'rodentDensity', 'stoatDensity', 'preyDensity']
#    nNames = len(arrNames)
    ## DEFAULT YEARS TO WRITE TO TIF: LAST FIVE YEARS
    if yearIndx is None:
        nTotalLayers = len(results[0].popAllYears_3D['stoatDensity'])
        if nTotalLayers >= 5:
            yearIndx = [nTotalLayers - 5, nTotalLayers]
        else:
            yearIndx = [0, nTotalLayers]
    nLayers = yearIndx[1] - yearIndx[0]
    ## DICTIONARY OF DATA TYPES
    gdt_DType = {'MastT' : gdal.GDT_Byte, 'ControlT' : gdal.GDT_Byte, 
        'rodentDensity' : gdal.GDT_Float32, 'stoatDensity' : gdal.GDT_Float32, 
        'preyDensity' : gdal.GDT_Float32}
    geoTrans_Dict = {'MastT' : data.rodentGeoTrans, 'ControlT' : data.rodentGeoTrans, 
        'rodentDensity' : data.rodentGeoTrans, 'stoatDensity' : data.stoatGeoTrans, 
        'preyDensity' : data.preyGeoTrans}

    for nn in arrNames:
        FName = nn + '_3D.tif'
        outFNamePath = os.path.join(params.outputDataPath, FName)
        print('FName', outFNamePath)
        raster_nn = results[0].popAllYears_3D[nn][yearIndx[0]:yearIndx[1]]

        (nrows, ncols) = np.shape(raster_nn[0])
        ds = gdal.GetDriverByName('GTiff').Create(outFNamePath, ncols,
            nrows, nLayers, gdt_DType[nn],
            options=['TILED=YES', 'COMPRESS=LZW', 'INTERLEAVE=BAND', 'BIGTIFF=IF_SAFER'])
        ds.SetGeoTransform(geoTrans_Dict[nn])
        ds.SetProjection(NZTM_WKT)
        # loop thru years (layers in tif)
        for nLay in range(nLayers):
            band = ds.GetRasterBand(nLay+1)
            band.WriteArray(raster_nn[nLay])
    del ds  # Flush



def writeTif(results, tempTifName, gdt_type, wkt, match_geotrans):
    """
    write single or multi-band tifs to directory
    """
    raster = results[0].popAllYears_3D['stoatDensity']

    # if single band
#    if not isinstance(raster, list):
    if len(np.shape(raster)) == 2:
        (nrows, ncols) = np.shape(raster)
        ds = gdal.GetDriverByName('GTiff').Create(tempTifName, ncols,
        nrows, 1, gdt_type,
        options=['TILED=YES', 'COMPRESS=LZW', 'INTERLEAVE=BAND', 'BIGTIFF=IF_SAFER'])
        ds.SetGeoTransform(match_geotrans)
        ds.SetProjection(wkt)
        band = ds.GetRasterBand(1)
        band.WriteArray(raster)
    # if multi-layered array
    else:
        (nrows, ncols) = np.shape(raster[0])
        nlayers = len(raster)
        ds = gdal.GetDriverByName('GTiff').Create(tempTifName, ncols,
            nrows, nlayers, gdt_type,
            options=['TILED=YES', 'COMPRESS=LZW', 'INTERLEAVE=BAND', 'BIGTIFF=IF_SAFER'])
        ds.SetGeoTransform(match_geotrans)
        ds.SetProjection(wkt)
        # loop thru years (layers in tif)
        for n in range(nlayers):
            band = ds.GetRasterBand(n+1)
            band.WriteArray(raster[n])
    del ds  # Flush


