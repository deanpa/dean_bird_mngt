#!/usr/bin/env python

import os
import shutil
import tempfile
import subprocess
import csv
import numpy as np
from osgeo import gdal
from osgeo import gdalconst
from scipy.stats.mstats import mquantiles
import pylab
from wand.image import Image
from wand.drawing import Drawing
from kiwimodelMthlyTS import calcresults
from kiwimodelMthlyTS import preProcessing


#########################################
#
#### ON NESI, HAVE TO: module load FFmpeg/3.4.5-GCCcore-7.4.0 
#
#########################################

x =1

# resize to this percent for each rodent image
# stoat and kiwi resizes are calculated from this
RODENT_RESIZE_PERCENT = 30 

# these may need tweaking
RODENT_DENSITY_RANGE = (0.0, 80.0)
STOAT_DENSITY_RANGE = (0.0, 8.0)
KIWI_DENSITY_RANGE = (0.0, 18.0)

# used for creating colour images of the densities
COLOUR_TABLE = 'colourtable.npy'

# labels
TEXT_LOCATION = (25, 10)

# output 
#OUTPUT_MOVIE = 'results.wmv'

FRAMERATE = 0.5

COLOUR_RAMP_FRACTION = 0.1 # of the image width

## NUMBER OF YEARS OVER WHICH CALC KIWI ANN GROWTH RATE
annGrowthYears = 5



def processResults():
    ### Output data path to write graphs, imagees and movies
    outputDataPath = os.path.join(os.getenv('KIWIPROJDIR', default='.'), 
            'KiwiProjResults', 'modC_Results')
    ## movie file path name
    movieFName = os.path.join(outputDataPath, 'results.wmv')

    resultsDataPath = os.path.join(outputDataPath, 'results.pkl')
    results = calcresults.KiwiResults.unpickleFromFile(resultsDataPath)



    # if results isn't already a list it came from testcalc.py (rather than testcalmulti.py)
    # so turn it into a list
    if not isinstance(results, list):
        results = [results]

    # TODO: we only need the preProcessing at present to 
    # get the names of the control areas. Maybe they should be
    # in the results so we don't have to??
    preProcessDataPath = os.path.join(outputDataPath, 'preProcData.pkl')
    data = preProcessing.KiwiData.unpickleFromFile(preProcessDataPath)

    # first, do the plots
    doTimeSeriesPlots(results, sorted(data.kiwiSpatialDictByMgmt.keys()), outputDataPath)
    doMthlyTimeSeriesPlots(results, sorted(data.kiwiSpatialDictByMgmt.keys()), outputDataPath)
###    # first, do the plots
###    doKiwiPlots(results, sorted(data.kiwiControlDictByMgmt.keys()), outputDataPath)

    ## MAKE TABLE OF CONTROL COUNTS; MEANS AND 95% CI
    controlCountTable(results, outputDataPath)




    # then the movie
    makeMovie(results, movieFName, outputDataPath)


    #########
    ###############################################################
    # OPTIONAL: MAKE 3-D TIFF OF DENSITIES OVER TIME
    tempTifName = os.path.join(outputDataPath, 'stoatDensity.tif')
    gdt_type = gdalconst.GDT_Float32

#    rasterS = results[0].popAllYears_3D['stoatDensity']

#    print('raster', np.shape(rasterS))
    writeTif(results, tempTifName, gdt_type, preProcessing.NZTM_WKT, 
        data.stoatGeoTrans)
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
    kiwiPNG = os.path.join(tempDir, 'kiwiDensity.png')

    for yearn in range(len(params.years)):
        yearName = params.years[yearn]

        # Note: Current directory
        thisFramePNG = 'frame_%02d.png' % yearn

        mastingMask = popMovie['Mastt'][yearn]
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

        kiwiDensity = popMovie['kiwiDensity'][yearn]
#        print('kiwiDensity', kiwiDensity.min(), kiwiDensity.max())

        # fudge
        kiwiResizePercent = ((rodentDensity.shape[0] / kiwiDensity.shape[0]) 
                    * RODENT_RESIZE_PERCENT)
        makeColourMapPNG(tempDir, kiwiDensity, kiwiPNG, 'Kiwi', 
                kiwiResizePercent, KIWI_DENSITY_RANGE)

        # make the frame with all the inputs
        frameDataPath = os.path.join(outputDataPath, thisFramePNG)
        # subprocess.check_call(['montage', mastingPNG, controlPNG, 
        #         rodentPNG, stoatPNG, kiwiPNG,'-geometry', 
        #         '+2+2', frameDataPath])
        filesList =  [mastingPNG, controlPNG, rodentPNG, stoatPNG, kiwiPNG]
        newfig = Image()
        for fname in filesList:        
            with Image(filename=fname) as img:
                newfig.sequence.append(img)
                newfig.smush(False, 0)
            newfig.save(filename=frameDataPath)        

    # now make the movie
    subprocess.check_call(['ffmpeg', '-framerate', str(FRAMERATE), '-i', 
        os.path.join(outputDataPath, 'frame_%02d.png'),
        '-vcodec', 'mpeg4', '-q:v', '1', '-y', 
        '-loglevel', 'error', movieFName])

    # tidy up
    shutil.rmtree(tempDir)

def doKiwiPlots(results, controlKeys, outputDataPath):

    nAreas, nYears = results[0].kiwiPropKMap_2D.shape
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
                val = results[iter].kiwiPropKMap_2D[i, year]
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
    nAreas, nYears = results[0].kiwiDensity_2D.shape
    mngtYears = nYears - burnin
    annGrowYearsStop = burnin + annGrowthYears - 1
    annGrowYears = np.min([annGrowthYears, mngtYears])
    nIterations = len(results)
    print('nIterations', nIterations)
    popChangeArray = np.empty((nAreas, ), dtype=[('Area', 'U32'), ('TotalMean', float), 
            ('TotalLow_CI', float), ('TotalHigh_CI', float),
            ('AnnualMean', float), ('AnnualLow_CI', float), 
            ('AnnualHigh_CI', float), 
            ('ProbIncrease', float)])
    # go through the management zones and plot
    for i, key in enumerate(controlKeys):
        # this ended up easier than stacking 
        kiwiMeansAllYears = []
        kiwiQuantsAllYears = []
        stoatMeansAllYears = []
        stoatQuantsAllYears = []
        rodentMeansAllYears = []
        rodentQuantsAllYears = []
        for year in range(nYears):
            kiwiDataThisYear = []
            stoatDataThisYear = []
            rodentDataThisYear = []
            for iter in range(nIterations):
                ## POPULATE 1-D LISTS
                val = results[iter].kiwiDensity_2D[i, year]
                kiwiDataThisYear.append(val)
                val = results[iter].stoatDensity_2D[i, year]
                stoatDataThisYear.append(val)
                val = results[iter].rodentDensity_2D[i, year] / nHectInRodent
                rodentDataThisYear.append(val)
            kiwiMeansAllYears.append(np.mean(kiwiDataThisYear))
            kiwiQuantsAllYears.append(mquantiles(kiwiDataThisYear, prob=[0.025, 0.975]))
            stoatMeansAllYears.append(np.mean(stoatDataThisYear))
            stoatQuantsAllYears.append(mquantiles(stoatDataThisYear, prob=[0.025, 0.975]))
            rodentMeansAllYears.append(np.mean(rodentDataThisYear))
            rodentQuantsAllYears.append(mquantiles(rodentDataThisYear, prob=[0.025, 0.975]))
            ## FOR CALCULATION OF POPULATION CHANGE FOR TABLE
            if year == burnin:
                kiwiChange = np.array(kiwiDataThisYear).copy()

            if year == annGrowYearsStop:
                ## ANNUAL PROPORTION GROWTH            
                annualPercGrow = np.log(kiwiDataThisYear / kiwiChange) / annGrowYears

            if year == (nYears - 1):
                deltaPop = kiwiDataThisYear - kiwiChange
                probIncrease = np.sum(deltaPop > 0.0) / nIterations
                kiwiChange = deltaPop / kiwiChange
    
        popChangeArray[i][0] = os.path.basename(key)
        popChangeArray[i][1] = np.mean(kiwiChange)
        quants = mquantiles(kiwiChange, prob=[0.025, 0.975])
        popChangeArray[i][2] = quants[0]
        popChangeArray[i][3] = quants[1]


        ## ADD IN ANNUAL PERCENT GROWTH    
        popChangeArray[i][4] = np.mean(annualPercGrow)
        quants = mquantiles(annualPercGrow, prob=[0.025, 0.975])
        popChangeArray[i][5] = quants[0]
        popChangeArray[i][6] = quants[1]

        popChangeArray[i][7] = probIncrease
        doAreaPlot(i, key, kiwiMeansAllYears, kiwiQuantsAllYears, 
                stoatMeansAllYears, stoatQuantsAllYears, 
                rodentMeansAllYears, rodentQuantsAllYears, outputDataPath, burnin)
#    print('popChange', popChangeArray)
    tableFilePathName = os.path.join(outputDataPath, 'percentKiwiChange.csv')
    np.savetxt(tableFilePathName, popChangeArray, 
        fmt=['%s', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'],
        comments = '', delimiter=',', 
        header='Area, TotalMean, TotalLow_CI, TotalHigh_CI, AnnualMean, AnnualLow_CI, AnnualHigh_CI, ProbIncrease')

def doMthlyTimeSeriesPlots(results, controlKeys, outputDataPath):   
    rodentResol = results[0].params.resolutions[0]
    # number of hectares in a rodent pixel
    nHectInRodent = (rodentResol / 100.0)**2
    nYrs = results[0].kiwiDensity_2D.shape[1]
    nAreas, nMths = results[0].kiwiDensity_2D_mth.shape
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
    MeanKiwi = np.empty((0),dtype=float)
    KiwiLow_CI = np.empty((0),dtype=float)
    KiwiHigh_CI = np.empty((0),dtype=float)
    
    # go through the management zones and plot
    for i, key in enumerate(controlKeys):
        # this ended up easier than stacking 
        kiwiMeansAllMths = []
        kiwiQuantsAllMths = []
        stoatMeansAllMths = []
        stoatQuantsAllMths = []
        rodentMeansAllMths = []
        rodentQuantsAllMths = []
        for mth in range(nMths):
            kiwiDataThisMth = []
            stoatDataThisMth = []
            rodentDataThisMth = []
            for iter in range(nIterations):
                ## POPULATE 1-D LISTS
                val = results[iter].kiwiDensity_2D_mth[i, mth]
                kiwiDataThisMth.append(val)
                val = results[iter].stoatDensity_2D_mth[i, mth]
                stoatDataThisMth.append(val)
                val = results[iter].rodentDensity_2D_mth[i, mth] / nHectInRodent
                rodentDataThisMth.append(val)
            kiwiMeansAllMths.append(np.mean(kiwiDataThisMth))
            kiwiQuantsAllMths.append(mquantiles(kiwiDataThisMth, prob=[0.025, 0.975]))
            stoatMeansAllMths.append(np.mean(stoatDataThisMth))
            stoatQuantsAllMths.append(mquantiles(stoatDataThisMth, prob=[0.025, 0.975]))
            rodentMeansAllMths.append(np.mean(rodentDataThisMth))
            rodentQuantsAllMths.append(mquantiles(rodentDataThisMth, prob=[0.025, 0.975]))
    
        doMthlyAreaPlot(i, key, kiwiMeansAllMths, kiwiQuantsAllMths, 
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
        MeanKiwi=np.append(MeanKiwi,kiwiMeansAllMths,axis=0)
        KiwiLow_CI=np.append(KiwiLow_CI,np.array(kiwiQuantsAllMths)[...,0],axis=0)
        KiwiHigh_CI=np.append(KiwiHigh_CI,np.array(kiwiQuantsAllMths)[...,1],axis=0)

    tableFilePathName = os.path.join(outputDataPath, 'monthlyDensities.csv')
    with open (tableFilePathName, 'w', newline='') as f:
        headers = ['Area', 'Year','Rodents_Mean','Rodents_LowCI', 'Rodents_HighCI',
                   'Stoats_Mean','Stoats_LowCI', 'Stoats_HighCI','Kiwi_Mean',
                   'Kiwi_LowCI', 'Kiwi_HighCI',]
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        writer.writerows(zip(area, decimYear, MeanRodents, RodentsLow_CI, 
                             RodentsHigh_CI, MeanStoats, StoatsLow_CI,
                             StoatsHigh_CI, MeanKiwi, KiwiLow_CI, KiwiHigh_CI))


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

    # # write text
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
    #             '+0+0', '-fill', 'yellow', '-pointsize', '10',
    #             '-draw', 'text %d, %d "%s"' % (TEXT_LOCATION[0] + rampWidth, TEXT_LOCATION[1], title),
    #             fname])
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


def doAreaPlot(i, key, kiwiMeansAllYears, kiwiQuantsAllYears, 
                stoatMeansAllYears, stoatQuantsAllYears, 
                rodentMeansAllYears, rodentQuantsAllYears, outputDataPath, burnin):
    kiwiQuantsAllYears = np.array(kiwiQuantsAllYears)
    stoatQuantsAllYears = np.array(stoatQuantsAllYears)
    rodentQuantsAllYears = np.array(rodentQuantsAllYears)

    pylab.figure(figsize = (11,9))
    pylab.subplot(3,1,3)
    pylab.plot(kiwiMeansAllYears, linewidth=4, color = 'k')
    pylab.plot(kiwiQuantsAllYears[..., 0], linewidth=1, color = 'k')
    pylab.plot(kiwiQuantsAllYears[..., 1], linewidth=1, color = 'k')
    YMAX = np.max(kiwiQuantsAllYears[..., 1])
    pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.xlabel('Years', fontsize = 12)
    pylab.ylabel('Kiwi density \n $(ind. * km^{-2})$', fontsize = 12,
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

def doMthlyAreaPlot(i, key, kiwiMeansAllMths, kiwiQuantsAllMths, 
                stoatMeansAllMths, stoatQuantsAllMths, 
                rodentMeansAllMths, rodentQuantsAllMths, outputDataPath):
    kiwiQuantsAllMths = np.array(kiwiQuantsAllMths)
    stoatQuantsAllMths = np.array(stoatQuantsAllMths)
    rodentQuantsAllMths = np.array(rodentQuantsAllMths)

    pylab.figure(figsize = (11,9))
    pylab.subplot(3,1,3)
    pylab.plot(kiwiMeansAllMths, linewidth=4, color = 'k')
    pylab.plot(kiwiQuantsAllMths[..., 0], linewidth=1, color = 'k')
    pylab.plot(kiwiQuantsAllMths[..., 1], linewidth=1, color = 'k')
    #YMAX = np.max(kiwiQuantsAllMths[..., 1])
    #pylab.vlines(x = (burnin), ymin = 0, ymax = YMAX, linestyles = 'dashed', colors='k')
    pylab.xlabel('Months', fontsize = 12)
    pylab.ylabel('Kiwi density \n $(ind. * km^{-2})$', fontsize = 12,
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
    pylab.ylabel('Kiwi proportion of kiwi_KMap')
    title = os.path.basename(key)
    pylab.title(title)

    
    outname = 'prop_plot_%s.png' % i

    filePathName = os.path.join(outputDataPath, outname)

    pylab.savefig(filePathName)
    pylab.cla()
    


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




    
if __name__ == '__main__':
    processResults()
