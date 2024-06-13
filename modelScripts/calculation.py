from osgeo import gdal
import numpy as np
from numba import njit
from scipy.stats.mstats import mquantiles
from scipy.special import logit

from .preProcessing import NZTM_WKT
from . import calcresults

rng = np.random.default_rng()

# for resampleRasterDown
RESAMPLE_SUM = 0
RESAMPLE_AVERAGE = 1

def round_up_to_odd(f):
    """
    Window sizes must be odd for calculations to work
    """
    f = int(np.ceil(f))
    return f + 1 if f % 2 == 0 else f

def inv_logit(x):
    """
    Inverse of the logit function
    """
    return np.exp(x) / (1 + np.exp(x))


def runModel(rawdata, params=None, loopIter=0):
    """
    Run the Prey model for a given iteration number.
    rawdata should be an instance of preProcessing.PreyData.
    params should be a params.PreyParams. If not given, the one
    in rawdata is used.
    loopIter is the iteration number that this model run is for.
    Returns an instance of calcresults.PreyResults.
    """
    if params is None:
        params = rawdata.params


#    print('in dict', 'All' in rawdata.rodentControlList)


#    print('loopIter in runModel', loopIter)

    # create results object.
    results = calcresults.PreyResults()
    # stash the params in case useful for plotting.
    results.params = params 

    nControlAreas = len(rawdata.preySpatialDictByMgmt)
    nYears = len(params.years)

    ### MAKE BURNIN AND KEEP MASK ARRAY
    totalYears = nYears + params.burnin
    keepMask = np.arange(totalYears) >= params.burnin

    # NOTE: data just for this iteration
    results.preyDensity_2D = np.zeros((nControlAreas, totalYears))
    results.stoatDensity_2D = np.zeros((nControlAreas, totalYears))
    results.rodentDensity_2D = np.zeros((nControlAreas, totalYears))
    results.preyDensity_2D_mth = np.zeros((nControlAreas, totalYears*12))
    results.stoatDensity_2D_mth = np.zeros((nControlAreas, totalYears*12))
    results.rodentDensity_2D_mth = np.zeros((nControlAreas, totalYears*12))


    ## COUNT NUMBER OF CONTROL OPERATIONS
    results.controlCount = 0

    stoatShp = np.shape(rawdata.stoatExtentMask)
    preyShp = np.shape(rawdata.preyExtentMask)
    # are we storing result.popAllYears_3D for this iteration?
    # if so, create the array.
    if loopIter == 0:
        rodentShp = np.shape(rawdata.rodentExtentMask)
        results.popAllYears_3D = {'MastT': np.zeros((nYears, rodentShp[0], rodentShp[1]), 
                                        dtype = bool),
                      'ControlT': np.zeros((nYears, rodentShp[0], rodentShp[1]), 
                                        dtype = bool),     
                      'rodentDensity': np.zeros((nYears, rodentShp[0], rodentShp[1]), 
                                        dtype = float),
                      'stoatDensity': np.zeros((nYears, stoatShp[0], stoatShp[1]), 
                                        dtype = float),
                      'preyDensity': np.zeros((nYears, preyShp[0], preyShp[1]), 
                                        dtype = float)}

    # For each year in params.years
    # 1. Masting Event (or not)
    # 2. Control (or not)
    # 3. Rodent Peak
    # 4. Rodent Dispersal
    # 5. Stoat Peak (using Rodent Population at control, and num dead toxic rodents)
    # 6. Stoat Dispersal
    # 7. Analysis of shape of Stoat Population Tail 
    # 8. Prey population decrease

    # Create our internal rasters
    rodentRasterShape = rawdata.DEM.shape

 
    ## KEEP THE SEEDS PER M^2 SAME FOR ALL RESOLUTIONS
    # rodent_kMap = rawdata.paramRodentCCLU[rawdata.kClasses] * 1
    # # modify rodent K when mast or crash occurs
    # rodent_kMapMast = rawdata.paramRodentMastCCLU[rawdata.kClasses] * 1
    # rodent_kMapCrash = rawdata.paramRodentCrashCCLU[rawdata.kClasses] * 1
    # print('k classes', rawdata.kClasses[1,1:10],'kMap', rodent_kMap[1,1:10], 
    #       'kMapMast', rodent_kMapMast[1,1:10],  'kMapCrash', rodent_kMapCrash[1,1:10] )
   
    # number of hectares in a rodent pixel; for getting stoat_K
    nHectInRodent = (params.resolutions[0] / 100.0)**2

    # update rawdata.rodentExtentMask so we don't include areas where kmap == 0
    #already done in preprocessing line#94
    #rawdata.rodentExtentMask = rawdata.rodentExtentMask & (rodent_kMap != 0)

    # update the rodent kMap to set k = 0 at elevation should have already been done in pre-proc
    # rodent_kMap[~rawdata.rodentExtentMask] = 0
    # rodent_kMapMast[~rawdata.rodentExtentMask] = 0
    # rodent_kMapCrash[~rawdata.rodentExtentMask] = 0

    # convert from window sizes in pixels to metres
    mastWindowSizePxls = round_up_to_odd(params.mastWindowSize / 
                params.resolutions[0])
    masthalfwinsize = int(mastWindowSizePxls / 2)

    rodentEmigrationWindowSizePxls = round_up_to_odd(
                params.emigrationWindowSize[0] / params.resolutions[0])

    stoatEmigrationWindowSizePxls = round_up_to_odd(
                params.emigrationWindowSize[1] / params.resolutions[1])

    preyEmigrationWindowSizePxls = round_up_to_odd(
                params.emigrationWindowSize[2] / params.resolutions[2])

    distArray = makeDistArray(mastWindowSizePxls, 
            params.mastRho, masthalfwinsize, params.resolutions[0])
    if np.isnan(distArray).any():
        raise ValueError("NaNs in distArray")

    beechMask = rawdata.mastingLU[rawdata.kClasses] & rawdata.rodentExtentForStoats
    

    #########
    ##  Generate initial populations at t-1
    ##  Assume no control, no mast and no dispersal
    #########
    # Eqn 8
    # random rodent population
    rodentPresence = rng.binomial(1, params.pRodentPres, (rodentRasterShape))
    rodentPresence[~rawdata.rodentExtentMask] = 0

    ## GET INITIAL RODENT DENSITY FOR EACH CELL AT RODENT RESOL (4 HA)
#    rodent_raster = initialPopWrapper(params.rodentProd,
#        params.rodentSurv, params.rodentSurvDecay, params.rodentRecDecay,
#        params.initialRodentN, params.rodentT,  
#        rodentPresence, rodentRasterShape, rodent_kMap)


    rodent_raster = rodentPresence * rng.poisson(params.initialRodentN,
            rodentRasterShape)


    rodent_raster[~rawdata.rodentExtentMask] = 0
    
    # create a time since control raster based on rodent extent and initialise
    # with months=rodentBouncePeriod so that doesn't have an effect yet
    timeSinceCtrl = np.zeros_like(rawdata.rodentExtentMask, dtype=(np.int16))
    timeSinceCtrl[rawdata.rodentExtentMask] = params.rodentBouncePeriod


    # get rodent population and toxic rodents at stoat resolution
    # sum rodents at rodent resolution, but raster is at stoat resol.

    ## RESCALE TO PER HA FOR STOAT CALCULATIONS
    haPixelRescale = (params.resolutions[1] / 100.0)**2

    islandPixRescale = (params.resolutions[1] / params.resolutions[0])**2 
    
    ## RODENTS DENSITY/HA AT STOAT RESOL
    rodent_raster_stoat = resampleRasterDown(rodent_raster, 
            rawdata.stoatExtentMask.shape, statMethod = RESAMPLE_SUM, 
            pixelRescale = haPixelRescale)

    # update for islands - stoatIslandPrp is a proportion
    stoatIslandPrp = resampleRasterDown(rawdata.islands, 
            rawdata.stoatExtentMask.shape, statMethod = RESAMPLE_SUM, 
            pixelRescale = islandPixRescale)
    stoatIslandMask = stoatIslandPrp > 0

    ## ARRAY OF RODENT DENSITY (PER HA) IN TRAPPED AREAS; SCALED FOR RODENT HABITAT
    stoatIslandKArray = params.islandK * (stoatIslandPrp[stoatIslandMask])

    ## ASSIGN NOTIONAL RODENT DENSITY TO ISLANDS
    rodent_raster_stoat[stoatIslandMask] = stoatIslandKArray
    # ## MAKE RODENT RASTER FOR t-lag
    # rodentRasterStoatT_lag = rodent_raster_stoat.copy()    
    
    #  Random stoat population and mask out 
    stoatExtentShape = rawdata.stoatExtentMask.shape
    # Eqn 22
    stoatPresence = rng.binomial(1, params.pStoatPres, stoatExtentShape)
    stoatPresence[~rawdata.stoatExtentMask] = 0

    ## GET INITIAL STOAT DENSITY FOR EACH CELL AT STOAT RESOL (1 KM)
    # Eqn 23
    stoat_raster = stoatPresence * rng.poisson(0.75 * params.initialStoatN,
                stoatShp) 
#    stoat_raster = stoatPresence * np.exp(np.random.normal(np.log(params.initialStoatN),
#        params.stoatPopSD), stoatShp).astype(int)

    ## REMOVE ALL NON-STOAT HABITAT
    stoat_raster[~rawdata.stoatExtentMask] = 0    

    ## Get initial prey population
    preyPresence = rng.binomial(1, params.pPreyPres, 
            rawdata.preyCorrectionK.shape)

    # Eqn 38 
    prey_raster = (rng.poisson(params.preyInitialMultiplier * params.preyK, 
                            rawdata.preyCorrectionK.shape) * preyPresence)
 

### TODO: CORRECT THIS. KEEP RESOL FOR STOATS AND PREY EQUAL FOR NOW

#    stoat_preyIslandMask = resampleRasterDown(stoatIslandMask, rawdata.preyExtentMask.shape, 
#                       statMethod = RESAMPLE_SUM, pixelRescale = 1)



    adjustPreyIsland = (prey_raster > 3) & stoatIslandMask  #stoat_preyIslandMask
    prey_raster[adjustPreyIsland] = 2
    prey_raster[~rawdata.preyExtentMask] = 0



    ## ADD AGE STRUCTURE TO EACH PIXEL IN PREY RASTER
    #split prey_raster into diff age classes
    ## EACH BAND IN 3D RASTER IS AN AGE GROUP (0-4)
    preyByAgeRasts = rng.multinomial(prey_raster, params.preyInitAgeStr)
    ## ADD THE TOTAL NUMBER IN EACH PIXEL AS A 6TH BAND IN 3D RASTER
    prey_raster = np.concatenate((preyByAgeRasts.transpose(2,0,1),prey_raster.reshape(1,
                        prey_raster.shape[0],prey_raster.shape[1])),axis=0)
    del preyByAgeRasts


    ## CORRECT PREY RECRUIT AND SURVIVAL PARAMETERS FOR HABITAT AREA
    preyRecDecay_1D = ((params.preyRecDDcoef* 
            rawdata.preyCorrectionK)[rawdata.preyExtentMask])
    preySurvDecay_1D = ((params.preySurvDDcoef * 
            rawdata.preyCorrectionK)[rawdata.preyExtentMask])


    # record when last control happens for each mask
    lastControlArray = np.zeros(len(rawdata.rodentControlList), dtype=int)

    # set initial no mast in year t-1
    mastT_1 = False
    oldMastingMask = np.zeros_like(beechMask, dtype=bool) 
    ## COLUMN 0 IS FOR CURRENT YEAR, COLUMN 1 IS FOR T-1
    propControlMaskMasting = np.zeros(nControlAreas, dtype=float) 
#    propControlMaskMasting = np.zeros((nControlAreas,2), dtype=float) 

    ## MONTHLY CONTROL SCHEDULE FOR FOLLOWING YEAR (T+1) IF ASSESS/CONTROL JUMPS YEAR LINE
    nextYearCtrlSched = np.full((nControlAreas), -1)

    ### SET INITIAL KEEP YEAR TO 0
    year = 0
    
    for year_all in range(totalYears):
        ### BOOLEAN INDICATOR OF A POST-BURNIN YEAR
        keepYear = keepMask[year_all]

        print('year_all', year_all, 'keepYear', keepYear, 'year', year)


        #reactiveControlMask = None # nothing by default

        mastingMask  = np.zeros_like(beechMask, dtype=bool)
        preyMastingMask = np.zeros_like(prey_raster[0,0:,0:], dtype=float)

        ## EMPTY OBJECT FOR STORING CONTROL MASK FOR ALL YEARS FOR RESULTS IN LOOP 0
        controlThisYear = None

        ##create an monthly control schedule for this year
        mthlyCtrlSched = nextYearCtrlSched
        ## RESET CONTROL SCHEDULE FOR NEXT YEAR (T+1) IF ASSESS/CONTROL JUMPS YEAR LINE
        nextYearCtrlSched = np.full((nControlAreas), -1)
        ## GET PRESCRIPTIVE MONTHLY CONTROL SCHEDULE
        if keepYear:
            for count, (dummask, startYear, revisit, 
                controlIndicator, shp) in enumerate(rawdata.rodentControlList):
                lastControl = lastControlArray[count]
                nYearsSinceControl = year - lastControl

                ## CONSIDER PRESCRIPTIVE IF CONTROL NOT SCHED FROM t-1 FROM MAST OR RODENTS
                if mthlyCtrlSched[count] < 0:
                    controlCondition =  ((year == startYear) or 
                        (year > startYear and nYearsSinceControl >= revisit))
                    if controlCondition:
                        ## add control area to schedule
                        mthlyCtrlSched[count] = params.prescrptCtrlMth
                        print('Prescriptive', 'Area', count, 'year', year,
                            'month', params.prescrptCtrlMth, 
                            'Yr Since', nYearsSinceControl, 'revisit', revisit,
                            'condit', controlCondition)

        # Masting affects rodents the same year now
        #but year starts in Sept, spring cos that's when beech flowering starts
        mastT = np.random.rand() < params.mastPrEvent
        # mastT = False
        #for testing purposes: have all iterations mast in same year 
#        mastYrarr=np.array([1,7,10,11,15,19,23,26,30])
#        if (year_all in mastYrarr):
#            mastT=True
#        else:
#            mastT=False   
            
        if mastT:
            print('Masting Year', year_all)
            mastingMask = doMasting(rawdata, params, distArray, masthalfwinsize, 
                        beechMask)
            preyMastingMask = resampleRasterDown(mastingMask*1.0,prey_raster[0,0:,0:].shape, 
                                                 statMethod = RESAMPLE_AVERAGE, pixelRescale=1)
            #if masting year and doing control reactive to masting then assess prop mgmt areas masting
            if (params.reactivePropMgmtMasting > 0):
                for count, (controlMask, startYear, revisit, 
                    controlIndicator, shp) in enumerate(rawdata.rodentControlList):
                    if not controlIndicator:
                        continue
                    propControlMaskMasting[count] = (np.count_nonzero(mastingMask & 
                        controlMask & rawdata.rodentExtentMask) / np.count_nonzero(controlMask & rawdata.rodentExtentMask))
                    # propControlMaskMasting[count] = (np.count_nonzero(mastingMask & 
                    #     controlMask) / np.count_nonzero(controlMask))
#                    propControlMaskMasting[count,0] = (np.count_nonzero(mastingMask & controlMask)
#                                / np.count_nonzero(controlMask))
                    if propControlMaskMasting[count] >= params.reactivePropMgmtMasting:                   
                        mthlyCtrlSched[count] = params.mastCtrlMth
                        print('mast prp > thres', np.round(propControlMaskMasting[count], 2),
                              'year', year,  'Area', count)

      
        ##age the prey
        if year_all > 0:
#            prey_raster[4,:,:] = prey_raster[4,:,:] + prey_raster[3,:,:]
#            prey_raster[3,:,:] = prey_raster[2,:,:] 
            prey_raster[3,:,:] = prey_raster[3,:,:] + prey_raster[2,:,:] 
            prey_raster[2,:,:] = prey_raster[1,:,:] 
            prey_raster[1,:,:] = prey_raster[0,:,:]
            prey_raster[0,:,:] = 0
        
        for mth in range(12):
            
            tMth = (year_all * 12) + mth
            
            #update time since control
            timeSinceCtrl[rawdata.rodentExtentMask] += 1
            #print('year=', year_all, 'month=', mth, 'time since ctrl=', timeSinceCtrl[1,1])
            
            ## DO RODENT GROWTH, first have to do some bodging to get seasonal stuff 
            ## use different kmaps depending on if crash or mast year, needs to be this
            #order so if double mast, the mast overrrides using the crash k map
            # rodent_kMap = rawdata.paramRodentCCLU[rawdata.kClasses] * 1
            # rodent_kMth = rodent_kMap.copy()
            rodent_kMth = rawdata.seasAdj[rawdata.kClasses,mth]

#            print('yr', year_all, 'mth', mth, 'Unique rodent kmth', len(np.unique(rodent_kMth)), 
#                'raw seasAdj', rawdata.seasAdj[:,mth], 'shp', np.shape(rodent_kMth))

             # print('t since:', timeSinceCtrl[1,15:25]) 
            # print('kvals b4:', rodent_kMth[1,15:25]) 
            #bounce adjustment shifted to before mast and crash adjustments are made so that these override the rat bounce effect
            rodent_kMth = np.where(timeSinceCtrl < params.rodentBouncePeriod,
                                rodent_kMth * params.rodentBounceMult * np.exp(
                                params.rodentBounceDecay * timeSinceCtrl), rodent_kMth)

            # print('kvals after:', rodent_kMth[1,15:25])
            if mastT_1:
                rodent_kMth = np.where(oldMastingMask==True,
                    rawdata.crashSeasAdj[rawdata.kClasses,mth], rodent_kMth)
            # #bounce adjustment shifted to before mast  adjustments are made so that masting overrides the rat bounce effect
            # rodent_kMth = np.where(timeSinceCtrl<params.rodentBouncePeriod,params.rodentBounceMult*rodent_kMth ,rodent_kMth)
            if mastT:
                rodent_kMth = np.where(mastingMask==True,rawdata.mastSeasAdj[rawdata.kClasses,mth],rodent_kMth)
            # #bounce adjustment shifted to after mast and crash adjustments are made
            # rodent_kMth = np.where(timeSinceCtrl<params.rodentBouncePeriod,params.rodentBounceMult*rodent_kMth ,rodent_kMth)

            rodent_raster = doRodentGrowth(rawdata, params,
                rodent_raster, rodent_kMth, mth)

            #assess need for control (mast or high tracking rates) in params.reactiveAssessMth 
            #but don't actually apply control until a few months later
            #all in the same month for now until figure out how to vary
            #only do assessment and control after burnin (i.e. in a 'keep year')
            if (mth==params.reactiveAssessMth) and keepYear:
                ## If control reactive to TT rates do a TT assessment:
                if (params.threshold_TT < 1.0):
                    for count, (controlMask, startYear, revisit, controlIndicator, 
                        shp) in enumerate(rawdata.rodentControlList):
                        ## DON'T ASSESS MZ IF NOT CONTROL ZONE
                        if not controlIndicator:
                            continue
                        ## TEST IF ALREADY SCHED TO DO CONTROL WITH PRESCRIBE OR MAST
                        no_TT_Test = mthlyCtrlSched[count] >=0
                        # ## ASSESS IF NOTHING SCHED OR IF MONTH IS NOT TT CONTROL
                        # ## TT CONTROL MONTH SHOULD BE PRIORITISED OVER PRESCRIBED.
                        # no_TT_Test = mthlyCtrlSched[count] == rawdata.reactiveCtrlMth
                        ## NO NEED TO ASSESS RODENTS IF ALREADY DOING CONTROL 
                        if no_TT_Test:
                            continue
                        ## ASSESS TT RATE
                        nPixelsZone = rawdata.pixelsInControlAndBeechMask[shp]
                        nRodentsinControlAndBeechMask = np.sum(rodent_raster[rawdata.controlAndBeechMask[shp]])
                        rodentsInControlAndBeechMaskDensity = (nRodentsinControlAndBeechMask/nPixelsZone)
                        TT_rate = trackingTunnelRate(rodentsInControlAndBeechMaskDensity, 
                                params.g0_TT, params.sigma_TT, 
                                params.nights_TT, params.nTunnels, params.resolutions[0])
                        ## IF TT RATE EXCEEDS THRESHOLD                            
                        if TT_rate > params.threshold_TT:
#                            mthlyCtrlSched[count] = rawdata.reactiveCtrlMth
                            ## DOES CONTROL JUMP THE YEAR LINE????
                            ## SCHEDULE CONTROL IN THIS YEAR (T)
                            if not rawdata.jumpYearCtrl:                      
                                mthlyCtrlSched[count] = rawdata.reactiveCtrlMth
                            ## else SCHEDULED CONTROL JUMPS INTO THE NEXT YEAR (T+1)
                            else:
                                nextYearCtrlSched[count] = rawdata.reactiveCtrlMth
                            print('TT rate > threshold', np.round(TT_rate, 2),
                                'year', year, 'assess mth', mth, 'ctrl mth', 
                                rawdata.reactiveCtrlMth, 'area', count, 
                                'jumpYear', rawdata.jumpYearCtrl)

            # Rodent Control. Get the combined control masks for this year
            # It's a bit hard to work out the mask ahead of time with the
            # revisits and reactive control, so we do it here.
            control_mth = False
            if (mth in mthlyCtrlSched):
                schedControlMask = None
                if keepYear:
                    areaList = np.where(mthlyCtrlSched == mth)[0]
#                    areaList = [i for i, mthlyCtrlSched in enumerate(mthlyCtrlSched) if mthlyCtrlSched == mth ]
                    for i in range(len(areaList)):
                        mArea=areaList[i]
                        cmask = rawdata.rodentControlList[mArea][0]
 
                        if schedControlMask is None:
                            schedControlMask = cmask.copy()
                        else:
                            schedControlMask |= cmask
                        ## ADD CONTROL ZONES TO CONTROLTHISYEAR FOR RESULTS.
                        if loopIter == 0:
                            if controlThisYear is None:
                                controlThisYear = cmask.copy()
                            else:
                                controlThisYear |= cmask
    
                        lastControlArray[mArea] = year
                        results.controlCount += 1
                control_mth = schedControlMask is not None and keepYear # and past burn in
        #        print('year', year, 'control_mth', control_mth)
                # if reactiveControlMask is not None:
                #     # we are defintely doing control due to reactive control
                #     control_mth = True
    
            if control_mth:
                
                # if schedControlMask is not None:
                #     # add in any reactive control
                #     if reactiveControlMask is not None:
                #         schedControlmMask = schedControlMask | reactiveControlMask
                # else:
                #     # the control is just due to reactive control
                #     schedControlMask = reactiveControlMask
    
                (rodent_raster, nToxicRodents) = doControl(rawdata, params, 
                    schedControlMask, rodent_raster)
                timeSinceCtrl[schedControlMask] = 0
                #print('time since ctrl =', timeSinceCtrl[1,1:40])

    
            else:
                # no control
                schedControlMask = np.zeros_like(beechMask)
                nToxicRodents = np.zeros_like(beechMask, dtype = np.uint16)
 
    
             ## stoat population grows
            (stoat_raster, rodent_raster_stoat) = calcStoatPopulation(stoat_raster, 
                    rodent_raster, nToxicRodents, rawdata.stoatExtentMask, params, 
                    mastT_1, schedControlMask, control_mth, keepYear, 
                    stoatIslandMask, haPixelRescale, stoatIslandKArray, mth)
    
    #        print('stoat', np.sum(stoat_raster[rawdata.stoatExtentMask]),
    #            'rodents', np.sum(rodent_raster[rawdata.rodentExtentMask]),
    #            'extmask', type(rawdata.rodentExtentMask), 
    #            rawdata.rodentExtentMask[100:120,100:105])
    
    
            ### Do Dispersal of three species
            # Rodent immigrat/emigrat occurs after stoats respond to rodents
            if params.rodentSeasDisp[mth]:
                rodent_raster = doRodentDispersal(rawdata, params, rodent_raster, 
                            rodent_kMth, rodentEmigrationWindowSizePxls)
    
    
            ## Do stoat dispersal follows rodent dispersal
            rodent_raster_stoat = resampleRasterDown(rodent_raster, 
                rawdata.stoatExtentMask.shape, statMethod = RESAMPLE_SUM, 
                pixelRescale = haPixelRescale)
            ## ASSIGN NOTIONAL RODENT DENSITY TO ISLANDS
            rodent_raster_stoat[stoatIslandMask] = stoatIslandKArray
    
    
    
    #        ## CONVERT TO RODENTS PER HECTARE
    #        rodent_raster_stoat = rodent_raster_stoat / nHectInRodent
    
    
    
            ## STOAT DISPERSAL
            if params.stoatSeasDisp[mth]:
                stoat_raster = doStoatDispersal(rawdata, params, stoat_raster,
                    rodent_raster_stoat, stoatEmigrationWindowSizePxls)
    
    
    
            ## TODO: WHAT RODENT COUNTS DO WE WANT? PER HA OR KM?
    
            # ## STORE COPY OF RODENT_RASTER FOR T-lag CALC FOR STOAT DYNAMICS
            # if mastT:
            #     rodentRasterStoatT_lag = rodent_raster_stoat.copy()
    
    
    
    
    
    
    
            # prey population growth after stoat and rodent dispersal
            ## IF PREY RESOL != STOAT RESOL, HAVE TO REPLACE "rodent_raster_stoat" in
            ## THE FOLLOWING FX AND GET 'rodent_raster_prey' WITHIN THE FX. THIS IS
            ## FOR THE COMPETITION EFFECT BETWEEN RODENTS AND PREY.
            prey_raster = doPreyGrowth(prey_raster, stoat_raster, params, 
                    rawdata.preyExtentMask, rodent_raster_stoat, nHectInRodent,
                    preyMastingMask, preyRecDecay_1D, preySurvDecay_1D, mth, 
                    rawdata.pLeadDeath3D, keepYear)
    
    
    
    
            ## Prey dispersal
            if params.preySeasDisp[mth]:            
                prey_raster = doPreyDispersal(rawdata, prey_raster, params, 
                            rawdata.preyExtentMask, preyEmigrationWindowSizePxls)
    
    
            #populate storage arrays with just densities each mth
            populateResultDensity(rodent_raster, rawdata.rodentExtentMask, 
                            rawdata.rodentControlList, rawdata.rodentAreaDictByMgmt, 
                            results.rodentDensity_2D_mth, stoat_raster, 
                            rawdata.stoatExtentMask, rawdata.stoatSpatialDictByMgmt, 
                            rawdata.stoatAreaDictByMgmt,results.stoatDensity_2D_mth, 
                            prey_raster, rawdata.preyExtentMask,rawdata.preySpatialDictByMgmt, 
                            rawdata.preyAreaDictByMgmt, results.preyDensity_2D_mth, tMth)

            ## Populate storage arrays with updated densities and maps - full output 
            #only once per year in Jan (mth 4)
            if (mth==4):
                populateResultArrays(loopIter, mastingMask, rodent_raster, 
                    rawdata.rodentExtentMask, rawdata.rodentControlList, 
                    rawdata.rodentAreaDictByMgmt, results.rodentDensity_2D,
                    stoat_raster, rawdata.stoatExtentMask, rawdata.stoatSpatialDictByMgmt, 
                    rawdata.stoatAreaDictByMgmt, results.stoatDensity_2D, prey_raster, 
                    rawdata.preyExtentMask, rawdata.preySpatialDictByMgmt, rawdata.preyAreaDictByMgmt, 
                    results.preyDensity_2D, results.popAllYears_3D, year, year_all, keepYear)


#                populateResultArrays(loopIter, mastingMask, schedControlMask, rodent_raster, 
#                    rawdata.rodentExtentMask, rawdata.rodentControlList, 
#                    rawdata.rodentAreaDictByMgmt, results.rodentDensity_2D,
#                    stoat_raster, rawdata.stoatExtentMask, rawdata.stoatSpatialDictByMgmt, 
#                    rawdata.stoatAreaDictByMgmt, results.stoatDensity_2D, prey_raster, 
#                    rawdata.preyExtentMask, rawdata.preySpatialDictByMgmt, rawdata.preyAreaDictByMgmt, 
#                    results.preyDensity_2D, results.popAllYears_3D, year, year_all, keepYear)

            if (mth == 11) and (loopIter == 0) and keepYear:
                populateResultControl(year, results.popAllYears_3D, controlThisYear)

        #print('time since ctrl =', timeSinceCtrl[1,1:40])
        ### IF BEYOND BURN IN PERIOD, STORE THE RESULTS
        if keepYear:
            ## update the year
            year += 1

        # update the masting status of T_1 for next year
        mastT_1 = mastT
        oldMastingMask = mastingMask.copy()
###        if mastT_1:
###            propControlMaskMasting[:,1] = propControlMaskMasting[:,0]

    return results
        

def calcStoatPopulation(stoat_raster, rodent_raster, nToxicRodents, stoatMask, params, 
        mastT_1, schedControlMask, control_mth, keepYear, stoatIslandMask, 
        haPixelRescale, stoatIslandKArray, mth):
    """
    ## Do the stoat processes at stoat resolution: control, growth, dispersal
    """
    # get rodent population and toxic rodents at stoat resolution
    # it is rodents per ha -> sum and divid by 100 
    # rodent_raster_stoat = resampleRasterDown(rodent_raster, 
    #             stoatMask.shape, statMethod = RESAMPLE_SUM, 
    #             pixelRescale = haPixelRescale)
    #lets try that per stoat home range so don't just get zero rodents all the time
    rodent_raster_stoat = resampleRasterDown(rodent_raster, 
                stoatMask.shape, statMethod = RESAMPLE_SUM, 
                pixelRescale = 1)
####    ## ASSIGN NOTIONAL RODENT DENSITY TO ISLANDS - PRE-CALCULATED
####    rodent_raster_stoat[stoatIslandMask] = stoatIslandKArray

#    print('rodent', np.sum(rodent_raster), 'rodentStoat', np.sum(rodent_raster_stoat))


    ## IF CONTROL IS DONE
    if control_mth:
        # Eqn 25        # resample toxic rodents at stoat resolution
        toxic_raster_stoat = resampleRasterDown(nToxicRodents, 
                stoatMask.shape, statMethod = RESAMPLE_SUM, pixelRescale = 1)
        
        ## stoats eat toxic rodents
        stoat_raster = eatToxicRodents(stoat_raster, toxic_raster_stoat, params)

    # ### Stoats decline more slowly than rodents following a mast,
    # ### Calc ave of current R_T and R_T_1
    # if mastT_1:
    #     ## GET AVERAGE FOR T AND T-1
    #     rodent_t = (rodent_raster_stoat + rodentRasterStoatT_lag) / 2.0
    #     ## USE CURRENT YEAR RODENT POP
    # else:
    #     rodent_t = rodent_raster_stoat.copy()
    rodent_t = rodent_raster_stoat.copy()    
    ###################################################### IS THIS DOUBLING FROM ABOVE?
    ## SET RODENT DENSITY LOW ON ISLANDS AND TRAPPED AREAS
    rodent_t[stoatIslandMask] = stoatIslandKArray       # params.islandK
    ######################################################
    ## MAKE 1-D ARRAYS FOR POPULATION UPDATE
    rodent_t = rodent_t[stoatMask]     ## RODENTS PER HA
    stoat_t = stoat_raster[stoatMask]


#    print('stoat grow 0', np.mean(stoat_t), np.min(stoat_t), np.max(stoat_t),
#        'sRas', stoat_raster[75:80, 75:80], rodent_raster_stoat[75:80, 75:80])


    ## WHERE HAVE ZERO RODENTS IN A KM CELL, SET TO 0.5 RODENT SO DON'T DIVIDE BY 0
    rodent_t = np.where(rodent_t < 0.25, 0.025, rodent_t)
    

    ## Stoat popn growth
    seasRec = params.stoatSeasRec[mth]
    pSurv = (params.stoatSurv * np.exp(-((stoat_t/(params.stoatSurvDDcoef*rodent_t))
                                         **params.stoatTheta)))
    recRate = ((np.exp(seasRec * params.stoatIRR)-1) * 
        np.exp(-((stoat_t/(params.stoatRecDDcoef*rodent_t))**params.stoatTheta)))
    stoat_t = stoat_t * pSurv * (1 + recRate) 


#    stoat_t = np.random.poisson(stoat_t)
    ## EQN 28: ADD STOCHASTICITY GAUSSIAN PROCESS
    stoat_t = (np.exp(rng.normal(np.log(stoat_t + 1.0), 
                params.stoatPopSD)) - 1.0)

#    print('stoat_t min, max', np.min(stoat_t),
#        np.max(stoat_t))


    ## FIX UP INAPPROPRIATE VALUES WHEN USING NORMAL DISTRIBUTION
    stoat_raster[stoatMask] = np.round(stoat_t, 0).astype(int)
    stoat_raster = np.where(stoat_raster < 0, 0, stoat_raster)
    stoat_raster[~stoatMask] = 0


#    print('stoat grow 11111', np.mean(stoat_t), np.min(stoat_t), np.max(stoat_t),
#        'sRas', stoat_raster[75:80, 75:80], rodent_raster_stoat[75:80, 75:80])

#    stoat_raster[stoatMask] = stoat_t
    ## RETURN STOAT RASTER
    return(stoat_raster, rodent_raster_stoat)        # Return stoat raster



def eatToxicRodents(stoat_raster, toxic_raster_stoat, params):
    """
    ## update stoat_raster by consuming toxic rodents
    """
    # Eqn 27        # probability of encounter 
    pEnc = 1.0 - np.exp(-params.pEncToxic * toxic_raster_stoat)
    # probability individ. stoat eating a toxic rodent
    pEat = params.pEatEncToxic * pEnc
    # Eqn 26     # update stoat_raster
    stoat_raster = rng.binomial(stoat_raster, (1.0 - pEat))
    return(stoat_raster)                   


def doRodentDispersal(rawdata, params, rodent_raster, rodent_kMth,
                rodentEmigrationWindowSizePxls):

    # First Emigrants out of each cell
    # Eqn 18
    pEm = params.gammaProbEmigrate[0]
    probRodentEmigrate = np.where(rawdata.rodentExtentMask, 
        (1.0 - np.exp(-rodent_raster * pEm)) *  
        np.exp(-rodent_kMth * params.tauImmigrate[0]), 0.0)

#    probRodentEmigrate[rawdata.rodentExtentMask] = ((1.0 - 
#        np.exp(-rodent_raster[rawdata.rodentExtentMask] * pEm)) * 
#        np.exp(-rodent_kMth[rawdata.rodentExtentMask] * params.tauImmigrate[0]))

    # Eqn 17
    rodent_emigrant_raster = rng.binomial(rodent_raster, probRodentEmigrate)



    # updates rodent_change_raster
    # Eqn 19-21
    rodent_immigrant_raster = calcImmigration(rodent_emigrant_raster, 
        rodentEmigrationWindowSizePxls, params.deltaImmigrate[0], 
        rawdata.rodentExtentMask, params.resolutions[0], 
        rodent_raster, pEm, rawdata.rodentPercentArea, 
        params.tauImmigrate[0], rodent_kMth)

    ### Eqn 16:  UPDATE TO rodent_raster by immigrants/emigrants
    rodent_raster += rodent_immigrant_raster
    rodent_raster -= rodent_emigrant_raster

    return(rodent_raster)


def doStoatDispersal(rawdata, params, stoat_raster, rodent_raster_stoat,
                stoatEmigrationWindowSizePxls):
    """
    ## calc emigration and immigration by stoats
    """
    # First Emigrants out of each cell
    pEm = params.gammaProbEmigrate[1]
    rodentRasterStoat = np.where(rodent_raster_stoat < 0.25, 0.25, rodent_raster_stoat) 
    # Eqn 34
    probstoatEmigrate = np.where(rawdata.stoatExtentMask,
        (1.0 - np.exp(-stoat_raster*pEm)) * 
        np.exp(-rodentRasterStoat*params.tauImmigrate[1]), 0.0)

    # Eqn 33
    stoat_emigrant_raster = rng.binomial(stoat_raster, probstoatEmigrate)

    # updates stoat_change_raster
    # Eqn 34-36
    stoat_immigrant_raster = calcImmigration(stoat_emigrant_raster, 
        stoatEmigrationWindowSizePxls, 
        params.deltaImmigrate[1], rawdata.stoatExtentMask, 
        params.resolutions[1], stoat_raster, pEm, rawdata.stoatPercentArea,
        params.tauImmigrate[1], rodentRasterStoat)

    ### EQN 32 Update stoat_raster by immigrants/emigrants
    stoat_raster += stoat_immigrant_raster
    stoat_raster -= stoat_emigrant_raster

    return(stoat_raster)



@njit
def calcImmigration(emigrant_raster, winsize, 
            deltaImmigrate, mask, resolution, raster, gammaPara,
            areaCorrection, tauPara, kMap):
    """
    Note: same for rodents and stoats, but for prey correct for area in cell
    """
    halfwinsize = int(winsize / 2)

    (ysize, xsize) = emigrant_raster.shape
    probEmigrate = np.empty((winsize, winsize))
    immigrant_raster = np.zeros(emigrant_raster.shape, dtype=np.uint32)

    for x in range(xsize):
        for y in range(ysize):
            if not mask[y, x]:
                continue

            if emigrant_raster[y, x] == 0:
                # no rodents to move out of here
                continue

            # offset into dist array - deal with top and left edges
            xdistoff = 0
            ydistoff = 0
            # top left x
            tlx = x - halfwinsize
            if tlx < 0:
                xdistoff = -tlx
                tlx = 0
            tly = y - halfwinsize
            if tly < 0:
                ydistoff = -tly
                tly = 0
            brx = x + halfwinsize
            if brx > xsize - 1:
                brx = xsize - 1
            bry = y + halfwinsize
            if bry > ysize - 1:
                bry = ysize - 1
            # calculate the relative probability of emigrating to each cell
            sumRelProbEm = 0.0
            for cx in range(tlx, brx):
                for cy in range(tly, bry):
                    if not mask[cy, cx]:
                        continue
                    distx = (x - cx) * resolution
                    disty = (y - cy) * resolution
                    dist = np.sqrt(distx*distx + disty*disty)


#####                    if dist < 50.0:
#####                        relProbEm = 0.0     # force emigrants to leave.
#####                    elif dist >= 50.0:
                    if tauPara == 0.0:
                        ## EQN 45: PREY DISPERSAL
                        relProbEm = (np.exp(-(raster[cy,cx] * gammaPara) /
                            areaCorrection[cy, cx]) *
                            np.exp(-dist / 1000.0 * deltaImmigrate))
                    else:
                        ## EQN 20 AND 36: RODENT AND STOAT DISPERSAL
                        relProbEm = (np.exp(-raster[cy,cx] * gammaPara / 
                            areaCorrection[cy, cx]) * 
                            (1.0 - np.exp(-kMap[cy, cx] * tauPara)) * 
                            np.exp(-dist / 1000.0 * deltaImmigrate))

                    sumRelProbEm += relProbEm
                    probEmigrate[ydistoff + cy - tly, xdistoff + cx - tlx] = relProbEm

            sumProbsUsed = 0.0
            emRemain = emigrant_raster[y, x] 
            for cx in range(tlx, brx):
                for cy in range(tly, bry):
                    if not mask[cy, cx]:
                        continue
                    if emRemain <= 0:
                        break
                    subx = xdistoff + cx - tlx
                    suby = ydistoff + cy - tly
                    # change probEmigrate to be the probability (sum to 1)
                    # Eqn 21, 37, 46

                    if sumRelProbEm <= 0.0:
                        print('bad sum relprobs =', sumRelProbEm)

                    absprob = probEmigrate[suby, subx] / sumRelProbEm

                    prob_fg = absprob / (1.0 - sumProbsUsed)
                    
                    if prob_fg > 1.0:
                        prob_fg = 1.0
                    elif prob_fg < 0.0:
                        prob_fg = 0.0

                    # Eqn 19, 35, 44    Multinomial draw 
                    X_fg = np.random.binomial(emRemain, prob_fg)

                    sumProbsUsed += absprob

                    immigrant_raster[cy, cx] += X_fg
                    emRemain -= X_fg
                if emRemain <= 0:
                    break

    return immigrant_raster


def doControl(rawdata, params, schedControlMask, rodent_raster):
    """
    ## calc number of rodents eating toxin and dieing
    """
    # Eqn 14
    nToxicRodents = np.where(schedControlMask, 
                rng.binomial(rodent_raster, params.rodentProbEatBait), 0)
    ## update rodent pop. with mortality
    # Eqn 15
    rodent_raster = rodent_raster - nToxicRodents
    return rodent_raster, nToxicRodents


def doRodentGrowth(rawdata, params, rodent_raster, rodent_kMth, mth):
    """
    ## calc rodent growth
    """
    ## mask by rawdata.rodentExtentMask
    mu = np.zeros_like(rodent_raster, dtype=float)
    #get seas adj to rec
    seasRec = params.rodentSeasRec[mth]
    ## RODENTS AND SEEDS FOR MTH MASKED BY EXTENT & seas adj
    rodent_t = rodent_raster[rawdata.rodentExtentMask]
    seed_t = rodent_kMth[rawdata.rodentExtentMask]

    ## SURVIVAL PROBABLITY
    pSurv = params.rodentSurv * np.exp(-((rodent_t/(seed_t*params.rodentSurvDDcoef))**params.rodentTheta))
    ## RECRUITMENT RATE
    recRate = ((np.exp(seasRec * params.rodentIRR)-1) * 
               np.exp(-((rodent_t/(seed_t*params.rodentRecDDcoef))**params.rodentTheta)))
    ## Eqn. 13 POPULATE MU RASTER
    mu[rawdata.rodentExtentMask] = rodent_t * pSurv * (1 + recRate) 

    # Eqn 12: Add stochasticity
    rodent_raster = rng.poisson(mu, rodent_raster.shape)
    # rodent_raster = np.round(mu).astype(int)
    ## RETURN RASTER
    return rodent_raster





def doPreyGrowth(prey_raster, stoat_raster, params, mask, 
        rodent_raster_prey, nHectInRodent, mastMsk, preyRecDecay_1D, 
        preySurvDecay_1D, mth, pLeadDeath3D, keepYear):
    """
    ## calc prey population growth by pixel
    """
###    s = np.exp(-params.preyPsi * stoat_raster)
    # get rodent population and toxic rodents at prey resolution
    # it is rodents per 4 ha because we take the average

    ###################################################################
    ###################################################################
    ## DO THE FOLLOWING ONLY IF STOAT AND PREY RESOL ARE NOT EQUAL
###    rodent_raster_prey = resampleRasterDown(rodent_raster, 
###                mask.shape, statMethod = RESAMPLE_AVERAGE)
###    rodentsPerHect = rodent_raster_prey / nHectInRodent

    ###################################################################
    ###################################################################
    ## REMOVE COMPETITION EFFECT FOR NOW
    # competition effect of rodents
###    c = np.exp(-params.competEffect * rodentsPerHect)
    # competition and predation survival 
###    CS = s * c
    ###################################################################
    ###################################################################

    ###################################################################
    ###################################################################
    ## DO THE FOLLOWING ONLY IF STOAT AND PREY RESOL ARE NOT EQUAL
#    stoat_raster_prey = resampleRasterDown(stoat_raster, mask.shape, 
#                       statMethod = RESAMPLE_SUM, pixelRescale = 1)

    ## MAKE 1-D ARRAYS FOR POPULATION UPDATE

## LEAD POISONING EFFECT ON SURVIVAL
    if pLeadDeath3D is not None:
        pLead_t = pLeadDeath3D[:, mask]
    
    prey_t = prey_raster[:,mask]
    stoat_t = stoat_raster[mask]
    rodent_t = rodent_raster_prey[mask]
    rodentSwitchMult_t = np.where(rodent_t<=params.rodentThresh,params.stoatMult,1)
    mastInd = np.where(mastMsk[mask]>0,1,0)
    
    for x in range(4):
        pSurv = (params.preySurv[x,mastInd] * np.exp(-((prey_t[4,:]/(preySurvDecay_1D))**params.preyTheta)) * 
                   np.exp(-params.preyEtaStoat[x] * stoat_t * rodentSwitchMult_t) * 
                   np.exp(-params.preyEtaRodent[x] * rodent_t)) 
        if pLeadDeath3D is not None:
            if params.removeLeadAction:
                if not keepYear:
                    if x<3:
                        pSurv = pSurv * (1.0 - pLead_t[0])
                    else:
                        pSurv = pSurv * (1.0 - pLead_t[1])
            else:
                if x<3:
                    pSurv = pSurv * (1.0 - pLead_t[0])
                else:
                    pSurv = pSurv * (1.0 - pLead_t[1])


        prey_t[x,:]=rng.binomial(prey_t[x,:],pSurv)
          
    ###  nb: total popn (prey_t[4,:]) is not updated b4 calculating dd-recruitment -
    #### this is cos clutch size/propn breeding is determined months earlier c.f. here
    #### which is when fledging/recruitment occurs
    for x in range(1,4):
          recRate =  (params.preySeasRec[mth,mastInd]*params.preyPropBreedpa[mastInd]*params.preyFec[x,mastInd] * 
                      np.exp(-((prey_t[4,:]/(preyRecDecay_1D))**params.preyTheta)) * 
                      np.exp(-params.preyPsiStoat * stoat_t) *
                      np.exp(-params.preyPsiRodent * rodent_t))         
          prey_t[0,:]=prey_t[0,:]+rng.poisson(prey_t[x,:]*recRate)  #new recruits added to zero age class            
    
    prey_t[4,:]=np.sum(prey_t[:4,:],axis=0)  #update total pop size
    
    # ## FIX UP INAPPROPRIATE VALUES
    # prey_raster[mask] = np.round(prey_t, 0).astype(int)
    # prey_raster = np.where(prey_raster < 0, 0, prey_raster)
    # prey_raster[~mask] = 0
    
    prey_raster[:,mask]=prey_t[:,:]
    prey_raster[:,~mask]=0
        
        # ## RETURN PREY RASTER
    
#     prey0_t = prey_raster[0,mask]
#     prey1_t = prey_raster[1,mask]
#     prey2_t = prey_raster[2,mask]
#     prey3_t = prey_raster[3,mask]
# #    prey4_t = prey_raster[4,mask]
#     preyN_t = prey_raster[4,mask]
# #    preyN_t = prey_raster[5,mask]
#     stoat_t = stoat_raster[mask]
# #   stoat_t = stoat_raster_prey[mask]
#     rodent_t = rodent_raster_prey[mask]
#     rodentSwitchMult_t = np.where(rodent_t<=params.rodentThresh,params.stoatMult,1)
#     #not sure why so many 3s in print output yet summary file "monthlyDensities.csv" doesn't show this???
#     #oh is prob cos monthlyDensities are not per ha - yes they are - maybe something to chk... 
#     #rodent raster is all integers?! understandable when at rodent resolution but thought 
#     #would be diff for stoat or pre resol rodent_raster_stoat cos averages???
#     # if (mth==1):
#     #     print(np.array2string(rodent_t[100:115], precision=2))
#     #     print(np.array2string(rodentSwitchMult_t[100:115], precision=2))
#     #     print('num <=0.5', np.sum(rodent_t <=0.5), 'num >0.5', np.sum(rodent_t >0.5))
#     #     print('num <=0.5', np.sum(rodentSwitchMult_t==3), 'num >0.5', np.sum(rodentSwitchMult_t==1))
#     mastEff_t = mastMsk[mask]
#     mastEff_t = np.where(mastEff_t>0,mastEff_t*params.preyMastMultFec,1)
    
    
#     ## PREY POPN DYNAMICS        
#     # pSurv = params.preySurv[0] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#     # prey0_t = rng.binomial(prey0_t, pSurv)
#     # pSurv = params.preySurv[1] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#     # prey1_t = rng.binomial(prey1_t, pSurv)
#     # pSurv = params.preySurv[2] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#     # prey2_t = rng.binomial(prey2_t, pSurv)
#     # pSurv = params.preySurv[3] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#     # prey3_t = rng.binomial(prey3_t, pSurv)
#     # pSurv = params.preySurv[4] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#     # prey4_t = rng.binomial(prey4_t, pSurv)
#     # seasRec = params.preySeasRec[mth]
#     # recRate = (seasRec * params.preyProd *np.exp(-((preyN_t/(preyRecDecay_1D))**params.preyTheta))
#     #            * np.exp(-params.preyPsi * stoat_t))
#     # prey0_t = prey0_t + rng.poisson(prey4_t * recRate) #fledglings get added to zero age class 
#     # preyN_t = prey0_t + prey1_t + prey2_t + prey3_t + prey4_t  #sum to get total popn
 

#     ## YEAR 0-1
#     pSurv = (params.preySurv[0] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#          * np.exp(-params.preyEtaStoatJuv * stoat_t * rodentSwitchMult_t)* np.exp(-params.preyEtaRodentJuv * rodent_t))
#     ## LEAD POISONING EFFECT ON SURVIVAL
#     if pLeadDeath3D is not None:
#         pSurv = pSurv * (1.0 - pLead_t[0])
#     prey0_t = rng.binomial(prey0_t, pSurv)

#     ## YEAR 1-2
#     pSurv = (params.preySurv[1] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#          * np.exp(-params.preyEtaStoatAd * stoat_t * rodentSwitchMult_t)* np.exp(-params.preyEtaRodentAd * rodent_t))    
#     ## LEAD POISONING EFFECT ON SURVIVAL
#     if pLeadDeath3D is not None:
#         pSurv = pSurv * (1.0 - pLead_t[0])
#     prey1_t = rng.binomial(prey1_t, pSurv)

#     ## YEAR 2-3
#     pSurv = (params.preySurv[2] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#          * np.exp(-params.preyEtaStoatAd * stoat_t * rodentSwitchMult_t)* np.exp(-params.preyEtaRodentAd * rodent_t))    
#     ## LEAD POISONING EFFECT ON SURVIVAL
#     if pLeadDeath3D is not None:
#         pSurv = pSurv * (1.0 - pLead_t[0])
#     prey2_t = rng.binomial(prey2_t, pSurv)

#     ## YEAR > 3
#     pSurv = (params.preySurv[3] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
#          * np.exp(-params.preyEtaStoatAd * stoat_t * rodentSwitchMult_t)* np.exp(-params.preyEtaRodentAd * rodent_t))    
#     ## LEAD POISONING EFFECT ON SURVIVAL
#     if pLeadDeath3D is not None:
#         pSurv = pSurv * (1.0 - pLead_t[1])
#     prey3_t = rng.binomial(prey3_t, pSurv)

# #    pSurv = (params.preySurv[4] * np.exp(-((preyN_t/(preySurvDecay_1D))**params.preyTheta))
# #         * np.exp(-params.preyEtaStoatAd * stoat_t * rodentSwitchMult_t)* np.exp(-params.preyEtaRodentAd * rodent_t))    
# #    prey4_t = rng.binomial(prey4_t, pSurv)
#     seasRec = params.preySeasRec[mth]
#     recRate = ((np.exp(seasRec * params.preyIRR * mastEff_t)-1) *
#                np.exp(-((preyN_t/(preyRecDecay_1D))**params.preyTheta)) * 
#                np.exp(-params.preyPsiStoat * stoat_t) *
#                np.exp(-params.preyPsiRodent * rodent_t))

#     prey0_t = prey0_t + rng.poisson(prey3_t * recRate) #fledglings get added to zero age class 
#     preyN_t = prey0_t + prey1_t + prey2_t + prey3_t   #sum to get total popn

#     #shouldn't need to do this if drawing from binomial and poisson??
#     ## ADD STOCHASTICITY GAUSSIAN PROCESS
#     # # Eqn. 38
#     # prey_t = np.exp(rng.normal(np.log(prey_t + 1.0), 
#     #             params.preyPopSD)) - 1.0
#     # ## FIX UP INAPPROPRIATE VALUES
#     # prey_raster[mask] = np.round(prey_t, 0).astype(int)
#     # prey_raster = np.where(prey_raster < 0, 0, prey_raster)
#     # prey_raster[~mask] = 0
    
#     prey_raster[0,mask]=prey0_t
#     prey_raster[1,mask]=prey1_t
#     prey_raster[2,mask]=prey2_t
#     prey_raster[3,mask]=prey3_t
#     prey_raster[4,mask]=preyN_t
#     prey_raster[:,~mask]=0
    
    # ## RETURN PREY RASTER
    return(prey_raster)


def doPreyDispersal(rawdata, prey_raster, params, mask, preyEmigrationWindowSizePxls):
    # First Emigrants out of each cell
    # Eqn 43
    ## EMIGRATION
    preyDispRaster = prey_raster[0,:,:]
    probPreyEmigrate = np.zeros_like(preyDispRaster)
    probPreyEmigrate[mask] = (1.0 - np.exp(-params.gammaProbEmigrate[2] * 
                preyDispRaster[mask] / rawdata.preyCorrectionK[mask]))
    # Eqn 42
    prey_emigrant_raster = rng.binomial(preyDispRaster, probPreyEmigrate)


    # Eqn 44-46
    ## IMMIGRATION
    prey_immigrant_raster = calcImmigration(prey_emigrant_raster, 
        preyEmigrationWindowSizePxls, params.deltaImmigrate[2], mask, 
        params.resolutions[2], preyDispRaster, params.gammaProbEmigrate[2], 
        rawdata.preyCorrectionK, params.tauImmigrate[2], rawdata.preyKDummy)


    ## POPULATE RASTERS
    # Eqn. 47
    preyDispRaster += prey_immigrant_raster
    preyDispRaster -= prey_emigrant_raster
    prey_raster[0,:,:]=preyDispRaster
    return prey_raster


def doMasting(rawdata, params, distArray, halfwinsize, beechMask):
    # risk of masting for a given cell
    a, b = params.mastCellParams
    mastingRisk = rng.gamma(a, b, (rawdata.rodentNrows, rawdata.rodentNcols))
    if np.isnan(mastingRisk).any():
        raise ValueError("NaNs in mastingRisk")

    # aggregate values weighted by distance
    absolMastingRisk = np.zeros((rawdata.rodentNrows, rawdata.rodentNcols))

    maxMastingRisk = np.zeros((rawdata.rodentNrows, rawdata.rodentNcols))
    spatialRandMast(mastingRisk, maxMastingRisk, absolMastingRisk, beechMask, 
        halfwinsize, params.mastSpatialSD, distArray)

    if np.isnan(absolMastingRisk).any():
        raise ValueError("NaNs in absolMastingRisk")

    a, b = params.mastProportionParams
    propMastingCells = rng.beta(a, b)

    # invert proportion so we can use greater than
    # Eqn 2
    quantLevels = mquantiles(absolMastingRisk[beechMask], 
            prob=[1.0 - propMastingCells])

    mastingMask = (absolMastingRisk >= quantLevels)
    mastingMask[~beechMask] = False

    # for debugging
#    driver = gdal.GetDriverByName('GTiff')
#    ds = driver.Create('mast_%d.tif' % year, rawdata.rodentNcols, rawdata.rodentNrows, 
#            1, gdal.GDT_Byte)
#    ds.SetProjection(NZTM_WKT)
#    ds.SetGeoTransform(rawdata.rodentGeoTrans)
#    band = ds.GetRasterBand(1)
#    band.WriteArray(np.where(mastingMask, 1, 0))
#    del ds

    return mastingMask

@njit
def makeDistArray(winsize, rho, halfwinsize, resolution):
    """
    Pre calculate this so we don't have to do it each time
    """
    distWt = np.empty((winsize, winsize))
    
    for x in range(winsize):
        for y in range(winsize):
            distx = (x - halfwinsize) * resolution
            disty = (y - halfwinsize) * resolution
            dist = np.sqrt(distx*distx + disty*disty)
            # Eqn 4
            distWt[y, x] = np.exp(-(dist**2.0) / 2.0 / rho**2)

    return distWt


@njit
def spatialRandMast(randRisk, maxRisk, absolRisk, mask, halfwinsize, 
        mastSpatialSD, distWt):
    """
    return a spatial effect raster se1 from the input raster and mask
    """
    (ysize, xsize) = absolRisk.shape
    for x in range(xsize):
        for y in range(ysize):
            if not mask[y, x]:
                absolRisk[y, x] = 0
                continue
            # offset into dist array - deal with top and left edges
            xdistoff = 0
            ydistoff = 0
            # top left x
            tlx = x - halfwinsize
            if tlx < 0:
                xdistoff = -tlx
                tlx = 0
            tly = y - halfwinsize
            if tly < 0:
                ydistoff = -tly
                tly = 0
            brx = x + halfwinsize
            if brx > xsize - 1:
                brx = xsize - 1
            bry = y + halfwinsize
            if bry > ysize - 1:
                bry = ysize - 1
            maxRisk_xy = 0.0
            for cx in range(tlx, brx):
                for cy in range(tly, bry):
                    if not mask[cy, cx]:
                        continue
                    risk_cxcy = randRisk[cy, cx]    # [ydistoff + cy - tly, xdistoff + cx - tlx]
                    if risk_cxcy > maxRisk_xy:
                        maxRisk_xy = risk_cxcy
            ## Eqn 3
            ## add random variation
            maxRisk[y, x] = np.exp(np.random.normal(np.log(maxRisk_xy), mastSpatialSD))
    ## loop thru again to get mean to smooth surface
    for x in range(xsize):
        for y in range(ysize):
            if not mask[y, x]:
                absolRisk[y, x] = 0
                continue
            # offset into dist array - deal with top and left edges
            xdistoff = 0
            ydistoff = 0
            # top left x
            tlx = x - halfwinsize
            if tlx < 0:
                xdistoff = -tlx
                tlx = 0
            tly = y - halfwinsize
            if tly < 0:
                ydistoff = -tly
                tly = 0
            brx = x + halfwinsize
            if brx > xsize - 1:
                brx = xsize - 1
            bry = y + halfwinsize
            if bry > ysize - 1:
                bry = ysize - 1
            ### distance weighted average
            sumWtRisk = 0.0
            sumWt = 0.0
            for cx in range(tlx, brx):
                for cy in range(tly, bry):
                    if not mask[cy, cx]:
                        continue
                    wt = distWt[ydistoff + cy - tly, xdistoff + cx - tlx]
                    sumWtRisk += maxRisk[cy, cx] * wt
                    sumWt += wt
            if sumWt != 0 and sumWtRisk != 0:
                absolRisk[y, x] = sumWtRisk / sumWt
               #if np.isnan(absolRisk[y, x]):
                #    print(x, y, sumWtRisk, sumWt, randRisk[cy, cx])
            else:
                absolRisk[y, x] = 0
    return absolRisk


@njit
def resampleRasterDown(inarray, outSize, statMethod, pixelRescale):
    """
    Resamples an input array down to the given resolution 
    using the resampling method specified by statMethod (either RESAMPLE_AVERAGE
    or RESAMPLE_SUM).
    """
    if outSize[0] >= inarray.shape[0] or outSize[1] >= inarray.shape[1]:
        raise ValueError('Array can only be reduced in size')

    # assume same x and y and that it needs to be rounded up - 
    # appears to be issues with the rasterisation that need to be addressed...
    oldPixPerNewPix = int(np.ceil(inarray.shape[0] / outSize[0]))

    outArray = np.empty(outSize, inarray.dtype)

    # go through each new pixel
    for newy in range(outSize[0]):
        for newx in range(outSize[1]):
            oldx = newx * oldPixPerNewPix
            oldy = newy * oldPixPerNewPix
            total = 0.0
            for x in range(oldPixPerNewPix):
                for y in range(oldPixPerNewPix):
                    total += inarray[oldy + y, oldx + x]
            if statMethod == RESAMPLE_AVERAGE:
                outArray[newy, newx] = total / (oldPixPerNewPix * oldPixPerNewPix)
            else:
                outArray[newy, newx] = total / pixelRescale
    return outArray

def getRodentMaskForFile(rodentControlList, shpFile):
    """
    Helper function. Goes through rodentControlList and returns the
    mask for the file specified in shpPath
    """
    mask = None
    for testmask, startYear, revisit, controlIndicator, testShpFile in rodentControlList:
        if testShpFile == shpFile:
            mask = testmask
            break
    return mask

def populateResultArrays(loopIter, mastingMask, rodent_raster, 
        rodentExtentMask, rodentControlList, rodentAreaDictByMgmt, rodentDensity_2D,
        stoat_raster, stoatExtentMask, stoatSpatialDictByMgmt,stoatAreaDictByMgmt, 
        stoatDensity_2D, prey_raster, preyExtentMask, preySpatialDictByMgmt, preyAreaDictByMgmt, 
        preyDensity_2D, popAllYears_3D, year, year_all, keepYear):
    """
    ## Populate storage arrays
    """
    ## loop thru management areas (includes a key for total area)
    # NB: using sorted() keys to ensure the relation between i and key
    # is always consistent (might already be - not sure)
    for i, key in enumerate(sorted(preySpatialDictByMgmt.keys())):
        # for each control area, and for the entire area, calculate the mean proportion of
            # of prey_KMap that is in prey_raster, excluding kmap values == 0.

        ### (1) DO PREY
        mgmtMask = preySpatialDictByMgmt[key]               #mask prey cells in mgmt zone
        ### PREY DENSITY
        sppMgmtMask = mgmtMask & preyExtentMask            # prey habitat in mgmt zone
        sppDensity = np.sum(prey_raster[4,sppMgmtMask]) / preyAreaDictByMgmt[key]
#        sppDensity = np.sum(prey_raster[5,sppMgmtMask]) / preyAreaDictByMgmt[key]
        preyDensity_2D[i, year_all] = sppDensity        

        ### (2) DO STOAT DENSITY 
        mgmtMask = stoatSpatialDictByMgmt[key]               #mask stoat cells in mgmt zone
        sppMgmtMask = mgmtMask & stoatExtentMask            # stoat habitat in mgmt zone
        sppDensity = np.sum(stoat_raster[sppMgmtMask]) / stoatAreaDictByMgmt[key]
        stoatDensity_2D[i, year_all] = sppDensity        

        ### (3) DO RODENT DENSITY
        mgmtMask = getRodentMaskForFile(rodentControlList, key)             # mask rodent cells in mgmt
        sppMgmtMask = mgmtMask & rodentExtentMask           # rodent habitat in mgmt zone
        sppDensity = np.sum(rodent_raster[sppMgmtMask]) / rodentAreaDictByMgmt[key]
        rodentDensity_2D[i, year_all] = sppDensity        

    if (loopIter == 0) & keepYear:
        # 2) populate 3-D array on given iteration
        # copying of arrays not needed since we are inserting into
        # a given layer of popAllYears_3D - copy done implicitly
        popAllYears_3D['MastT'][year] = mastingMask
#        popAllYears_3D['ControlT'][year] = schedControlMask
        popAllYears_3D['rodentDensity'][year] = rodent_raster
        popAllYears_3D['stoatDensity'][year] = stoat_raster
        popAllYears_3D['preyDensity'][year] = prey_raster[4,:,:]







def populateResultControl(year, popAllYears_3D, controlThisYear):
    """
    ## POPULATE THE ANNUAL CONTROL RASTER FOR RESULTS
    """
    popAllYears_3D['ControlT'][year] = controlThisYear



def populateResultDensity(rodent_raster, rodentExtentMask, rodentControlList, 
        rodentAreaDictByMgmt, rodentDensity_2D_mth, stoat_raster, stoatExtentMask, 
        stoatSpatialDictByMgmt, stoatAreaDictByMgmt, stoatDensity_2D_mth, prey_raster, 
        preyExtentMask, preySpatialDictByMgmt, preyAreaDictByMgmt, preyDensity_2D_mth, tMth):
    """
    ## Populate storage arrays
    """
    ## loop thru management areas (includes a key for total area)
    # NB: using sorted() keys to ensure the relation between i and key
    # is always consistent (might already be - not sure)
    for i, key in enumerate(sorted(preySpatialDictByMgmt.keys())):
        # for each control area, and for the entire area, calculate the mean proportion of
            # of prey_KMap that is in prey_raster, excluding kmap values == 0.

        ### (1) DO PREY
        mgmtMask = preySpatialDictByMgmt[key]               #mask prey cells in mgmt zone
        ### PREY DENSITY
        sppMgmtMask = mgmtMask & preyExtentMask            # prey habitat in mgmt zone
        sppDensity = np.sum(prey_raster[4,sppMgmtMask]) / preyAreaDictByMgmt[key]
        preyDensity_2D_mth[i, tMth] = sppDensity        

        ### (2) DO STOAT DENSITY  
        mgmtMask = stoatSpatialDictByMgmt[key]               #mask stoat cells in mgmt zone
        sppMgmtMask = mgmtMask & stoatExtentMask            # stoat habitat in mgmt zone
        sppDensity = np.sum(stoat_raster[sppMgmtMask]) / stoatAreaDictByMgmt[key]
        stoatDensity_2D_mth[i, tMth] = sppDensity        

        ### (3) DO RODENT DENSITY
        mgmtMask = getRodentMaskForFile(rodentControlList, key)             # mask rodent cells in mgmt
        sppMgmtMask = mgmtMask & rodentExtentMask           # rodent habitat in mgmt zone
        sppDensity = np.sum(rodent_raster[sppMgmtMask]) / rodentAreaDictByMgmt[key]
        rodentDensity_2D_mth[i, tMth] = sppDensity        







def trackingTunnelRate(densityRodent, g0_TT, sigma_TT, 
        nights_TT, nTunnels, rodentResolution):
    """
    # Estimate the expected tracking tunnel rate in a management zone
    # Use rodent density at the rodent resolution, ie 200 m (4 ha)
    """
    # Area for single tracking tunnel
    TTArea = np.pi * ((4.0 * sigma_TT)**2)
    # density per m2
    den_m2 = densityRodent / (rodentResolution**2)         
    # number of rats in the tracking tunnel area (TTArea)
    nRats = np.int8(np.round((den_m2 * TTArea), 0))
    ## Eqn. 7
    # generate some random distances
    dist = rng.uniform(0.0, (4.0 * sigma_TT), nRats)
    # probability of detection by a single tracking tunnel over 1 night
    pdect = g0_TT * np.exp(-(dist**2) / 2.0 / (sigma_TT**2))
    ## Eqn. 6
    # probability of detection over multiple nights
    PD = 1.0 - np.prod((1.0 - pdect)**nights_TT)
    ## Eqn. 5
    trackingRate = rng.binomial(nTunnels, PD, 1) / nTunnels
    return(trackingRate)



def initialPopWrapper(sppProd, sppSurv, sppSurvDecay,
        sppRecDecay, sppInitialN, sppT, sppPresence, sppRasterShape, sppKMap):
    """
    ## calc stoat growth by pixel
    """
    mu = np.zeros(sppRasterShape)
    ## NUMBA LOOPING FUNCTION
    mu = growthLoop(mu, sppRasterShape, sppProd, 
        sppSurv, sppSurvDecay, sppRecDecay, sppInitialN, sppT, sppKMap, sppPresence)
    # Eqn 9: Add stochasticity
    spp_raster = rng.poisson(mu * 0.7, sppRasterShape)
    return(spp_raster)


@njit
def growthLoop(mu, sppRasterShape, growthRate, 
        sppSurv, sppSurvDecay, sppRecDecay,
        initialN, sppT, sppKMap, sppPresence):
    ## LOOP ROWS AND COLS
    for row in range(sppRasterShape[0]):
        for col in range(sppRasterShape[1]):
            N = initialN
            if ((sppPresence[row, col] == 0) | (sppKMap[row,col] == 0)):
                continue
            for t in range(sppT):
                ## Eqn. 10
                surv_i = sppSurv * (np.exp(-N**2 / 
                    sppKMap[row, col]**sppSurvDecay))
                NStar = N * surv_i
                ## Eqn 11
                pMaxRec = np.exp(-NStar**2 / 
                    sppKMap[row, col]**sppRecDecay)
                recRate = growthRate * pMaxRec
                ## Eqn 12
                N = (1 + recRate) * NStar
            mu[row,col] = N
    return(mu)

