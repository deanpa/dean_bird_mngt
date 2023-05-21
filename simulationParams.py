#!/usr/bin/env python

import os
import numpy as np
from pathlib import Path


class PreyParams(object):
    def __init__(self, species=None, scenario=None):
        """
        ## PARAMETER CLASS FOR SIMULATION
        """
        #########   SET THE SCENARIO    ######################
        self.species = species      # Species set in command line
        self.scenario = scenario    # 1,2,3,4 number set in command line
        ######################################################

        ### SET YEARS AND BURN IN YEARS

        self.burnin = 1
        self.years = np.arange(3)
        # self.burnin = 4
        # self.years = np.arange(6)

        ### SET ITERATIONS
        self.iter = 2
        ## IS FIRST RUN; IF FALSE IT WON'T RUN PREPROCESSING TO SAVE TIME
        self.firstRun = True        # True or False
        ## DO WE SUMMARISE RESULTS FOR FULL EXTENT? TRUE OR FALSE
        self.summariseFullExtent = False


        print('############################')
        print('Species:      ', self.species)
        print('Scenario:     ', self.scenario) 
        print('Iterations:   ', self.iter)
        print('Burnin:       ', self.burnin)
        print('Years:        ', len(self.years))
        print('First run:    ', self.firstRun)
        print('Sum Full ext: ', self.summariseFullExtent)

        ## GET BASE DIRECTORY LOCALLY OR ON NESI
        baseDir = os.getenv('BIRDPROJDIR', default='.')
        if baseDir == '.':
            baseDir = Path.cwd()
        # SET PATH TO DATA - SHOULD EXIST ALREADY 
        self.inputDataPath = os.path.join(baseDir, 'SpeciesProjects', self.species, 'Data')
        ## RESULTS DIRECTORY FOR THIS SCENARIO AND SPECIES
        scenDir = 'Scen' + str(self.scenario) + self.species
        resultsPath = os.path.join('SpeciesProjects', self.species, 'Results', scenDir)
        ## PUT TOGETHER THE BASE DIRECTORY AND PATH TO RESULTS DIRECTORY 
        self.outputDataPath = os.path.join(baseDir, resultsPath)
        ## MAKE NEW RESULTS DIRECTORY IF DOESN'T EXIST
        if not os.path.isdir(self.outputDataPath):
            (baseDir / resultsPath).mkdir(parents = True)

        print('Results directory:', self.outputDataPath)
        print('############################')

        # ### SET DATA AND PATHS TO DIRECTORIES
        self.extentShp = os.path.join(self.inputDataPath, 'fullExtent.shp')
        self.AOIShp = os.path.join(self.inputDataPath, 'Kea_Model_Region3.shp')
        ##########################################
        ## TEST CONTROL ##########################
#        # self.extentShp = os.path.join(self.inputDataPath, 'test_fullExtent.shp')
#        # self.AOIShp = os.path.join(self.inputDataPath, 'test_AOI.shp')
#        self.extentShp = os.path.join(self.inputDataPath, 'extentDummy.shp')
#        self.AOIShp = os.path.join(self.inputDataPath, 'AOIDummy.shp')        
        ##########################################
        ##########################################


#        self.kClasses = os.path.join(self.inputDataPath, 'seed_Kea2.img')    
        #self.kClasses = os.path.join(self.inputDataPath, 'resourcesDummyNewK.grd')    
        self.kClasses = os.path.join(self.inputDataPath, 'seed_KeaTemp.img')    
#        self.kClasses = os.path.join(self.inputDataPath, 'resourcesDummy.grd')    


        ### Area trapped in recent times.
        self.islands = os.path.join(self.inputDataPath, 'stoatTrappingRaster.img')
        self.DEM = os.path.join(self.inputDataPath, 'dem200_kea.img')
        self.preyHabitatShp = os.path.join(self.inputDataPath, 'Kea_Habitat.shp')
#        self.islands = os.path.join(self.inputDataPath, 'trapsDummy.tif')
#        self.DEM = os.path.join(self.inputDataPath, 'DEMDummy.tif')
#        self.preyHabitatShp = os.path.join(self.inputDataPath, 'KeaHabDummy.shp')

        ##########################################
        ## TEST CONTROL ##########################
#        self.controlFile = os.path.join(self.inputDataPath, 'testControl.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, 'noCtrlDummy.csv') 
        #self.controlFile = os.path.join(self.inputDataPath, 'oneOffCtrlDummy.csv') 
        #self.controlFile = os.path.join(self.inputDataPath, 'dblNovCtrlDummy.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, '3yrlyCtrlDummy.csv') 

        # "3yrlyCtrlDummy.csv" or "noCtrlDummy.csv" or "oneOffCtrlDumm.csv"
        ##########################################
        ##########################################
        self.controlFile = os.path.join(self.inputDataPath, 'reactControl_kea1.csv')
#        self.controlFile = os.path.join(self.inputDataPath, 'control_kea1.csv') # control3 is effectively no control (st yr set to 100)

        ## LEAD POINT DATA
        self.leadPointData = os.path.join(self.inputDataPath, 'leadPtsRegion3.csv')
        #self.leadPointData = os.path.join(self.inputDataPath, 'test_LeadPoints.csv')
#        self.leadPointData = os.path.join(self.inputDataPath, 'dummyLeadPoints.csv')
#        self.leadPointData = None


        self.seasAdjResFile = os.path.join(self.inputDataPath, 'mastLUpTable.csv')

        ## SET PICKLE FILE NAMES FOR PRE-PROCESSING AND RESULTS
        preProcFName = 'preProc_' + scenDir + '.pkl'
        self.preProcFName = os.path.join(self.outputDataPath, preProcFName)
        resultsFName = 'results_' + scenDir + '.pkl'
        self.resultsFName = os.path.join(self.outputDataPath, resultsFName)

        ## DATA DICTIONARY TO ASSOCIATE CALENDAR MONTHS WITH NUMBERS
        self.monthDict = {'Sep' : 0, 'Oct' : 1, 'Nov' : 2, 'Dec' : 3, 'Jan' : 4, 
            'Feb' : 5, 'Mar' : 6, 'Apr' : 7, 'May' : 8, 'Jun' : 9, 'Jul' : 10, 'Aug' : 11}

        ## RESOLUTION (RATS, STOATS, PREY)
        self.resolutions = (200.0, 1000.0, 1000.0)

        # Control parameters
        # proportion of zone in mast required for reactive control
        # model for control that is reactive to masting
        self.reactivePropMgmtMasting = 0   #0.5 # set > 0 to enable
        self.reactiveAssessMth = self.monthDict['Apr']  #what month to do a mast prop or tracking tunnel assessment, mth7=Apr
        self.reactiveCtrlDelay = 2  #delay implement reactive control, mth7+2=9=June 
        ## PRESCRIPTIVE CONTROL MONTH WILL BE SAME AS REACTIVE
        ## THIS IS CALCULATED IN PREPROCESSING (ASSESS MONTH PLUS DELAY)
        ## REACTIVE COULD JUMP THE YEAR LINE (SEPT) DUE TO DELAY, IN WHICH CASE,
        ## REACTIVE CONTROL COULD BE IN YEAR T+1, AND PRESCRIPTIVE IN YEAR T.
        ## THE REVISIT (IN CONTROL FILE) PARAMETER WOULD SET MIN TIME BETWEEN CONTROL.
#        self.prescrptCtrlMth = self.monthDict['Jun']    ## MONTH OF PRESCRIPTIVE CONTROL

        ### Masting parameters
        self.mastCellParams = (0.001, 1000.0)
        self.mastWindowSize = self.resolutions[0] * 150 # in metres
        self.mastRho = 16000.0
        self.mastPrEvent = 1.0 / 5.1  # p(mast) = 1 out of 5.1 years
        self.mastProportionParams = findBeta(0.5, 0.4)  # ALPHA AND BETA PARAMETERS
        self.mastSpatialSD = 2.1 



        ## RODENT PARAMETERS
        self.islandK = 2.0
        self.initialRodentN = 10
        #self.rodentProd = 0.3916
        rodentFec= 5.5  #num offspring recruited per generation
        self.rodentIRR = np.log(1+rodentFec/2)/4 #convert to monthly Instantaneous Rec Rate 
                        #divided by generation time (4 mths to sexual mturity) in this case
        #seasRec describes what proportion of population is breeding in each month
        # where month 0= Sept, month 11 = Aug BUT is actually time of recruitment - 
        # when indep young recruited into population
        #could also use as a proxy for age structure even if breeding season by setting
        # <1 because proportion of population is still juvenile and not breeding 
        self.rodentSeasRec = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5,0.5]) 
                        #can potentially breed year round - can reign this in...
        #Dispersl: bodge this for now boolean on or off (in future could have sine function rather than step): 
        #Maybe should just have dispersal in a pulse/only 1 month to save running dispersal algorithm multiple times?
        self.rodentSeasDisp = np.array([1,1,1,1,1,1,1,1,1,0,0,0], dtype=bool) #can disp most months of yr, don't in winter?? 
                            #Juv males tend to make up most of disp. popn. but we don't have age or sex structure so yeah...
        self.rodentSurv = 0.958 #per month 0.958
        self.rodentSurvDDcoef = 2
        self.rodentRecDDcoef = 0.1
        self.rodentTheta = 1  #1 gives Ricker model
        self.prpGrowRateControl = 1.0  # proportion of rodent growth before control is applied
        self.rodentProbEatBait = 0.7 # pT
        self.pRodentPres = 0.95
        self.rodentInitialMultiplier = 0.8     
        self.rodentMaxAltitude = 1000.0  # metres
        
        ##rat bounce parameters###
        self.rodentBouncePeriod = 28  #28 in months - time since control within which K/resources are multiplied to drive rat bounce
        self.rodentBounceMult = 2  #2 how much to multiply resources/Kmap by to drive rat bounce

        ######## TRACKING TUNNEL PARAMETERS
        self.threshold_TT = 1            #(1 = no reac) Thres prop of TT with detections
        self.g0_TT = 0.02                   # Tracking tunnel g0
        self.sigma_TT = 22.0                # Rat sigma
        self.nights_TT = 4                  # Tracking tunnel nights
        self.nTunnels = 120                 # number of tracking tunnels

        ######## STOAT PARAMETERS
        #self.stoatProd = 0.71
        stoatFec= 9.2  #num offspring recruited per generation
        self.stoatIRR = np.log(1+stoatFec/2) #convert to Instantaneous Rec Rate  
                        #Stoat breeding highly synchronised (daylength dependent) 
                        #so breed in one pulse/time step
        self.stoatSeasRec = np.array([0,0,0,0,1.,0,0,0,0,0,0,0]) 
                        #most births in Oct but dependent on mum until families  
                        #break up when young are 12-14 wks so say Jan recruitment
        #make sure dispersal happens after recuritment                        
        self.stoatSeasDisp = np.array([0,0,0,0,0,1,0,0,0,0,0,0], dtype=bool) #disperse Feb
        self.stoatRecDDcoef = 8
        self.stoatSurv = 0.94
        self.stoatSurvDDcoef = 9.5
        self.stoatTheta = 1
        #self.stoatRecLag = 3 #calc recruitment based on rat numbers 3 mths before young stoats become in
        self.stoatPopSD = 0.22
        self.pEncToxic = 0.004          # operates at stoat scale
        self.pEatEncToxic = 0.8          # operates at stoat scale
        self.stoatInitialMultiplier = .85
        self.pStoatPres = 0.75
        self.initialStoatN = 4.0
        self.stoatMaxAltitude = 1100.0
        
        ## PREY SPECIES PARAMETERS
        self.preyK = 20.0
        self.pPreyPres = 0.68
        self.initialpreyN = 5.0    #not used anymore?
        self.preyInitialMultiplier = 0.3  #Binomial(pPreyPres)*preyK*preyInitialMultiplier 
                                          #used to initialise kea densities (@t=0)
        #self.preyPsi = 0.7  # Eqn 32
        self.preyPsiStoat = 0.05  #  Effect of stoats on kea recruitment
        self.preyEtaStoatJuv = 0.08  #Effect of stoats on juvenile (age class 0) on kea survival
        self.preyEtaStoatAd = 0.02  #Effect of stoats on adult (age class 1-4) on kea survival
        self.preyPsiRodent = 0.0  #Effect of rodents on kea recruitment
        self.preyEtaRodentJuv = 0.0  #Effect of rodents on juvenile (age class 0) on kea survival
        self.preyEtaRodentAd = 0.0  #Effect of rodents on adult (age class 1-4) on kea survival
        self.competEffect = 0.0  #not used anymore?
        self.preyPopSD = .12     #not used anymore? 
        
        self.rodentThresh = 0.5 #0.5 Threshold rat density per ha at which stoat prey switching kicks in
        self.stoatMult = 3 #3 Multiplier for stoat offtake of prey once prey switch kicks in

        self.preySurv = np.array([0.982, 0.992, 0.994, 0.998]) #mthly max surv rats for age class 0-4
#        self.preySurv = np.array([0.982,0.992,0.994,0.997,0.998]) #mthly max surv rats for age class 0-4
        self.preySurvDDcoef = 110.0 #this is effectively a carrying capacity (per 1km2) for Kea survival, 
                                    #since so large here effectively no density dependence in surv
        #self.preyProd = 0.1067
        preyFec= 2  #num chicks fledged per clutch
        self.preyIRR = np.log(1+preyFec/2) #convert to annual Instantaneous Rec Rate  
                        #only have one clutch per season will try again later if lose first clutch 
                        #this prod is spread out across season by preySeasRec param                 
        self.preySeasRec = np.array([0,0,0,0.3,0.4,0.3,0,0,0,0,0,0]) 
                        #egg laying pks Aug-Oct theefore recruitment peaks (+4mths=
                        #incubation 22-24 days + 13 weeks in nest b4 fledge)in Dec-Feb
        self.preySeasDisp = np.array([0,0,0,0,0,0,1,0,0,0,0,0], dtype=bool) 
                        #dispersasl in autumn Mar after recruit.
                        #oh what age class does dispersal act on? <<<<chk this - should be juvs[0]
        self.preyMastMultFec = 1 #multiplies up recruitment rate if masting is occuring #not needed for kea but will be for kaka
        self.preyRecDDcoef = 10.00  #this is effectively a carrying capacity (per 1km2) for Kea recruitment
        self.preyTheta = 2          #theta is how density dependence scales if <1 dd kicks in early, 
                                    #if >1 rate inc remains close to rm until K nearly reached
        self.preyInitAgeStr = np.array([0.3, 0.15, 0.15, 0.4], dtype=float)
#        self.preyInitAgeStr = np.array([0.3,0.1,0.1,0.1,0.4], dtype=float)
        self.preyMaxAltitude = 2000.0  # metres
        ## PREY SIGMA FOR HOME RANGE STANDARD DEVIATION OF BIVARIATE NORMAL KERNEL
        self.preySigma = 5000
        self.pLeadMax = {'preAdult': 0.10, 'adult' : 0.07}
        

        ## IMMIGRATION AND EMIGRATION PARAMETERS
        self.gammaProbEmigrate = np.array([0.1, 0.2, 0.4])   # gamma for rodent, stoats, 
                                                    # Eqn 18, 20 and others 
        self.deltaImmigrate = np.array([0.8, 0.4, 0.05])    # delta (rodents, stoats, keas)
                                                        # Eqn 20 and others
        self.tauImmigrate = np.array([0.003, 0.25, 0.0])      # Eqn 18, 20 and others
                                                 # rate parameter Imm (rodent, stoat, kea)
        self.emigrationWindowSize = (self.resolutions[0] * 13, 
                    self.resolutions[1] * 20, self.resolutions[2] * 20) # in metres




        ## NUMBER OF YEARS OVER WHICH CALC PREY ANN GROWTH RATE
        self.annGrowthYears = 10


### data = preProcessing.preyData(pars)
### data.pickleSelf(os.path.join(pars.outputDataPath, 'preProcData.pkl'))

def findBeta(mu, sdev):
    """
    Find a and b of a Beta distribution given mean and standard deviation
    """
    sdevsq = sdev * sdev;
    a = mu * ( ( mu * (1.0 - mu) ) / sdevsq - 1.0)
    b = ( 1.0 - mu ) * ( ( mu * ( 1.0 - mu ) ) / sdevsq - 1.0)

    return a, b



## DIAGNOSTIC SIMULATION; plot and save graph
#simDiag = diagnosticSim.OnePixelSim(pars)



