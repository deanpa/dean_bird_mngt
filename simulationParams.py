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

        self.burnin = 40
        self.years = np.arange(40)

        ### SET ITERATIONS
        self.iter = 200
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
            os.makedirs(self.outputDataPath)
#        if not os.path.isdir(self.outputDataPath):
#            (baseDir / resultsPath).mkdir(parents = True)

        print('Results directory:', self.outputDataPath)
        print('############################')

        # ### SET DATA AND PATHS TO DIRECTORIES
        self.extentShp = os.path.join(self.inputDataPath, 'fullExtent.shp')
        self.AOIShp = os.path.join(self.inputDataPath, 'Kea_Model_Region3.shp')
        ##########################################
        ## TEST CONTROL ##########################
        # self.extentShp = os.path.join(self.inputDataPath, 'test_fullExtent.shp')
        # self.AOIShp = os.path.join(self.inputDataPath, 'test_AOI.shp')
        #self.extentShp = os.path.join(self.inputDataPath, 'extentDummy.shp')
        #self.AOIShp = os.path.join(self.inputDataPath, 'AOIDummy.shp')        
        # self.extentShp = os.path.join(self.inputDataPath, 'EglintonExtent.shp')
        # self.AOIShp = os.path.join(self.inputDataPath, 'EglintonAOI.shp')        
        # #########################################

        self.kClasses = os.path.join(self.inputDataPath, 'eco5_RAT.kea')
#        self.kClasses = os.path.join(self.inputDataPath, 'seed_Kea2.img')    
        # self.kClasses = os.path.join(self.inputDataPath, 'test_resource_Kea.grd')    
#        self.kClasses = os.path.join(self.inputDataPath, 'seed_KeaTemp.img')    
#        self.kClasses = os.path.join(self.inputDataPath, 'resource_Kea.img')    
        #self.kClasses = os.path.join(self.inputDataPath, 'resourcesDummyPureBeech.grd')    
        # self.kClasses = os.path.join(self.inputDataPath, 'resourceEglinton.tif')    
##        self.kClasses = os.path.join(self.inputDataPath, 'resourcesDummy.grd')    


        ### Area trapped in recent times.
        self.islands = os.path.join(self.inputDataPath, 'stoatTrapping750m.kea')
        self.DEM = os.path.join(self.inputDataPath, 'demKea200m.kea')
        self.preyHabitatShp = os.path.join(self.inputDataPath, 'Kea_Habitat.shp')
        # self.islands = os.path.join(self.inputDataPath, 'stoatTrappingRaster.img')
        # self.DEM = os.path.join(self.inputDataPath, 'dem200_kea.img')
        # self.preyHabitatShp = os.path.join(self.inputDataPath, 'Kea_Habitat.shp')
        # self.islands = os.path.join(self.inputDataPath, 'trapsDummy.tif')
        # self.DEM = os.path.join(self.inputDataPath, 'DEMDummy.tif')
        # self.preyHabitatShp = os.path.join(self.inputDataPath, 'KeaHabDummy.shp')
        # self.islands = os.path.join(self.inputDataPath, 'EglintonTraps.tif')
        # self.DEM = os.path.join(self.inputDataPath, 'EglintonDEM.tif')
        # self.preyHabitatShp = os.path.join(self.inputDataPath, 'EglintonKeaHabitat.shp')
##        self.islands = os.path.join(self.inputDataPath, 'trapsDummy.tif')
##        self.DEM = os.path.join(self.inputDataPath, 'DEMDummy.tif')
##        self.preyHabitatShp = os.path.join(self.inputDataPath, 'KeaHabDummy.shp')


        ##########################################
        ## TEST CONTROL ##########################
#        self.controlFile = os.path.join(self.inputDataPath, 'testControl.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'noCtrlTestControl.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'noCtrlAllMZsControl.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'oneoffTestControl.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, '3yrlyTestControl.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, '3yrlyMZ3Control.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, '3yrlyMZ4Control.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, 'reactiveTestControl.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, 'reactiveMZ3Control.csv') 
#        self.controlFile = os.path.join(self.inputDataPath, 'reactiveMZ4Control.csv') 
##        self.controlFile = os.path.join(self.inputDataPath, 'noCtrlDummy.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'oneOffCtrlDummy.csv') 
        #self.controlFile = os.path.join(self.inputDataPath, 'dblNovCtrlDummy.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, '3yrlyCtrlDummy.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'noCtrlEglinton.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'oneOffCtrlEglinton.csv') 
        # self.controlFile = os.path.join(self.inputDataPath, 'reactiveCtrlEglinton.csv') 
        ## "3yrlyCtrlDummy.csv" or "noCtrlDummy.csv" or "oneOffCtrlDumm.csv"
        # self.controlFile = os.path.join(self.inputDataPath, 'reactControl_kea1.csv')
    
        ### DEAN KEA SCENARIOS
#        self.controlFile = os.path.join(self.inputDataPath, 'control_kea1.csv') # control3 is effectively no control (st yr set to 100)
        self.controlFile = os.path.join(self.inputDataPath, 'control_keaSc2_3.csv') # control3 is effectively no control (st yr set to 100)

        ## LEAD POINT DATA
        self.leadPointData = os.path.join(self.inputDataPath, 'leadHutsVillages.csv')
        # self.leadPointData = os.path.join(self.inputDataPath, 'dummyLeadPoints.csv')
#        self.leadPointData = None

        ##Table for monthly resource/reodent K-values
        # self.seasAdjResFile = os.path.join(self.inputDataPath, 'testRes5.csv')
        self.seasAdjResFile = os.path.join(self.inputDataPath, 'testRes4.csv')

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

        ##Control parameters
        ##Mast-reactive control: proportion of zone in mast required for reactive control and mth control applied
        self.reactivePropMgmtMasting = 0.0  # 0.25    # 0.5  #0.5 # set > 0 to enable
        self.mastCtrlMth = self.monthDict['Nov']   ## 'Nov' is default
        ##Tracking Tunnel reactive control:
        self.threshold_TT = 0.2   #(1 = no reac) Thres prop of TT with detections
        self.reactiveAssessMth = self.monthDict['Jan']  #what month to do a mast prop or tracking tunnel assessment, mth7=Apr
                                                        ## DEFAULT IS 'Jan'
        self.reactiveCtrlDelay = 2  #delay in months from TT assessment to implementation of control
        #in preProcessing, getReactCtrlMth fn calcs reactiveCtrlMth and if jumps year line (Sept) due to delay
        ## Prescriptive Control: month prescriptive control is applied (can now be different to mast and TT- reactive control)
        self.prescrptCtrlMth = self.monthDict['Jul']    ## DEFAULT IS 'Jan'
        
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
        rodentFec= 8  #num offspring recruited per generation
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
        self.rodentRecDDcoef = 1 #0.2
        self.rodentSurv = 0.958 #per month 0.958
        self.rodentSurvDDcoef = 2 #2
        self.rodentTheta = 0.8  #1 gives Ricker model
        self.prpGrowRateControl = 1.0  # proportion of rodent growth before control is applied
        self.rodentProbEatBait = 0.98    #0.88  # 0.7 # pT
        self.pRodentPres = 0.9
        self.rodentInitialMultiplier = 0.8     
        self.rodentMaxAltitude = 1100.0  # metres based on Christie Mt Misery paper (inc slightly from orig 1000m)
        
        ##rat bounce parameters###
        self.rodentBouncePeriod = 36  #36 in months - time since control within which K/resources are multiplied to drive rat bounce
        self.rodentBounceMult = 1.5 #2.0  #2 how much to multiply resources/Kmap by to drive rat bounce, set to 1 to turn off rat bounce
        self.rodentBounceDecay = np.log(1/self.rodentBounceMult)/self.rodentBouncePeriod  #decay rate of bounce effect > mult goes to 1 at end of bounce period
        
        ######## TRACKING TUNNEL PARAMETERS
        self.g0_TT = 0.02                   # Tracking tunnel g0
        self.sigma_TT = 22.0                # Rat sigma
        self.nights_TT = 4                  # Tracking tunnel nights
        self.nTunnels = 120                 # number of tracking tunnels

        ######## STOAT PARAMETERS
        #self.stoatProd = 0.71
        stoatFec= 8.5   # 9.2  #num offspring recruited per generation
        self.stoatIRR = np.log(1+stoatFec/2) #convert to Instantaneous Rec Rate  
                        #Stoat breeding highly synchronised (daylength dependent) 
                        #so breed in one pulse/time step
        self.stoatSeasRec = np.array([0,0,0,0,1.,0,0,0,0,0,0,0]) 
                        #most births in Oct but dependent on mum until families  
                        #break up when young are 12-14 wks so say Jan recruitment
        #make sure dispersal happens after recuritment                        
        self.stoatSeasDisp = np.array([0,0,0,0,0,1,0,0,0,0,0,0], dtype=bool) #disperse Feb
        self.stoatRecDDcoef = 0.06  #0.08 #8
        self.stoatSurv = 0.94
        self.stoatSurvDDcoef = 0.095 #9.5
        self.stoatTheta = 1
        #self.stoatRecLag = 3 #calc recruitment based on rat numbers 3 mths before young stoats become in
        self.stoatPopSD = 0.22
        self.pEncToxic = 0.012   #.008 or .007          # operates at stoat scale
        self.pEatEncToxic = 0.95          # operates at stoat scale
        self.stoatInitialMultiplier = .85
        self.pStoatPres = 0.8
        self.initialStoatN = 2.0 #4
        self.stoatMaxAltitude = 1200.0 #changed from 1100 to 2000 based on Foster paper but is moot point in this model
                                        #cos rats are limited to <1100 m so can't get stoats at high elev 
                                        #like do in reality (stoats in alpine due to mice)
        
        ## PREY SPECIES PARAMETERS
        self.preyK = 20.0
        self.pPreyPres = 0.75
        self.initialpreyN = 5.0    #not used anymore?
        self.preyInitialMultiplier = 0.5  #Binomial(pPreyPres)*preyK*preyInitialMultiplier 
                                          #used to initialise kea densities (@t=0)
        #self.preyPsi = 0.7  # Eqn 32
        self.preyPsiStoat = 0.18  #  Effect of stoats on kea recruitment
        self.preyEtaStoat = np.array([0.11, 0.04, 0.04, 0.04]) #Effect of stoats on kea survival in each age class (0-3)
        self.preyPsiRodent = 0.0  #Effect of rodents on kea recruitment
        self.preyEtaRodent = np.array([0.0, 0.0, 0.0, 0.0]) #Effect of rodents on kea survival in each age class (0-3), none here
        
        self.rodentThresh = 0.5 #0.5 Threshold rat density per ha at which stoat prey switching kicks in
        self.stoatMult = 1.5 #3 Multiplier for stoat offtake of prey once prey switch kicks in

        #self.preySurv = np.array([0.982, 0.992, 0.994, 0.998]) #mthly max surv prey for age class 0-4
        self.preySurv = np.array([[0.982, 0.982],  #annual survival for each age class
                                 [0.992,  0.992],  #dim 0 = age class 0-3, dim 1 = nonmast, mast Surv
                                 [0.994,  0.994],  #here no difference between nonmast and mast years
                                 [0.998,  0.998]]) #adults potentially very high survival
        self.preySurvDDcoef = 110.0 #this is effectively a carrying capacity (per 1km2) for Kea survival, 
                                    #since so large here effectively no density dependence in surv
        #self.preyProd = 0.1067
        # preyFec= 2  #num chicks fledged per clutch
        # self.preyIRR = np.log(1+preyFec/2) #convert to annual Instantaneous Rec Rate  
        #                 #only have one clutch per season will try again later if lose first clutch 
        #                 #this prod is spread out across season by preySeasRec param                 
        self.preyFec = np.array([[0.0, 0.0],   #number of chicks fledged per breeding female per annum =2*0.5 (assumed even sex ratio)
                                 [0.0,  0.0],  #dim 0 = age class 0-3, dim 1 = nonmast, mast Fec 
                                 [0.0,  0.0],  #here no difference between nonmast and mast years
                                 [2.5, 2.5]])  #only adults >3 years breed with 2 female chicks per female per annum /2 to get per capita fec
        #self.preySeasRec = np.array([0,0,0,0.3,0.4,0.3,0,0,0,0,0,0]) 
        SeasRec = np.array([[0,   0.0], #how recruitment is spread out across season here peaks in Jan
                            [0.0, 0.0], #dim 0 = mth 0-11, dim 1 = nonmast, mast rec
                            [0.0, 0.0], #here no difference between nonmast and mast years
                            [0.3, 0.3], #egg laying pks Aug-Oct theefore recruitment peaks (+4mths=Dec-Feb)
                            [0.4, 0.4],
                            [0.3, 0.3],
                            [0,   0.0],
                            [0,   0],
                            [0,   0],
                            [0,   0],
                            [0,   0],
                            [0,   0]])         
        self.preySeasRec = SeasRec/np.sum(SeasRec,0)  #normalise to make sure adds to 1  
        self.preyPropBreedpa = np.array([1.0, 1.0])   #added in to account for higher prop breeding in mast year (here no diff)  
        #self.preyMastMultFec = 1 #multiplies up recruitment rate if masting is occuring #not needed for kea but will be for kaka
        self.preyRecDDcoef = 10.00  #this is effectively a carrying capacity (per 1km2) for Kea recruitment
        self.preyTheta = 2          #theta is how density dependence scales if <1 dd kicks in early, 
                                    #if >1 rate inc remains close to rm until K nearly reached
        self.preySeasDisp = np.array([0,0,0,0,0,0,1,0,0,0,0,0], dtype=bool) 
                        #dispersasl in autumn Mar after last recruitment 
                        #oh what age class does dispersal act on? <<<<chk this - should be juvs[0]
        self.preyInitAgeStr = np.array([0.3, 0.15, 0.15, 0.4], dtype=float)
#        self.preyInitAgeStr = np.array([0.3,0.1,0.1,0.1,0.4], dtype=float)
        self.preyMaxAltitude = 1400.0  # metres  dec from 2000 cos Kemp et al max elev breeding=1350m
                                #I guess that means creating an artificial refuge since kea above 1100
                                #are not exposed to predation (cos no rats above 1100m =>no stoats)
        ## PREY SIGMA FOR HOME RANGE STANDARD DEVIATION OF BIVARIATE NORMAL KERNEL
        self.preySigma = 5000
        self.pLeadMax = {'preAdult': 0.2, 'adult' : 0.15}  # .15 and 0.1
        self.removeLeadAction = False        

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
        self.annGrowthYears = 5


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



