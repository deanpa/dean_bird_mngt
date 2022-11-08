#!/usr/bin/env python

import os
import numpy as np
#from modelScripts  import preProcessing

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
        self.years = np.arange(1)
        ### SET ITERATIONS
        self.iter = 1
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
        print('############################')


        # SET PATHS TO DATA AND RESULTS
        self.inputDataPath = os.path.join(os.getenv('KIWIPROJDIR', default='.'), 
                'SpeciesProjects', self.species, 'Data')
        scenDir = 'Scen' + str(self.scenario) + self.species
        self.outputDataPath = os.path.join(os.getenv('KIWIPROJDIR', default='.'), 
                'SpeciesProjects', self.species, 'Results', scenDir)

        # ### SET DATA AND PATHS TO DIRECTORIES
        self.extentShp = os.path.join(self.inputDataPath, 'fullExtent.shp')
        self.AOIShp = os.path.join(self.inputDataPath, 'Kea_Model_Region3.shp')
        self.kClasses = os.path.join(self.inputDataPath, 'seed_Kea.img')    

        ### Area trapped in recent times.
        self.islands = os.path.join(self.inputDataPath, 'stoatTrappingRaster.img')
        self.DEM = os.path.join(self.inputDataPath, 'dem200_kea.img')
        self.preyHabitatShp = os.path.join(self.inputDataPath, 'Kea_Habitat.shp')
        self.controlFile = os.path.join(self.inputDataPath, 'control_kea1.csv') # control3 is effectively no control (st yr set to 100)
        self.seasAdjResFile = os.path.join(self.inputDataPath, 'mastLUpTable.csv')
#        ### SET DATA AND PATHS TO DIRECTORIES - dummy/simplified landscape
#        self.extentShp = os.path.join(self.inputDataPath, 'fullExtent.shp')
#        self.AOIShp = os.path.join(self.inputDataPath, 'Kea_Model_Region3.shp')
#        self.KClasses = os.path.join(self.inputDataPath, 'kClassesDummy.grd')       #'seeds_RmIslands.img'))    
#        ### Area trapped in recent times.
#        self.islands = os.path.join(self.inputDataPath, 'trapsDummy.tif')
#        self.DEM = os.path.join(self.inputDataPath, 'dem_region3.tif')
#        self.resolutions = (200.0, 1000.0, 2000.0)
#        self.controlFile = os.path.join(self.inputDataPath, 'controlDummy.csv') # control3 is effectively no control (st yr set to 100)
#        self.seasAdjResFile = os.path.join(self.inputDataPath, 'mastLUpTable.csv')

        ## RESOLUTION (RATS, STOATS, PREY)
        self.resolutions = (200.0, 1000.0, 1000.0)

        # Control parameters
        # proportion of zone in mast required for reactive control
        # model for control that is reactive to masting
        self.reactivePropMgmtMasting = 0.5 # set > 0 to enable
        self.reactiveAssessMth = 7  #what month to do a mast prop or tracking tunnel assessment, mth7=Apr
        self.reactiveCtrlMth = 11  #what month to implement reactive control, mth11=August 
        #nb: ctrl mth must be greater than assessment mth which is bit tricky when non-std year
        
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
        self.rodentProd = 0.3916
        self.rodentSurv = 0.9707
        self.rodentSurvDDcoef = 3.0
        self.rodentRecDDcoef = 1.8
        self.rodentTheta = 0.6  #1 gives Ricker model
        #bodge this for now (in future could have sine function rather than step): 
        #seasRec describes what proportion of population
        #is breeding in each month where month 0= Sept, month 11 = Aug
        #actually is recruitment - when indep young recruited into population
        #could also use as a proxy for age structure even if breeding season set
        # <1 because proportion of population is still juvenile and not breeding 
        self.rodentSeasRec = [1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5,0.5] #can potentially breed year round -
        self.rodentSeasDisp = np.array([1,1,1,1,1,1,1,1,1,0,0,0], dtype=bool) #can disperse most months o
        self.prpGrowRateControl = 1.0  # proportion of rodent growth before control is applied
        self.rodentProbEatBait = 0.7 # pT
        self.pRodentPres = 0.95
        self.rodentInitialMultiplier = 0.8     
        self.rodentMaxAltitude = 1000.0  # metres

        ######## TRACKING TUNNEL PARAMETERS
        self.threshold_TT = 1.0             # Thres prop of TT with detections
        self.g0_TT = 0.02                   # Tracking tunnel g0
        self.sigma_TT = 22.0                # Rat sigma
        self.nights_TT = 4                  # Tracking tunnel nights
        self.nTunnels = 120                 # number of tracking tunnels

        ######## STOAT PARAMETERS
        self.stoatProd = 0.71
        self.stoatRecDDcoef = 8
        self.stoatSurv = 0.9439
        self.stoatSurvDDcoef = 10
        self.stoatTheta = 1
        self.stoatSeasRec = [0,0,0,1.,1.,1.,0,0,0,0,0,0] #born Sep-Nov, become indep Dec-Feb
        self.stoatSeasDisp = np.array([0,0,0,1,1,1,0,0,0,0,0,0], dtype=bool) #disperse Dec-Feb
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
        self.initialpreyN = 5.0
        self.preyInitialMultiplier = 0.3
        self.preyPsi = 0.7  # Eqn 32
        self.competEffect = 0.0
        self.preyPopSD = .12

        self.preySurv = np.array([0.982,0.992,0.995,0.995,0.997])

        self.preySurvDDcoef = 110.0
        self.preyProd = 0.1067
        self.preyRecDDcoef = 10.00
        self.preyTheta = 2
        self.preySeasRec = np.array([0,0,0.5,0.8,1.,0.8,0.5,0.3,0.1,0,0,0]) #egg laying July-Jan with peak in Sept therefore peak recruitment (+4mths) in Jan
        self.preySeasDisp = np.array([0,0,0,1,1,1,1,0,0,0,0,0], dtype=bool) #disperse Dec-Mar
        self.preyInitAgeStr = np.array([0.3,0.1,0.1,0.1,0.4], dtype=float)
        self.preyMaxAltitude = 2000.0  # metres


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



