import numpy as np

class KeaParams(object):
    """
    The parameters to the Kea model
    """
    def __init__(self):
        self.extentShp = None
        self.kClasses = None
        self.islands = None
        self.islandK = 1.3
        self.DEM = None
        # resolution for rodent, stoats, and kea, respectively.
        self.resolutions = (200.0, 1000.0, 2000.0)
        self.controlFile = None
        self.controlPathPrefix = None
        self.seasAdjResFile = None
        self.years = []
        self.burnin = 0
        self.iter = 3

        self.mastCellParams = (0.001, 1000.0)
        self.mastWindowSize = self.resolutions[0] * 150 # in metres
        self.mastRho = 16000.0
        self.mastPrEvent = 1.0 / 5.1  # p(mast) = 1 out of 5.1 years
        self.mastProportionParams = findBeta(0.5, 0.4)
        self.mastSpatialSD = 2.1  

        self.islandK = 2.0
        self.initialRodentN = 10
        self.rodentProd = 0.3916
        self.rodentSurv = 0.9707
        self.rodentSurvDDcoef = 3.8
        self.rodentRecDDcoef = 2
        self.rodentTheta = 0.6  #1 gives Ricker model
        #bodge this for now (in future could have sine function rather than step): 
        #seasRec describes what proportion of population
        #is breeding in each month where month 0= Sept, month 11 = Aug
        #actually is recruitment - when indep young recruited into population
        #could also use as a proxy for age structure even if breeding season set
        # <1 because proportion of population is still juvenile and not breeding 
        self.rodentSeasRec = [1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5,0.5] #can potentially breed year round - further controlled by DD fn
        self.rodentSeasDisp = np.array([1,1,1,1,1,1,1,1,1,0,0,0], dtype=bool) #can disperse most months of year

        self.prpGrowRateControl = 1.0  # proportion of rodent growth before control is applied
        self.rodentProbEatBait = 0.7 # pT

        self.pRodentPres = 0.95
        self.rodentInitialMultiplier = 0.8     
        
        
        ######## TRACKING TUNNEL PARAMETERS
        self.threshold_TT = 1.0             # Thres prop of TT with detections
#        self.prpGrowAssessRodent = 1.0      # prop of rat growth before apply control
        self.g0_TT = 0.02                   # Tracking tunnel g0
        self.sigma_TT = 22.0                # Rat sigma
        self.nights_TT = 4                  # Tracking tunnel nights
        self.nTunnels = 120                 # number of tracking tunnels


        self.gammaProbEmigrate = np.array([0.1, 0.2, 0.4])   # gamma for rodent, stoats, 
                                                    # Eqn 18, 20 and others 
        self.deltaImmigrate = np.array([0.8, 0.4, 0.05])    # delta (rodents, stoats, keas)         
                                                        # Eqn 20 and others
        self.tauImmigrate = np.array([0.003, 0.25, 0.0])      # Eqn 18, 20 and others
                                                 # rate parameter Imm (rodent, stoat, kea)
        self.emigrationWindowSize = (self.resolutions[0] * 13, 
                    self.resolutions[1] * 20, self.resolutions[2] * 20) # in metres

        self.rodentMaxAltitude = 1000.0  # metres
        self.keaMaxAltitude = 2000.0  # metres
        self.stoatMaxAltitude = 1100.0
        ############
        ############    Stoat parameters
    
        self.stoatProd = 0.7652
        self.stoatRecDDcoef = 5
        self.stoatSurv = 0.9672
        self.stoatSurvDDcoef = 6
        self.stoatTheta = 1
        self.stoatSeasRec = [0,0,0,1.,1.,1.,0,0,0,0,0,0] #born Sep-Nov, become indep Dec-Feb
        self.stoatSeasDisp = np.array([0,0,0,1,1,1,0,0,0,0,0,0], dtype=bool) #disperse Dec-Feb
        #self.stoatRecLag = 3 #calc recruitment based on rat numbers 3 mths before young stoats become independent
        self.stoatPopSD = 0.22
        self.pEncToxic = 0.004          # operates at stoat scale
        self.pEatEncToxic = 0.8          # operates at stoat scale
        self.stoatInitialMultiplier = .85
        self.pStoatPres = 0.75
        self.initialStoatN = 3.0
        
        ############
        ############    Kea parameters

        self.keaK = 20.0
        self.pKeaPres = 0.68
        self.initialkeaN = 5.0
        self.keaInitialMultiplier = 0.5
        self.keaPsi = 2.0  # Eqn 32
        self.competEffect = 0.0
        self.keaPopSD = .12

        self.keaSurv = [0.982,0.992,0.995,0.995,0.997]

        self.keaSurvDDcoef =120.0
        self.keaProd = 0.1067
        self.keaRecDDcoef =20
        self.keaTheta = 2
        self.keaSeasRec = [0,1.,1.,1.,1.,1.,1.,0,0,0,0,0] #kea become indep Oct-Mar
        self.keaSeasDisp = np.array([0,0,0,1,1,1,1,0,0,0,0,0], dtype=bool) #disperse Dec-Mar
        self.keaInitAgeStr = np.array([0.3,0.1,0.1,0.1,0.4], dtype=float)
        
        # model for control that is reactive to masting
        self.reactivePropMgmtMasting = 0.5 # set > 0 to enable
        self.reactiveAssessMth = 7  #what month to do a mast prop or tracking tunnel assessment, mth7=April
        self.reactiveCtrlMth = 11  #what month to implement reactive control, mth11=August 
        #nb: ctrl mth must be greater than assessment mth which is bit tricky when non-std year
        

        ## MODEL FOR CONTROL REACTIVE TO RODENT TRACKING TUNNEL RATE
        
        ## NUMBER OF YEARS OVER WHICH CALC KEA ANN GROWTH RATE
        self.annGrowthYears = 10


    def setExtentShapeFile(self, shpFile):
        self.extentShp = shpFile


    def setAOIShapeFile(self, shpFile):
        self.AOIShp = shpFile

    def setKClasses(self, kClasses):
        self.kClasses = kClasses

    def setIslands(self, islands):
        self.islands = islands

    def setIslandK(self, islandK):
        self.islandK = islandK

    def setDEM(self, dem):
        self.DEM = dem

    def setResolutions(self, resolutions):
        """
        Resolutions for rodent, stoats, and kea, respectively.
        """
        self.resolutions = resolutions

    def setControlFile(self, controlFile):
        self.controlFile = controlFile

    def setControlPathPrefix(self, prefix):
        self.controlPathPrefix = prefix
        
    def setSeasAdjResFile(self, resFile):
        self.seasAdjResFile = resFile        

    def setBurnin(self, burnin):
        self.burnin = burnin

    def setYears(self, years):
        self.years = years

    def setIterations(self, iter):
        self.iter = iter

    def setMastCellParams(self, a, b):
        self.mastCellParams = (a, b)

    def setMastWindowSize(self, winsize):
        "in metres"
        self.mastWindowSize = winsize

    def setMastRho(self, rho):
        self.mastRho = rho

    def setMastSpatialSD(self, mastSpatialSD):
        """
        ## Spatial masting
        """
        self.mastSpatialSD = mastSpatialSD 

    def setMastPrEvent(self, pr):
        self.mastPrEvent = pr

    def setMastProportionParams(self, mu, sdev):
        self.mastProportionParams = findBeta(mu, sdev)

    def setPRodentPresence(self, pRodentPres):
        self.pRodentPres = pRodentPres

    def setRodentInitialMultiplier(self, rodentInitialMultiplier):
        self.rodentInitialMultiplier = rodentInitialMultiplier

    def setRodentGrowthRate(self, rate):
        self.rodentGrowthRate = rate

    def setPrpGrowRateControl(self, prpGrowRateControl):
        self.prpGrowRateControl = prpGrowRateControl

    def setRodentProbEatBait(self, probEat):
        self.rodentProbEatBait = probEat

    def setRodentEmigrationWindowSize(self, winsize):
        "in metres"
        self.rodentEmigrationWindowSize = winsize

    ### NEW FUNCTIONS FOR RODENT GROWTH
    def setRodentT(self, rodentT):
        self.rodentT = rodentT

    def setRodentSurv(self, rodentSurv):
        self.rodentSurv = rodentSurv

    def setRodentSurvDDcoef(self, rodentSurvDDcoef):
        self.rodentSurvDDcoef = rodentSurvDDcoef         

    def setInitialRodentN(self, initialRodentN):
        self.initialRodentN = initialRodentN

    def setRodentProd(self, rodentProd):
        self.rodentProd = rodentProd

    def setRodentRecDDcoef(self, rodentDDcoef):
        self.rodentDDcoef = rodentDDcoef 
        
    def setRodentTheta(self, rodentTheta):
        self.rodentTheta = rodentTheta 
        
    def setRodentSeasRec(self, seasRec):
        self.rodentSeasRec= seasRec    
        
    def setRodentSeasDisp(self, seasDisp):
        self.rodentSeasDisp = np.array(seasDisp, dtype=bool) 


    def setThreshold_TT(self, threshold_TT):
        """
        Threshold proportion of tracking tunnels with detections
        """
        self.threshold_TT = threshold_TT

    def setG0_TT(self, g0_TT):
        """
        Tracking tunnel g0
        """
        self.g0_TT = g0_TT

    def setSigma_TT(self, sigma_TT):
        """
        rodent Sigma
        """
        self.sigma_TT = sigma_TT

    def setNights_TT(self, nights_TT):
        """
        Tracking tunnel nights
        """
        self.nights_TT = nights_TT

    def setNTunnels(self, nTunnels):
        """
        number of tracking tunnels
        """
        self.nTunnels = nTunnels


    ### STOAT FUNCTIONS AND PARAMETERS
    def setPStoatPresence(self, pStoatPres):
        self.pStoatPres = pStoatPres

    def setStoatInitialMultiplier(self, stoatInitialMultiplier):
        self.stoatInitialMultiplier = stoatInitialMultiplier

    def setStoatGrowthRate(self, stoatGrowthRate):
        self.stoatGrowthRate = stoatGrowthRate

    def setStoatPopSD(self, stoatPopSD):
        self.stoatPopSD = stoatPopSD

    def setStoatKAsymptotes(self, loAsym, hiAsym):
        self.kAsymptotes = (loAsym, hiAsym)

    def setPEncToxic(self, pEncToxic):
        self.pEncToxic = pEncToxic

    def setPEatEncToxic(self, pEatEncToxic):
        self.pEatEncToxic = pEatEncToxic

    ### NEW FUNCTIONS FOR STOAT GROWTH
    def setStoatT(self, stoatT):
        self.stoatT = stoatT

    def setStoatSurv(self, stoatSurv):
        self.stoatSurv = stoatSurv

    def setStoatSurvDDcoef(self, stoatSurvDDcoef):
        self.stoatSurvDDcoef = stoatSurvDDcoef

    def setInitialStoatN(self, initialStoatN):
        self.initialStoatN = initialStoatN

    def setStoatProd(self, stoatProd):
        self.stoatProd = stoatProd

    def setStoatRecDDcoef(self, stoatRecDDcoef):
        self.stoatRecDDcoef = stoatRecDDcoef 

    def setStoatTheta(self, stoatTheta):
        self.stoatTheta = stoatTheta 

    def setStoatSeasRec(self, seasRec):
        self.stoatSeasRec= seasRec         

    def setStoatSeasDisp(self, seasDisp):
        self.stoatSeasDisp = np.array(seasDisp, dtype=bool)



    ### KEA FUNCTIONS AND PARAMETERS
    def setPKeaPresence(self, pKeaPres):
        self.pKeaPres = pKeaPres

    def setKeaInitialMultiplier(self, keaInitialMultiplier):
        self.keaInitialMultiplier = keaInitialMultiplier

    def setKeaPsi(self, keaPsi):
        self.keaPsi = keaPsi

#    def setKeaGrowthRate(self, keaGrowthRate):
#        self.keaGrowthRate = keaGrowthRate

    def setKeaK(self, keaK):
        self.keaK = keaK

    def setKeaPopSD(self, keaPopSD):
        self.keaPopSD = keaPopSD

#    def setMinPredationSurv(self, minPredationSurv):
#        self.minPredationSurv = minPredationSurv

    def setCompetitionEffect(self, competEffect):
        self.competEffect = competEffect

    ### NEW FUNCTIONS FOR KEA GROWTH
    def setKeaSurv(self, keaSurv):
        self.keaSurv = np.array(keaSurv, dtype=float)    

    def setKeaSurvDDcoef(self, keaSurvDDcoef):
        self.keaSurvDDcoef = keaSurvDDcoef

    def setInitialKeaN(self, initialkeaN):
        self.initialkeaN = initialkeaN

    def setKeaProd(self, keaProd):
        self.keaProd = keaProd

    def setKeaRecDDcoef(self, keaRecDDcoef):
        self.keaRecDDcoef = keaRecDDcoef 

    def setKeaTheta(self, keaTheta):
        self.keaTheta = keaTheta

    def setKeaSeasRec(self, seasRec):
        self.keaSeasRec= seasRec       
        
    def setKeaSeasDisp(self, seasDisp):
        self.keaSeasDisp = np.array(seasDisp, dtype=bool)
        
    def setKeaInitAgeStr(self, initAgeStr):
        self.keaInitAgeStr = np.array(initAgeStr, dtype=float)    

    def setReactiveMode(self, prop):
        """
        Put into reactive mode by specifying the proportion
        of pixels in a control area that mast before the whole
        control area is controlled in a given year.
        """
        self.reactivePropMgmtMasting = prop
 
    def setReactiveAssessMth(self, assessMth):
        self.reactiveAssessMth = int(assessMth)
        
    def setReactiveCtrlMth(self, ctrlMth):
        self.reactiveCtrlMth = int(ctrlMth)


def findBeta(mu, sdev):
    """
    Find a and b of a Beta distribution given mean and standard deviation
    """
    sdevsq = sdev * sdev;
    a = mu * ( ( mu * (1.0 - mu) ) / sdevsq - 1.0)
    b = ( 1.0 - mu ) * ( ( mu * ( 1.0 - mu ) ) / sdevsq - 1.0)

    return a, b
