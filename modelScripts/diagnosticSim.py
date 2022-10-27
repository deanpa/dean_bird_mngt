import numpy as np
import pylab as P


class OnePixelSim(object):
    """
    Simulate a single pixel at resol of keas for diagnosis/evaluation of parameters
    """
    def __init__(self, params):
        self.params = params
        self.nYears = 100
        self.rodent_K = np.array([5.0, 20.0])
        self.controlFreq = 5
#        self.reactiveControl = False
        self.probImm = 0.5
        self.controlType = 'Regular'   # Options are 'Both', 'Regular', 'Reactive'


        ########################
        ## RUN FUNCTIONS
        self.setInitialPars()
        self.simYears()
        self.plotResults()
        ## END CALLING FUNCTIONS
        #########################

    def setInitialPars(self):
        """
        set initial pop for rodents, stoats and keas
        """
        ## RODENTS
        self.nHectInRodent = (self.params.resolutions[0] / 100.0)**2
        self.nRodentPixInStoat = (self.params.resolutions[1] / self.params.resolutions[0])**2
        ## UPDATE RODENT K - for set resol in params
        self.rodent_K = self.rodent_K * self.nHectInRodent  
        self.rodent = np.random.poisson(self.rodent_K[0] * self.params.rodentInitialMultiplier)
        ## STOATS
        self.getStoatK()
        # Eqn 24    Random stoat population 
        self.stoat = np.random.poisson(self.stoat_K)
#        self.stoat = np.random.poisson(self.stoat_K * self.params.stoatInitialMultiplier)
        ## KEAS
        ## Get initial kea population
        self.kea = np.random.poisson(self.params.keaInitialMultiplier * self.params.keaK)
        ## YEARS, CONTROL, AND MAST ARRAYS
        self.popArray = np.zeros((3, self.nYears))
        if self.controlType == 'Reactive':
            self.controlArray = np.zeros(self.nYears, dtype = bool)
        else:
            controlYearsTmp = np.arange(self.controlFreq - 1, self.nYears, self.controlFreq)
            self.controlArray = np.in1d(np.arange(self.nYears), controlYearsTmp)

###        if self.reactiveControl:
###            self.controlArray = np.zeros(self.nYears, dtype = bool)
###        else:
###            controlYearsTmp = np.arange(self.controlFreq - 1, self.nYears, self.controlFreq)
###            self.controlArray = np.in1d(np.arange(self.nYears), controlYearsTmp)
        self.mastArray = np.zeros(self.nYears)

        print('self.nRodentPixels', self.nRodentPixInStoat, 'ha in rod', self.nHectInRodent)
        print('Initial rodent', self.rodent, 'rodent per ha', self.rodent / 4.0, 'stoat K', self.stoat_K, 
            'stoat', self.stoat, 'kea', self.kea)

    def getStoatK(self):
        """
        use rodents to gt stoat K; generalised logistic function
        """
        # Eqn 17
        rodentPerHa = self.rodent / self.nHectInRodent
        stoatNumerator = self.params.kAsymptotes[1] - self.params.kAsymptotes[0]
        denomArr = (1.0 + self.params.Q * 
            np.exp(-self.params.kOmega*(rodentPerHa - self.params.L)))**(1.0/self.params.kNu)
        self.stoat_K = self.params.kAsymptotes[0] + (stoatNumerator / denomArr)

    def simYears(self):
        """
        simulate dynamics over all years
        """
        mast_T2 = False
        for i in range(self.nYears):
            ## MASTING
            mast_T1 = np.bool(np.random.binomial(1, self.params.mastPrEvent))
            if mast_T1:
                self.mastArray[i] = .33
            else:
                self.mastArray[i] = -10.
            ## CONTROL
            if (self.controlType != 'Regular') & mast_T1:
                self.controlArray[i] = True
            if self.controlArray[i]:
                self.doControl()




#            if self.reactiveControl & mast_T1:
#                self.controlArray[i] = True
#            if self.controlArray[i]:
#                self.doControl()



            ## RODENT GROWTH
            self.doRodentGrowth(mast_T1, i)
            
#            print('#######################     control', self.controlArray[i])


            ## STOAT GROWTH
            self.doStoatGrowth(mast_T1, mast_T2, i)
#            print('i', i, 'mast_T1', mast_T1, 'mast_T2', mast_T2, 'rodent_N', self.rodent, 
#                    'stoat_K', self.stoat_K, 'stoat N', self.stoat)
            ## KEA GROWTH
            self.doKeaGrowth(i)

            ## UPDATE MASTING HISTORY
            mast_T2 = mast_T1


    def doControl(self):
        """
        ## calc number of rodents eating toxin and dieing
        """
        # Eqn 8 - done at rodent resolution
        nToxicRodents = np.random.binomial(self.rodent, self.params.rodentProbEatBait)

#        print('ntoxic', nToxicRodents)

        ## update rodent pop. with mortality
        # Eqn 9
        self.rodent -= nToxicRodents
        ## update stoat_raster by consuming toxic rodents
        # Eqn 23        # probability of encounter
        nToxicAtStoatResol = nToxicRodents * self.nRodentPixInStoat 
        pEnc = 1.0 - np.exp(-self.params.pEncToxic * nToxicAtStoatResol)
        # probability individ. stoat eating a toxic rodent
        pEat = self.params.pEatEncToxic * pEnc
        print('stoat before eat', self.stoat, 'pEat', pEat, 'nToxicRodents', nToxicRodents)
        # Eqn 20 and 22     # update stoat_raster
        self.stoat = np.int(np.round(self.stoat * (1.0 - pEat), 0))
        print('stoat after eat', self.stoat, 'nToxicAtStoatResol', nToxicAtStoatResol)



    def doRodentGrowth(self, mast_T1, i):
        """
        ## calc rodent growth and new pop size at rodent resol
        """
        ## IF RODENT = 0; ALLOW FOR IMMIGRATION
        if self.rodent == 0:
            self.rodent = np.random.binomial(5, .5)

        ## if masting in t-1, update kmap
        if mast_T1:
            kRodent = self.rodent_K[1]
        else:
            kRodent = self.rodent_K[0]
        # Eqn 11
#        mu = (self.rodent * np.exp(self.params.rodentGrowthRate) *  
#            np.exp(-self.rodentDampen * self.rodent / kRodent)) 
        mu = (self.rodent * np.exp(self.params.rodentGrowthRate * (1.0 - 
            (self.rodent / kRodent)))) 
        # Eqn 10: Add stochasticity
        self.rodent = np.random.poisson(mu)
#        print('###############')
#        print('mast_T1', mast_T1, 'rodent_k', self.rodent_K, 'rodent', self.rodent)

        self.popArray[0, i] = self.rodent / self.nHectInRodent


    def doStoatGrowth(self, mast_T1, mast_T2, i):
        """
        let stoat pop growth according to K
        """
        ## IF STOAT = 0; ALLOW FOR IMMIGRATION
        if self.stoat == 0:
            self.stoat = np.random.binomial(1, .3)
#        print('pre growth stoat', self.stoat)
            
        ### Stoats decline more slowly than rodents following a mast,
        ### Calc ave of current stoat_kMap and stoat_kMapT_1
        if mast_T2:
            stoatK_T2 = self.stoat_K
            self.getStoatK()
            self.stoat_K = (self.stoat_K + stoatK_T2) / 2.0
        ### if no mast T-2 years
        else:
            self.getStoatK()
        # Eqn 25 stoat population growth 
#        mu = (self.stoat * np.exp(self.params.stoatGrowthRate) * 
#            np.exp(-self.stoatDampen * self.stoat / self.stoat_K)) 

        mu = self.stoat * np.exp(self.params.stoatGrowthRate * (1.0 - 
            (self.stoat / self.stoat_K))) 

        # Eqn 24
        self.stoat = np.random.poisson(mu)

#        print('stoat K', self.stoat_K, 'stoat', self.stoat)
        self.popArray[1, i] = self.stoat


    def doKeaGrowth(self, i):
        """
        ## KEA POPULATION GROWS
        """
        ## IF KEA = 0; ALLOW FOR IMMIGRATION
        if self.kea == 0:
            self.kea = np.random.binomial(2, .3)

        # Eqn 32
###        r_kt = (self.params.keaGrowthRange[0] + np.exp(-self.params.keaPsi * self.stoat) * 
###                    (self.params.keaGrowthRange[1] - self.params.keaGrowthRange[0]))
        s = np.exp(-self.params.keaPsi * self.stoat)

        # Eqn 34
        mu_kt = (self.kea * np.exp(0.1 * (1.0 - (self.kea / self.params.keaK))) * s)
#        mu_kt = (self.kea * np.exp(r_kt * (1.0 - (self.kea / self.params.keaK))))
#        mu_kt = (self.kea * np.exp(r_kt) * 
#                np.exp(-self.keaDampen * self.kea / self.params.keaK))
        
#        print('Before Kea growth', self.kea, 'stoat', self.stoat, 'r', r_kt, 'mu', mu_kt)

        # Eqn 31
        self.kea = np.random.poisson(mu_kt)
#        print('i', i, 's', s, 'mu_kt', mu_kt, 'After Kea growth', self.kea)
        self.popArray[2, i] = self.kea

    def plotResults(self):
        """
        plot results
        """
        years = np.arange(1, self.nYears+1)
        maxy = np.max(self.popArray) + 5


        P.figure(figsize=(11,9))

        controlpoints = np.where(self.controlArray, maxy - 6., -10.)

        ax1 = P.gca()

        lns1 = ax1.plot(years, self.popArray[0], color = 'k', linewidth = 3.0, 
                label = 'Rats density $(individuals * ha^{-1})$')
        lns1 = tuple(lns1)        
        lns2 = ax1.plot(years, self.popArray[1], color = 'r', linewidth = 3.0, 
                label = 'Stoats density $(individuals * km^{-2})$')        
        lns2 = tuple(lns2)        
        lns3 = ax1.plot(years, self.popArray[2], color = 'b', linewidth = 3.0, 
                label = 'Keas density $(individuals * km^{-2})$')
        lns3 = tuple(lns3)        
        lns4 = ax1.plot(years, self.mastArray , 'gD',  
            markerfacecolor = 'g', ms = 10, mew = 1.5)
        lns4 = tuple(lns4)        
        lns5 = ax1.plot(years, controlpoints , 'r*',  
            markerfacecolor = 'r', ms = 15, mew = 2.0)
        lns5 = tuple(lns5)
        ax1.legend([lns1, lns2, lns3, lns4, lns5],
            ['Rats density $(individuals * ha^{-1})$',
            'Stoats density $(individuals * km^{-2})$',
            'Keas density $(individuals * km^{-2})$',
            'Masting event',
            '1080 control applied'],
            loc = 'upper left')        


#        lns = lns1 + lns2 + lns3 
#        labs = [l.get_label() for l in lns]
#        ax1.legend(lns, labs, loc = 'upper left', fontsize = 12)



#        maxy = np.max(self.rodent)
        miny = -.33
        ax1.set_ylim([miny, maxy])

        P.xlabel('Years', fontsize = 17)
        P.ylabel('Local population density', fontsize = 17)


        P.show()
