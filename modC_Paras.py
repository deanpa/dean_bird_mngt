#!/usr/bin/env python

import os
import numpy as np
from kiwimodelMthlyTS import params
from kiwimodelMthlyTS  import preProcessing
#from kiwimodel import diagnosticSim

pars = params.KiwiParams()

# set paths to scripts and data
pars.inputDataPath = os.path.join(os.getenv('KIWIPROJDIR', default='.'), 'kiwi_data')
pars.outputDataPath = os.path.join(os.getenv('KIWIPROJDIR', default='.'), 
                'KiwiProjResults', 'modC_Results')
if not os.path.isdir(pars.outputDataPath):
    os.mkdir(pars.outputDataPath)

### SET DATA AND PATHS TO DIRECTORIES
pars.setExtentShapeFile(os.path.join(pars.inputDataPath, 'FiordlandConArea.shp'))
pars.setKClasses(os.path.join(pars.inputDataPath, 'seeds_RmIs_EqualBeech.img'))       #'seeds_RmIslands.img'))    
### Area trapped in recent times.
pars.setIslands(os.path.join(pars.inputDataPath, 'EmsTraps.tif'))
pars.setDEM(os.path.join(pars.inputDataPath, 'dem.tif'))
pars.setResolutions((200.0, 1000.0, 1000.0))
pars.setControlFile(os.path.join(pars.inputDataPath, 'control8.csv')) # control3 is effectively no control (st yr set to 100)
pars.setControlPathPrefix(pars.inputDataPath)
pars.setSeasAdjResFile(os.path.join(pars.inputDataPath, 'mastLUpTable.csv'))

### SET YEARS AND BURN IN YEARS
pars.setBurnin(10)
pars.setYears(np.arange(20))
### SET ITERATIONS
pars.setIterations(10)
print('iterations:', pars.iter)
print('Burnin:', pars.burnin)
print('Years:', len(pars.years))

# Control parameters
# proportion of zone in mast required for reactive control
pars.setReactiveMode(0.0)   #if = 0.0, then no reactive control to masting
pars.setReactiveAssessMth(0)
pars.setReactiveCtrlMth(2)
pars.setThreshold_TT(1.0) #set to 1.0 to have no reactive control to rat tracking rate

### Masting parameters
pars.setMastRho(16000.0)
#pars.setMastWindowSize(130)
pars.setMastCellParams(0.001, 1000.0)
#pars.setMastProportionParams(0.5, 0.4)
pars.setMastSpatialSD(2.1)
pars.setMastPrEvent(1.0 / 5.1)          #5.3)

## rodent parameters
pars.setPRodentPresence(0.95)
pars.setRodentInitialMultiplier(.80)
pars.setRodentProbEatBait(0.7)
###pars.setRodentGrowthRate(1.2)
pars.setPrpGrowRateControl(0.25) # 1 if before growth


## NEW RODENT PARAMETERS
#pars.setRodentT(10)     # 10 time years to K
pars.setRodentSurv(0.9707)
pars.setRodentSurvDDcoef(3.0)    
pars.setInitialRodentN(10)   
pars.setRodentProd(0.3916)   #2.5
pars.setRodentRecDDcoef(1.8)   
pars.setRodentTheta(0.6)    
pars.setRodentSeasRec([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])
pars.setRodentSeasDisp([1,1,1,1,1,1,1,1,1,0,0,0]) #can disperse most months of year


## rodent tracking tunnel parameters
pars.setG0_TT(0.02)
pars.setSigma_TT(22.0)
pars.setNights_TT(4)
pars.setNTunnels(120)

## stoat parameters
pars.setPStoatPresence(0.75)
#pars.setStoatInitialMultiplier(3)
pars.setIslandK(2.0)                    # rodents per ha
#pars.setStoatGrowthRate(1.2)
#pars.setStoatKAsymptotes(1.0, 8.0)
pars.setPEncToxic(0.004)  #.006        # operates at stoat scale
pars.setPEatEncToxic(0.8)
pars.setStoatPopSD(0.22)


## NEW STOAT PARAMETERS
#pars.setStoatT(5)     
pars.setStoatSurv(0.9439)
pars.setStoatSurvDDcoef(10)
pars.setInitialStoatN(4.0)   
pars.setStoatProd(0.7100)
pars.setStoatRecDDcoef(8)
pars.setStoatTheta(1)
pars.setStoatSeasRec([0,0,0,1.,1.,1.,0,0,0,0,0,0])
pars.setStoatSeasDisp([0,0,0,1,1,1,0,0,0,0,0,0]) #disperse Dec-Feb

### Kiwi parameters
#pars.setKiwiK(20)
pars.setPKiwiPresence(0.68)
pars.setKiwiInitialMultiplier(0.3)
pars.setKiwiPsi(0.7)  
pars.setCompetitionEffect(0.000)    #0.004    
#pars.setKiwiGrowthRate(.10)
pars.setKiwiPopSD(.12)



## NEW KIWI PARAMETERS
pars.setKiwiSurv(0.9913)
pars.setKiwiSurvDDcoef(110)
pars.setInitialKiwiN(5.0)   
pars.setKiwiProd(0.1067)  #0.6
pars.setKiwiRecDDcoef(10.0)
pars.setKiwiTheta(2)
pars.setKiwiSeasRec([0,1.,1.,1.,1.,1.,1.,0,0,0,0,0])
pars.setKiwiSeasDisp([0,0,0,1,1,1,1,0,0,0,0,0]) #disperse Dec-Mar

data = preProcessing.KiwiData(pars)
data.pickleSelf(os.path.join(pars.outputDataPath, 'preProcData.pkl'))



## DIAGNOSTIC SIMULATION; plot and save graph
#simDiag = diagnosticSim.OnePixelSim(pars)



