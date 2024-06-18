#!/usr/bin/env python

import os
import numpy as np
from scipy import stats
import pylab as P
from scipy.special import gamma
from scipy.special import gammaln
from scipy.stats.mstats import mquantiles
from numba import njit
import datetime



class survTest(object):
    def __init__(self):
        self.pSurv = {'Rats' : 0.68, 'Stoats' : 0.9, 'Prey' : 0.98}


class leadTest(object):
    def __init__(self):
        self.preySigma = 5000
        self.pLeadMax = {'preAdult': 0.22, 'adult' : 0.18}        
        self.preySurv = np.array([0.98, 0.98, 0.982, 0.988])**12
        self.nYears = 25
        self.years = np.arange(1, self.nYears + 1)
        self.plotPSurv()
        self.plotDistDecay()

    def plotPSurv(self):
        P.figure(figsize=(8,8))
        pS = np.repeat(self.preySurv[-1], self.nYears)
        pS[:3] = self.preySurv[:3]
        pS = np.cumprod(pS)
        print('cumprod', pS, (1-pS))
        P.plot(self.years, pS, color='k', linewidth = 4, label = 'No lead risk')
        pSLead = np.repeat(1 - self.pLeadMax['adult'], self.nYears)
        pSLead[:3] = 1 - self.pLeadMax['preAdult']
        pSLead[:3] = pSLead[:3] * self.preySurv[:3]
        pSLead[3:] = pSLead[3:] * self.preySurv[-1]
        pSLead = np.cumprod(pSLead)
        P.plot(self.years, pSLead, color='r', linewidth = 4, label = 'Lead risk')
        P.xlabel('Age (years)', fontsize = 14)
        P.ylabel('Probability of survival', fontsize = 14)
        P.legend(loc = 'upper right')
        P.savefig('psurviveLead.png', dpi = 150)
        P.show()

    def plotDistDecay(self):
        dist = np.arange(0, self.preySigma * 3.0, 1)
        pS0 = (1-self.pLeadMax['preAdult']) * self.preySurv[0]
        pS = pS0 * np.exp(-dist**2 / 2.0 / self.preySigma**2)
        P.figure(figsize=(8,8))
        P.plot(dist, pS, 'k', linewidth = 4)
        P.xlabel('Distance from lead point to HR centre (m)', fontsize = 14)
        P.ylabel('First year probability of survival', fontsize = 14)
        P.savefig('leadDecay.png', dpi = 150)
        P.show()

        



########            Main function
#######
def main():

#    survTest()
    leadTest()

if __name__ == '__main__':
    main()

