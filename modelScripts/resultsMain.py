#!/usr/bin/env python

import os
import pickle
from modelScripts import resultsFX

def main(params):

    ## READ IN PRE-PROCESSING PICKLE
    fileobj = open(params.preProcFName, 'rb')
    data = pickle.load(fileobj)
    fileobj.close()

    ## READ IN PICKLE HOLDING DATA CREATED IN SIMULATION
    fileobj = open(params.resultsFName, 'rb')
    results = pickle.load(fileobj)
    fileobj.close()

    print('preProcess data:', params.preProcFName)    
    print('simulation results data:', params.resultsFName)    



#    resultsFX.writeTmpArray(params, data, results)




    ## RUN RESULTS PROCESSING MODULE
    results = resultsFX.processResults(params, data, results)


if __name__ == '__main__':
    main()























