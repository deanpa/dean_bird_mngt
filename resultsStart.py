#!/usr/bin/env python

import os
import sys
import importlib
import optparse
from modelScripts import resultsMain

class CmdArgs(object):
    def __init__(self):
        p = optparse.OptionParser()
        p.add_option("--species", dest="species", help="Set species")
        p.add_option("--scenario", dest="scenario", help="Set scenario")
        (options, args) = p.parse_args()
        self.__dict__.update(options.__dict__)

if __name__ == '__main__':
    cmdargs = CmdArgs()
    print('Species results start', cmdargs.species)
    print('Scenario results start', int(cmdargs.scenario), type(int(cmdargs.scenario)))

    #####################################################
    ######################
    ##  RESULT DIRECTORY    -    USER DIRECTORY FOR RESULTS (SPECIES AND SCENARIO)
    scen = cmdargs.scenario      #'1'
    species = cmdargs.species                      #'Cats'
    ######################
    #####################################################

    ## MAKE PATH TO RESULTS
    resScenPath = os.path.join('SpeciesProjects', species, 'Results', 'Scen' + scen + species)
    paramsMod = 'params_' + species + 'Scen' + scen
    print('paramsMod', paramsMod)

    ## APPEND SYS DIRECTORY TO IMPORT PARAMS MODULE
    basepath = os.getenv('KIWIPROJDIR', default = '.')
    if basepath == '.':
        pathDir = os.path.join(os.getcwd(), resScenPath)
    else:
        pathDir = os.path.join(basepath, resScenPath)
    sys.path.append(pathDir)

    print('path to results', pathDir)

    ##  IMPORT PARAMS MODULE   -   USER MODIFY MODULE NAME
    scenParams = importlib.import_module(paramsMod)
    params = scenParams.PreyParams(cmdargs.species, int(cmdargs.scenario))

    resultsMain.main(params)

##########  COMMAND LINE:
## ./resultsStart.py --species='Kea' --scenario=1


