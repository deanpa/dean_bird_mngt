#!/usr/bin/env python

import os
import sys
import importlib
import optparse


class CmdArgs(object):
    def __init__(self):
        p = optparse.OptionParser()
        p.add_option("--species", dest="species", help="Set species")
        p.add_option("--scenario", dest="scenario", help="Set scenario")
        (options, args) = p.parse_args()
        self.__dict__.update(options.__dict__)

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
## APPEND SYS DIRECTORY TO IMPORT PARAMS MODULE
basepath = os.getenv('KIWIPROJDIR', default = '.')
if basepath == '.':
    pathDir = os.path.join('/home/dean/workfolder/projects/dean_bird_mngt', resScenPath)
else:
    pathDir = os.path.join(basepath, resScenPath)
sys.path.append(pathDir)

print('path to results', pathDir)


