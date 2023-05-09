#!/usr/bin/env python

import sys
import optparse
from modelScripts import simulationMain
import simulationParamsDean


class CmdArgs(object):
    def __init__(self):
        p = optparse.OptionParser()
        p.add_option("--species", dest="species", help="Set species")
        p.add_option("--scenario", dest="scenario", help="Set scenario")
        (options, args) = p.parse_args()
        self.__dict__.update(options.__dict__)

if __name__ == '__main__':


    cmdargs = CmdArgs()

    params = simulationParamsDean.PreyParams(cmdargs.species, int(cmdargs.scenario))
    #params = simulationParams.PreyParams("Kea", int(8))

    simulationMain.main(params)



## ON COMMAND LINE, RUN SCRIPT LIKE THIS FOR KEA AND SCENARIO 1:
### ./simulationStart.py --species Kea --scenario 1
