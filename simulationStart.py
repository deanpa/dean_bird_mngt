#!/usr/bin/env python

import sys
import optparse
import simulationMain
import simulationParams


class CmdArgs(object):
    def __init__(self):
        p = optparse.OptionParser()
        p.add_option("--species", dest="species", help="Set species")
        p.add_option("--scenario", dest="scenario", help="Set scenario")
        (options, args) = p.parse_args()
        self.__dict__.update(options.__dict__)

cmdargs = CmdArgs()

params = simulationParams.PreyParams(cmdargs.species, int(cmdargs.scenario))

simulationMain.main(params)

