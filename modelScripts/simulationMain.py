#!/usr/bin/env python

"""
Set the following variables for a SLURM run on Pan:

export RIOS_DFLT_JOBMGRTYPE=slurm
export RIOS_SLURMJOBMGR_SBATCHOPTIONS="--job-name=preyTest --account=landcare00074 --time=00:06:00 --mem-per-cpu=10000"
export RIOS_SLURMJOBMGR_INITCMDS="export PYTHONPATH=$PWD;module load Python-Geo"

"""

import os
import multiprocessing
import pickle
import shutil
#import calculation
#import preProcessing
from modelScripts import calculation
from modelScripts import preProcessing
from modelScripts import calcresults
from rios.parallel import jobmanager

# Use the same environment variable as RIOS to define the type of
# parallel processing.
# Default to the multiprocessing type.
JOBMGR_TYPE = os.getenv('RIOS_DFLT_JOBMGRTYPE', default='multiprocessing')


# from mahuika transition page
#SCRATCH_DIR=$(mktemp -d --tmpdir=/nesi/nobackup/$SLURM_JOB_ACCOUNT "scratch_${SLURM_JOB_ID}_XXX.tmp")

## TEMP SCRATCH DIRECTORY 
#NESI_TMP_DIR=$(mktemp -d --tmpdir=/nesi/nobackup/landcare00074)

## TEMP SCRATCH DIRECTORY 
NESI_TMP_DIR=os.path.join(os.sep, 'nesi', 'nobackup', 'landcare00074')



def parallelRunModel(data, iteration, results):
    """
    A slight variation on calculation.runModel
    which makes the results a parameter so it 
    can be used with rios.parallel.
    """
    newresults = calculation.runModel(data, loopIter=iteration)
    results.params = newresults.params
    results.rodentDensity_2D = newresults.rodentDensity_2D
    results.stoatDensity_2D = newresults.stoatDensity_2D
    results.preyDensity_2D = newresults.preyDensity_2D
    results.popAllYears_3D = newresults.popAllYears_3D
    results.rodentDensity_2D_mth = newresults.rodentDensity_2D_mth
    results.stoatDensity_2D_mth = newresults.stoatDensity_2D_mth
    results.preyDensity_2D_mth = newresults.preyDensity_2D_mth

    results.controlCount = newresults.controlCount


class PreyJobInfo(jobmanager.JobInfo):
    """
    Contains an implementation of RIOS's jobmanager.JobInfo
    for the prey model.
    """
    def __init__(self, data, iteration):
        self.data = data
        self.iteration = iteration

    def getFunctionParams(self):
        "make input suitable for parallelRunModel"
        results = calcresults.PreyResults()
        return self.data, self.iteration, results

    def getFunctionResult(self, params):
        "output was the last parameter"
        return params[-1]
    
def runMultipleJobs(data, resultsFName):
    # if using multiprocessing, run a job per cpu
    # otherwise (assume SLURM) run a job per iteration
    # not sure if this is correct
    if JOBMGR_TYPE == 'multiprocessing':
        nThreads = multiprocessing.cpu_count()
    else:
        nThreads = data.params.iter
    jobmgrClass = jobmanager.getJobManagerClassByType(JOBMGR_TYPE)
    jobmgr = jobmgrClass(nThreads)

    # Home dir runs out of quota
    # I couldn't find a cluster-wide temp var created on Pan
    # so simply use this dir if it exists (assume we are on Pan)
    # otherwise leave as default.
    if os.path.isdir(NESI_TMP_DIR):
        jobmgr.setTempdir(NESI_TMP_DIR)

    jobInputs = []
    for i in range(data.params.iter):

#        print('iter', i)

        jobInfo = PreyJobInfo(data, i)
        jobInputs.append(jobInfo)
    # run all in parallel and collect results
    results = jobmgr.runSubJobs(parallelRunModel, jobInputs)

    # pickle results
    fileobj = open(resultsFName, 'wb')
    pickle.dump(results, fileobj, protocol=4) # so we get large file support
    fileobj.close()



########  MAIN FUNCTION  ########
#######
def main(params):
    ## PREPROCESSING AND RESULTS PICKLE NAMES

    preProcFName = params.preProcFName

    resultsFName = params.resultsFName

    ## MAKE PARAMS FILE NAME TO COPY TO RESULTS DIRECTORY FOR REFERENCE
    paramsFName = 'params_' + params.species + 'Scen' + str(params.scenario) + '.py'
    paramsFName = os.path.join(params.outputDataPath, paramsFName)
    shutil.copy('simulationParams.py', paramsFName)

    ## IF FIRST RUN, FROM PARAMS SCRIPT, RUN PREPROCESSING, ELSE UNPICKLE
    if params.firstRun:
        data = preProcessing.PreyData(params)
        data.pickleSelf(preProcFName)
    else:
        data = preProcessing.PreyData.unpickleFromFile(preProcFName)

    runMultipleJobs(data, resultsFName)


if __name__ == '__main__':
    main()    
