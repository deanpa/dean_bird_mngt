#!/bin/bash

#SBATCH --job-name=Scen1Kea
#SBATCH --account=landcare00074 
#SBATCH --mail-type=end
#SBATCH --mail-user=andersond@landcareresearch.co.nz
#SBATCH --time=1:00:00

#SBATCH --mem=2000  
#SBATCH --cpus-per-task=1

export RIOS_DFLT_JOBMGRTYPE=slurm
export RIOS_SLURMJOBMGR_SBATCHOPTIONS="--job-name=subMod1 --account=landcare00074 --time=00:40:00 --mem-per-cpu=2000"
export RIOS_SLURMJOBMGR_INITCMDS="export PYTHONPATH=$PWD;module load Python-Geo/3.8.2-gimkl-2020a"

module load Python-Geo/3.8.2-gimkl-2020a

./simulationStart.py --species Kea --scenario 2