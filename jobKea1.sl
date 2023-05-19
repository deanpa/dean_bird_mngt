#!/bin/bash

#SBATCH --job-name=Scen1Kea
#SBATCH --account=landcare00074 
#SBATCH --mail-type=end
#SBATCH --mail-user=deanpa@pm.me
#SBATCH --time=00:45:00

#SBATCH --mem=6000  
#SBATCH --cpus-per-task=1

export RIOS_DFLT_JOBMGRTYPE=slurm
export RIOS_SLURMJOBMGR_SBATCHOPTIONS="--job-name=subMod1 --account=landcare00074 --time=00:35:00 --mem-per-cpu=6000"
export RIOS_SLURMJOBMGR_INITCMDS="export PYTHONPATH=$PWD;module load Python-Geo/3.8.2-gimkl-2020a"

module load Python-Geo/3.8.2-gimkl-2020a

./simulationStart.py --species Kea --scenario 1