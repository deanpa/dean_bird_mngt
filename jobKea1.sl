#!/bin/bash

#SBATCH --job-name=Scen3
#SBATCH --account=landcare00074 
#SBATCH --mail-type=end
#SBATCH --mail-user=andersond@landcareresearch.co.nz
#SBATCH --time=09:30:00

#SBATCH --mem=2500  
#SBATCH --cpus-per-task=1

export RIOS_DFLT_JOBMGRTYPE=slurm
export RIOS_SLURMJOBMGR_SBATCHOPTIONS="--job-name=Sc4_Sub --account=landcare00074 --time=04:45:00 --mem-per-cpu=2500"
export RIOS_SLURMJOBMGR_INITCMDS="export PYTHONPATH=$PWD;module load Python-Geo/3.8.2-gimkl-2020a"

module load Python-Geo/3.8.2-gimkl-2020a

./simulationStart.py --species Kea --scenario 3