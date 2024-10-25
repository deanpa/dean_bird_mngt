#!/bin/bash

#SBATCH --job-name=ScKeaMod2
#SBATCH --account=landcare00074 
#SBATCH --mail-type=end
#SBATCH --mail-user=andersond@landcareresearch.co.nz
#SBATCH --time=15:30:00

#SBATCH --mem=3000  
#SBATCH --cpus-per-task=1

export RIOS_DFLT_JOBMGRTYPE=slurm
export RIOS_SLURMJOBMGR_SBATCHOPTIONS="--job-name=sub_Sc_2 --account=landcare00074 --time=10:00:00 --mem-per-cpu=3000"
export RIOS_SLURMJOBMGR_INITCMDS="export PYTHONPATH=$PWD;module load Python-Geo/3.8.2-gimkl-2020a"

module load Python-Geo/3.8.2-gimkl-2020a

./simulationStart.py --species Kea --scenario 2