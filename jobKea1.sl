#!/bin/bash

#SBATCH --job-name=Scen8
#SBATCH --account=landcare00074 
#SBATCH --mail-type=end
#SBATCH --mail-user=andersond@landcareresearch.co.nz
#SBATCH --time=00:40:00

#SBATCH --mem=2500  
#SBATCH --cpus-per-task=1

export RIOS_DFLT_JOBMGRTYPE=slurm
export RIOS_SLURMJOBMGR_SBATCHOPTIONS="--job-name=Sc8_sub --account=landcare00074 --time=00:30:00 --mem-per-cpu=2500"
export RIOS_SLURMJOBMGR_INITCMDS="export PYTHONPATH=$PWD;module load Python-Geo/3.8.2-gimkl-2020a"

module load Python-Geo/3.8.2-gimkl-2020a

./simulationStart.py --species Kea --scenario 8