#!/bin/bash

#SBATCH --job-name=Scen1Kea
#SBATCH --account=landcare00074 
#SBATCH --mail-type=end
#SBATCH --mail-user=deanpa@pm.me
#SBATCH --time=00:05:00

#SBATCH --mem=2000  
#SBATCH --cpus-per-task=1

module load TuiView/1.2.6-gimkl-2020a-Python-3.8.2

./simulationStart.py --species Kea --scenario 1