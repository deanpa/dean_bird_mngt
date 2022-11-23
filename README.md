# dean_bird_mngt
Broadscale simulation model on monthly time step to forecast predator control and bird outcomes

## 5 FILES ASSOCIATED WITH RUNNING THE SIMULATION
1) ./simulationParams.py
2) ./simulationStart.py
3) ./modelScripts/simulationMain.py
4) ./modelScripts/preProcessing.py
5) ./modelScripts/calculation.py

## BRIEF DETAILS ON THE 5 FILES FOR RUNNING THE SIMULATION
1) ./simulationParams.py
This is the only file that the user should modify. The user identifies the species, scenario number  and all parameters. The setting of the species and scenario number will create the Results directory, if it doesn't already exist.

2) ./simulationStart.py
The users starts the simulation on the command line by running this script. An example for kea and scenario 1 is as follows:

./simulationStart.py --species Kea --scenario 1

This will create a directory where the results will be saved, for example:  ./SpeciesProjects/Kea/Results/Scen1Kea

3) ./modelScripts/simulationMain.py
This script is located in the modelScripts directory. It sets up the parallel processing and calls the processing scripts.

4) ./modelScripts/preProcessing.py
This script is called if it is the first time the scenario is run (see self.firstRun in simulationParams.py). It does the initial processing of GIS data. If it is not the first run, it saves time to set the self.firstRun = True.

5) ./modelScripts/calculation.py
This script runs loops through the years of the simulation and has the processing functions.

## 3 FILES ASSOCIATED WITH PROCESSING RESULTS

1) ./resultsStart.py
The users start processing the results on the command line by running this script. It calls the params module specified by the species and scenario. An example for kea and scenario 1 is as follows:

./resultsStart.py --species Kea --scenario 1

2) ./modelScripts/resultsMain.py
Reads in the saved pickles from preProcessing and calculation.py, and calls the functions to make graphs, tables and images.

3) ./modelScripts/resultsFX.py
This script holds the functions for producing graphs, tables and images, which are written to the results directory

