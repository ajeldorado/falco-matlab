README
FRNcalculator Version_2.0 (updated 5/13/19)
B. Nemati and S. Miller


INTRODUCTION
------------
This package takes a set of scenario parameters and returns the top level 
numbers in the flux ratio noise (FRN) error budget. The numbers are 
expressed in parts per billion (ppb), and the topmost number (total FRN) will be 
used for mask design optimization.


USAGE
-----
To run the core function calculateFRN.m, use one of the scripts located within a 
directory labeled with the letter s (eg. FRNscript_Best_Exo_Zodi_IMG_Band1.m 
within directory s_Best_Exo_Zodi_IMG_Band1).

Each script contains the scenario and planet inputs for a different scenario 
case. Alongside each script is its corresponding coronagraph parameters, stored 
in the directory labeled CGdata.


REQUIREMENTS
------------
All files provided must be in csv format. Until future updates, any new calls to 
calculateFRN.m must be created by mimicking the exact directory format of the 
current example scripts. When creating a new script/function call to 
calculateFRN.m, the user must provide the function inputs for calculateFRN.m, as 
well as the following directories:

	- CGdata
	- FRNscenarioData

CGdata must contain the following files, with the following filenames, generated 
by the user for the specific coronagraph in question:

	- AnnZoneList.csv
	- InitialRawContrast.csv
	- KristTable.csv
	- Sensitivities.csv
	- SensitivityMUF.csv

FRNscenarioData must contain the following files, with the following filenames. 
The user may use the files provided, although the Disturbance.csv table can be 
switched out for other disturbance cases:

	- Disturbance.csv
	- DisturbanceUnits.csv
	- HSTsynPhotSpectra.csv
	- PlanetTable.csv
	- ThroughputTable.csv
	- QE_e2v.mat


CONTENT
-------
The functions called by calculateFRN.m are as listed:

	- CStability.m
	- readDataTables.m
	- readDesignFile.m
	- readHSTsynPhotSpectra.m
	- readPlanetTable.m
	- readThroughputTable
	- unitsConstants.m
	- setDetector (within calculateFRN.m)
	- trapCTE (within calculateFRN.m)

These should remain as they are, bundled alongside calculateFRN.m inside 
FRNcode.

ForReconciliation.xlsm is a copy of WFIRST_CGI_Performance.xlsm. It is the 
spreadsheet to which all outputs of calculateFRN.m must be reconciled. In its 
current state, the error budget outputs of calculateFRN.m match the outputs of 
ForReconciliation.xlsm to 12 decimal places. This is as precise a match as 
possible when considering the floating point arithmetic of Excel versus the 
fixed point arithmetic of Matlab, and therefore should be considered an exact 
match.

