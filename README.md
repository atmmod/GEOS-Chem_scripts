# GEOS-Chem_scripts
Contains python scripts for visualizing and comparing GEOS-Chem to GEOS-Chem-hyd

## About GEOS-Chem-hyd

GEOS-chem-hyd is the hyperdual enabled version of GEOS-Chem. Below is a detailed readme for changes made in the development of GEOS-Chem-hyd.

------------------------------------------------------------------------------------------

### README CONTENTS

1.0 Contact Information

1.1 Release Contents

1.2 Requirements for Running the Model

1.3 GCClassic-hyd

1.4 geos-chem-hyd

1.5 HEMCO-hyd

1.6 GEOS-Chem_scripts

1.7 Contributing

1.8 Revision History

------------------------------------------------------------------------------------------

### 1.0 Contact Information

Name 				|			Email
--------------------|------------------------
Samuel Akinjole     |  soa52@drexel.edu
Shannon Capps		|  sc3623@drexel.edu               

------------------------------------------------------------------------------------------

### 1.1 Release Contents

The contents of the GEOS-Chem-hyd include the following repositories;
- GCClassic-hyd
- geos-chem-hyd
- HEMCO-hyd
- GEOS-Cehm_scripts
  
All changes are in the branch version_14.0.0-rc.3 of the default GEOS-Chem. GClassic-hyd is the main superwrapper with additional submodules geos-chem-hyd and HEMCO-hyd which are present in the checked out branch - version_14.0.0-rc.3, which may serve as a feature submodule when merging to GEOS-Chem. The links for GEOS-Chem and HEMCO submodules have also been set to the hyperdual repositories were developments are being made.

------------------------------------------------------------------------------------------

### 1.2 Requirements for Running the Model

The requirements for running GEOS-Chem-hyd follow the same procedure with compiling and running the default GEOS-Chem, the relevant changes to makefiles and library dependencies have been made to ensure HDMod is also compiled with the other default files.

------------------------------------------------------------------------------------------

### 1.3 GClassic-hyd

GClassic-hyd is the main superwrapper with additional submodules geos-chem-hyd and HEMCO-hyd which are present in the checked out branch - version_14.0.0-rc.3, which may serve as a feature submodule when merging to GEOS-Chem. The links for GEOS-Chem and HEMCO submodules have also been set to the hyperdual repositories were developments have been made
All changes are in the branch version_14.0.0-rc.3 

------------------------------------------------------------------------------------------

### 1.4 geos-chem-hyd

The development of geos-chem-hyd has been multifaceted, here the current active branches are described below:

- full/hyd: This is the leading hyd branch were all sub-branches are merged into.

- hetp-hyd: This branch contains all development from the full/hyd branch, in addition to the hetp implementation in GEOS-Chem from Sebastian, and the concurrent changes of hetp for hyperdual calculations. Current tests are being made on this branch using ISORROPIA for aerosol thermodynamics

- hyd/isotest: This is a dedicated branch for Isorropia testing

------------------------------------------------------------------------------------------

### 1.5 HEMCO-hyd

Developments on HEMCO have been similar to that in GEOS-Chem, there is one main hyperdual branch in HEMCO shown below:

- hemco/hyd: This contains all current developments to HEMCO files and routines with the same development version as the GEOS-Chem-hyd version. In this case, development has been on HEMCO version 3.5

------------------------------------------------------------------------------------------

### 1.6 GEOS-Chem_scripts

Here are the run scripts used in visualizing outputs from GEOS-Chem-hyd. They include the following files as listed below:
	- sensitivity.py
	- igcplot.py
	
These files are located in my main GEOS-Chem run directory; gc_4x5_merra2_fullchem. To avoid history conflicts and ease of management with subsequent changes across all run directory branches, the local repository is in a subdirectory named 'python' with symbolink links to the files in the run directory. This ensures changes made to the python scripts in the outer directory across branches are reflected in the local repository.
The branches in the remote repository are:

- main: This has all the present changes to the python script files
- scripts: This was a test run with default rundir repository - Currently obsolete.

------------------------------------------------------------------------------------------

### 1.7 Contributing

If you wish to contribute to this project, please clone or fork the repository and create a pull request with your changes. Ensure your code adheres to the existing style and include appropriate tests.

------------------------------------------------------------------------------------------

### 1.8 Revision History

Currently the version developments follow that of the default GEOS-Chem source code. The current implemented version is 14.0.0-rc.3 with some up-to-date changes to the thermodynamics file that includes hetp

------------------------------------------------------------------------------------------


<end of file>
