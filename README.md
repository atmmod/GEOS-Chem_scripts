# GEOS-Chem_scripts
Contains python scripts for visualizing and comparing GEOS-Chem to GEOS-Chem-hyd

## About GEOS-Chem-hyd

GEOS-chem-hyd is the hyperdual enabled version of GEOS-Chem. Below is a detailed readme for changes made in the development of GEOS-Chem-hyd.

------------------------------------------------------------------------------------------

### README CONTENTS

1.0 Contact Information

1.1 Release Contents

1.2 Requirements for Running the Model

1.3 Git Branches

1.4 GEOS-Chem Python Scripts

1.5 HEMCO

1.6 Contributing

1.7 Revision History

------------------------------------------------------------------------------------------

### 1.0 Contact Information

Name: Samuel Akinjole                       
Email: soa52@drexel.edu

------------------------------------------------------------------------------------------

### 1.1 Release Contents

The contents of the GEOS-Chem-hyd include all files present in the duplicated GEOS-Chem repository and additional testing repositories and files as shown below
	GEOS-Chem
	HDMod.f90
	GEOS-Chem_scripts
	hydtesting
	hydfortran
	
------------------------------------------------------------------------------------------

### 1.2 Requirements for Running the Model

The requirements for running GEOS-Chem-hyd follow the same procedure with compiling and running the default GEOS-Chem, the relevant changes to makefiles and library dependencies have been made to ensure HDMod is also compiled with the other default files.

------------------------------------------------------------------------------------------

### 1.3 Git branches

The development of GEOS-Chem-hyd has been multifaceted, here the current active branches are described below:

- full/hyd: This is the leading hyd branch were all sub-branches are merged into. This includes the separate branches for handling and testing each modular operation such as transport, chemistry, cloud convection etc. These merged branches include:
	- hyd/chemistry
	- hyd/convection
	- hyd/pbl
	- hyd/transport

	The current state of the branch is just before changes made for the hetp-hyd branch, which includes the inclusion of the hetp module and the relevant changes to set hetp as a method. Currently the hetp-hyd branch is were all developments are currently being made on.

- hetp-hyd: This branch contains all development from the full/hyd branch, in addition to the hetp implementation in GEOS-Chem from Sebastian, and the concurrent changes of hetp for hyperdual calculations. Current tests are being made on this branch

- default/hetp: Includes the default implementation in GEOS-Chem by Sebastian with changes to reflect the final updated code file in the published hetp paper.

- test/hetp: This is the default implementation in GEOS-Chem by Sebastian with no changes made. This serves as a reference branch, which will soon be obselete

- test/hetp2: Similar repository to test/hetp but for my local mac repository. This is an obsolete branch.

Overall, for hetp, the implementation will have to be changed to that present in the latest GEOS-Chem classic release. However, since hetp presently is not being used for current test, this development will be done much later.

------------------------------------------------------------------------------------------

### 1.4 GEOS-Chem Python Scripts

Here are the run scripts used in visualizing outputs from GEOS-Chem-hyd. They include the following files as listed below:
	- sensitivity.py
	- igcplot.py
	
These files are located in my main GEOS-Chem run directory; gc_4x5_merra2_fullchem. To avoid history conflicts and ease of management with subsequent changes across all run directory branches, the local repository is in a subdirectory named 'python' with symbolink links to the files in the run directory. This ensures changes made to the python scripts in the outer directory across branches are reflected in the local repository.
The branches in the remote repository are:

- main: This has all the present changes to the python script files
- scripts: This was a test run with default rundir repository - Currently obsolete.

------------------------------------------------------------------------------------------

### 1.5 HEMCO

Developments on HEMCO have been similar to that in GEOS-Chem, there is one main branch in HEMCO shown below:

- hemco/hyd: This contains all current developments to HEMCO files and routines with the same development version as the GEOS-Chem-hyd version. In this case, development has been on HEMCO version 3.5



------------------------------------------------------------------------------------------

### 1.6 Contributing

If you wish to contribute to this project, please clone or fork the repository and create a pull request with your changes. Ensure your code adheres to the existing style and include appropriate tests.

------------------------------------------------------------------------------------------

### 1.7 Revision History

Currently the version developments follow that of the default GEOS-Chem source code. The current implemented version is 14.0.0-rc.3 with some up-to-date changes to the thermodynamics file that includes hetp


------------------------------------------------------------------------------------------


<end of file>