# FALCO

Fast Linear Least-Squares Coronagraph Optimization (FALCO) software package

Refer to the SPIE conference paper "Fast Linearized Coronagraph Optimizer (FALCO) I: A software toolbox for rapid coronagraphic design and wavefront correction" for an overview of FALCO and its uses. 
DOI: 10.1117/12.2313812

Developed by A.J. Riggs at the Jet Propulsion Laboratory, California Institute of Technology.
Major contributions and testing were done by Garreth Ruane, Erkin Sidick, Carl Coker, Navtej Saini, and Jorge Llop-Sayson.

 
# What is FALCO?

The Fast Linearized Coronagraph Optimizer (FALCO) is an open-source toolbox of routines and example scripts for coronagraphic focal plane wavefront correction. The goal of FALCO is to provide a free, modular framework for the simulation or testbed operation of several common types of coronagraphs, and the design of coronagraphs that use wavefront control algorithms to shape deformable mirrors (DMs) and masks. FALCO includes routines for pair-wise probing estimation of the complex electric field and Electric Field Conjugation (EFC) control, and we ask the community to contribute other wavefront correction algorithms and optical layouts. FALCO utilizes and builds upon PROPER, an established optical propagation library. The key innovation in FALCO is the rapid computation of the linearized response matrix for each DM, which facilitates re-linearization after each control step for faster DM-integrated coronagraph design and wavefront correction experiments. FALCO is freely available as source code in MATLAB at github.com/ajeldorado/falco-matlab and is in development in Python 3.

# Documentation and Support

FALCO is provided as-is and has no guarantee of performance. Nevertheless, reasonable attempts have been made to debug and troubleshoot the code, and the developers are still using and improving the software.

#### DOCUMENTATION:  
For now, the documentation is available at the Github Wiki at https://github.com/ajeldorado/falco-matlab/wiki.


# Matlab Versions and Libraries

FALCO was built and tested in Matlab 2017a and b. It may still work on older versions, but functionality is not guaranteed.

No Matlab toolboxes should be required for FALCO. However, the Parallel Computing Toolbox or Distributed Computing Toolbox can be used to parallelize some repetitive calculations by changing the value of a flag, *mp.flagParfor = true;*. 

Please email the developer if you find that any other toolboxes are accidentally and/or unnecessarily used or called. There may be some instances of *imresize.m* still in the code. Calls to *rms.m* and *sinc.m* have been removed since those simple functions exist in the Signal Processing Toolbox. Thank you to Jason Kay for reporting the rms issue.

FALCO was written primarily on the MacOS operating system and used to a lesser extent on Windows and Linux systems. Please report any operating system-related FALCO bugs to the developer.


# Installation Instructions

1) You need a MATLAB license and an install of MATLAB. Multi-wavelength simulations may need a desktop computer or server instead of a laptop to run. Monochromatic and/or lower-resolution trials usually run quickly on a laptop with a few GB of RAM.

2) Linking to PROPER (in MATLAB): Download the latest MATLAB version of the PROPER optical propagation library from SourceForge at [https://sourceforge.net/projects/proper-library/](https://sourceforge.net/projects/proper-library/). Wherever you decide to unzip and place the PROPER library on your machine, add the path to the PROPER folder to the MATLAB path. 
  A) You can temporarily do that by defining the variable `mp.path.proper` in each main script of config file you use, or 
  B) (better, but may not work on servers with restricted write permissions) Permanently adding PROPER to the MATLAB path with the commands 
       `addpath(path/to/proper); savepath;`
     where _path/to/proper_ is the absolute file path on your computer for the PROPER directory. The command `savepath` will keep the directory you included in the *pathdef.m* file that MATLAB uses to look for the functions is expects. The *pathdef.m* file might not be writeable on a server without admin privileges. 

3) Try to run one the example script file _EXAMPLE_try_running_FALCO.m_ in the sub-folder *falco-matlab/main*. The only needed modifications are the FALCO and PROPER file path definitions listed near the very top of the script to be for your computer. FALCO must know where the FALCO library resides (given by `mp.path.falco`) and where the PROPER library resides (given by `mp.path.proper` if it is not in MATLAB's path already). If the example script runs through with no errors, then all the file paths are set correctly.

4) Now go ahead and try some of the other example scripts in _falco-matlab/main_, which start with "EXAMPLE_", again adjusting the path definitions for the FALCO and PROPER libraries. For this initial functionality test, in the main script you should set the number of wavelengths to 1 (`mp.Nsbp = 1;`) and turn off the parallel computing  flag (`mp.flagParfor = false;`) for it to run quickly. I recommend finding the script closest to your intended purpose and starting to make changes from a copy of that.

5) (Optional--not used for regular functionality.) Download CVX from [cvxr.com](cvxr.com). Wherever you decide to unzip and place the PROPER library on your machine, add the path to the CVX directory to the MATLAB path (similar to how it was done for PROPER). Then, perform the CVX installation instructions listed on the CVX website.


# Version History

Upcoming: In Version 3.1, there will be the first "full models" that differ from the compact model and look like actual testbeds or instrument layouts. More coronagraph types will be added back in that were not fully converted in the Version 3.0 release.

February 28, 2018: Version 3.0 released. Version 3 should be the first stable version of FALCO and be supported for some time. Version 3 was re-designed to be expandable for many people's different needs and uses without them having to overwrite each other's features.
  - Wavefront estimation added (pairwise probing with batch process and Kalman filter).
  - propcustom_dm.m added to allow use of different DM actuator influence functions.
  - Config file loaded first in the main script rather than afterward. Makes everything much easier to read and modify.
  - Optical models top-level changed from coronagraph type to layout (i.e., which bench or instrument). This makes FALCO more expandable for different people's uses.
  - Scalar vortex model added from Garreth Ruane.
  - Many changes to syntax and variable names to make code more uniform, simpler, and easier to expand upon.

October 10, 2018:  Version 2.0 released. HLC design code added. Many more features added and syntax changes.

April 11, 2018:    Version 1.0 released. Wavefront control functionality for LC, SPLC, and VC coronagraphs.


# Legal Notices

Copyright 2018-2019. California Institute of Technology ("Caltech"). This software, including source and object code, and any accompanying documentation ("Software") is owned by Caltech. Caltech has designated this Software as Technology and Software Publicly Available ("TSPA"), which means that this Software is publicly available under U.S. Export Laws. With the TSPA designation, a user may use and distribute the Software on a royalty-free basis with the understanding that:

(1) THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE; and

(2) Caltech is under no obligation to provide technical support for the Software; and

(3) all copies of the Software released by user must be marked with this marking language, inclusive of the copyright statement, TSPA designation and user understandings.

