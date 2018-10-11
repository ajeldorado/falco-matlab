# FALCO
Fast Linear Least-Squares Coronagraph Optimization (FALCO) software package

Refer to the SPIE conference paper "Fast Linearized Coronagraph Optimizer (FALCO) I: A software toolbox for rapid coronagraphic design and wavefront correction" for an overview of FALCO and its uses. 
DOI: 10.1117/12.2313812

Developed by A.J. Riggs at the Jet Propulsion Laboratory, California Institute of Technology.
Significant contributions and testing were done by Garreth Ruane, Erkin Sidick, and Carl Coker.

Copyright 2018, by the California Institute of Technology. ALL RIGHTS RESERVED. 
United States Government Sponsorship acknowledged. 
Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
%----------------------------------------------------------------------------------%


 

%----------------------------------------------------------------------------------%
* Overview *
%----------------------------------------------------------------------------------%

The Fast Linearized Coronagraph Optimizer (FALCO) is an open-source toolbox of routines for coronagraphic focal plane wavefront correction. The goal of FALCO is to provide a free, modular framework for the simulation or testbed operation of several common types of coronagraphs, and the design of coronagraphs that use wavefront control algorithms to shape deformable mirrors (DMs) and masks. FALCO includes routines for pair-wise probing estimation of the complex electric field and Electric Field Conjugation (EFC) control, and we ask the community to contribute other wavefront correction algorithms. FALCO utilizes and builds upon PROPER, an established optical propagation library. The key innovation in FALCO is the rapid computation of the linearized response matrix for each DM, which facilitates re-linearization after each control step for faster DM-integrated coronagraph design and wavefront correction experiments. FALCO is freely available as source code in MATLAB at github.com/ajeldorado/falco-matlab and will be available later this year in Python 3 at github.com/ajeldorado/falco-python.


%----------------------------------------------------------------------------------%
* Matlab Versions and Libraries *
%----------------------------------------------------------------------------------%

FALCO was built and tested in Matlab 2017a and b. It may still work on older versions, but functionality is not guaranteed.

No Matlab toolboxes should be required for FALCO. However, the Parallel Computing Toolbox or Distributed Computing Toolbox can be used to parallelize some large calculations by changing the value of a flag, "mp.flagParfor = true;". Please email the developer if you find that any other toolboxes are accidentally and/or unnecessarily used or called. There may be some instances of “imresize.m” still in the code.

FALCO was written primarily on the MacOS operating system and used to a lesser extent on Windows and Linux systems. Please report any operating system-related FALCO bugs to the developer.





%----------------------------------------------------------------------------------%
* Installation Instructions *
%----------------------------------------------------------------------------------%

1) Linking to PROPER (in MATLAB): Download the latest MATLAB version of the PROPER optical propagation library from https://sourceforge.net/projects/proper-library/. Wherever you decide to unzip and place the PROPER library on your machine, add the path to the PROPER folder to the MATLAB path. 
  A) You can temporarily do that by defining the variable "mp.path.proper" in each main script you use, or 
  B) (better, but may not work on servers with restricted write permissions) Permanently adding PROPER to the MATLAB path with the commands 
       addpath(path/to/proper); savepath;
     where "path/to/proper" is the absolute file path on your computer for the PROPER directory. The command "savepath" will keep the directory you included in the "pathdef.m" file that MATLAB uses to look for the functions is expects. The pathdef.m file might not be writeable on a server without admin privileges. 

2) (Optional--not used for regular functionality.) Download CVX from cvxr.com. Wherever you decide to unzip and place the PROPER library on your machine, add the path to the CVX directory to the MATLAB path (similar to how it was done for PROPER). Then, perform the CVX installation instructions listed on the cvxd.com website.

3) Try to run one of the template script files, which starts with "TEMPLATE_", in the folder “main” as it is, except for adjusting the file path definitions listed near the very top. FALCO must know where the FALCO library resides (given by "mp.path.falco") and where the PROPER library resides (given by "mp.path.proper" if it is not in MATLAB's path already). If the template script runs through with no errors, then all the file paths are set correctly.

%----------------------------------------------------------------------------------%
*  *
%----------------------------------------------------------------------------------%
