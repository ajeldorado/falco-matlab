% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

# FALCO
Fast Linear Least-Squares Coronagraph Optimization code

Please cite FALCO if you use it as part of the work for a publication. For now, the paper to cite is the SPIE Conference Proceeding for the 2018 SPIE Astronomical Telescopes + Instrumentation conference. The citation will not be available until July 2018. The talk will be:
“Fast linearized coronagraph optimizer (FALCO) I: a software toolbox for rapid coronagraphic design and wavefront correction”
Paper 10698-101 by A.J. Riggs


%----------------------------------------------------------------------------------%
* Matlab Versions and Libraries *
%----------------------------------------------------------------------------------%

FALCO was built and tested for Matlab 2016a. It may still work on older versions, but functionality is not guaranteed.

For ease of use, no Matlab toolboxes should be required for FALCO. Please notify the developer if this is not the case.

FALCO was written primarily on the MacOS operating system and tested on Windows and Unix/Linux systems. Please report any FALCO bugs on Windows or Unix/Linux systems to the developers.



%----------------------------------------------------------------------------------%
* Installation Instructions *
%----------------------------------------------------------------------------------%

1) Install the PROPER optical propagation library from https://sourceforge.net/projects/proper-library/ in the folder “lib/PROPER”. 

2) All other file paths should be known internally in the “FALCO” folder. Try running as-is one of the template files (starting with "TEMPLATE_") in the folder “main”. If it works, then all the file paths are set correctly.

3) (Optional, and not recommended until the constrained EFC controller is tested and characterized further.) Install CVX from cvxr.com in the folder “cvx”. If you install CVX elsewhere, you will need to add the actual path of the “cvx” folder to the file “pathdef.m” that defines your Matlab path. (To find your pathdef.m file, type “which pathdef” at the Matlab command prompt.) Perform the CVX installation instructions listed on the cvxd.com website. Note that the professional version of CVX with Gurobi is required; all the free solvers are too slow.

%----------------------------------------------------------------------------------%
*  *
%----------------------------------------------------------------------------------%
