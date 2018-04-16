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

No Matlab toolboxes should be required for FALCO. If it is called anywhere, the only functions used are “imresize.m” and “padarray.m”.

FALCO was written primarily on the MacOS operating system and tested on Windows and Unix/Linux systems. Please report any FALCO bugs on Windows or Unix/Linux systems to the developers.



%----------------------------------------------------------------------------------%
* Installation Instructions *
%----------------------------------------------------------------------------------%

1) Install the PROPER optical propagation library from https://sourceforge.net/projects/proper-library/ in the folder “PROPER”. 

2) (Optional) Install CVX from cvxr.com in the folder “cvx”. If you install CVX elsewhere, you will need to add the actual path of the “cvx” folder to the file “pathdef.m” that defines your Matlab path. (To find your pathdef.m file, type “which pathdef” at the Matlab command prompt.) Perform the CVX installation instructions listed on the cvxd.com website.

3) All other file paths should be known internally in the “FALCO” folder. Try running one of the files in the folder “run_templates” as is. If it works, then all the file paths are set correctly.

%----------------------------------------------------------------------------------%
*  *
%----------------------------------------------------------------------------------%
