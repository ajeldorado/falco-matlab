# FALCO
Fast Linear Least-Squares Coronagraph Optimization code


%----------------------------------------------------------------------------------%
* Matlab Versions and Libraries *
%----------------------------------------------------------------------------------%

FALCO was built and tested for Matlab 2016a. It may still work on older versions, but functionality is not guaranteed.

No Matlab toolboxes should be required for FALCO. However, the Parallel Computing Toolbox or Distributed Computing Toolbox can be used to parallelize some large calculations by changing the value of a flag, "mp.flagParfor = true;". Please email the developer if you find that any other toolboxes are accidentally and/or unnecessarily used or called. There may be some instances of “imresize.m” still in the code.

FALCO was written primarily on the MacOS operating system and used to a lesser extent on Windows and Linux systems. Please report any operating system-related FALCO bugs to the developer.





%----------------------------------------------------------------------------------%
* Installation Instructions *
%----------------------------------------------------------------------------------%

1) Install the latest MATLAB version of the PROPER optical propagation library from https://sourceforge.net/projects/proper-library/ in the folder "lib/PROPER". 

2) (Optional but recommended) Install CVX from cvxr.com in the folder "lib/cvx". If you install CVX elsewhere, you will need to add the actual path of the “cvx” folder to the file “pathdef.m” that defines your Matlab path. (To find your pathdef.m file, type “which pathdef” at the Matlab command prompt.) Perform the CVX installation instructions listed on the cvxd.com website.

3) All other file paths should be known internally in the “FALCO” folder. Try running one of the template script files, which starts with "TEMPLATE_", in the folder “main” as is. If it works, then all the file paths are set correctly.

%----------------------------------------------------------------------------------%
*  *
%----------------------------------------------------------------------------------%
