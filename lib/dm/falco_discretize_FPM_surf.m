% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to discretize the surface profiles of the metal and dielectric
% layers of the FPM.
%
% OUTPUT:
% - DMtransInd: array index of the metal/dielectric layer for each DM8/9 surface value in the complex transmission matrix 
%
% NOTE:
% - DM8 is the metal layer
% - DM9 is the dielectric layer
%
% REVISION HISTORY:
% - Created on 2018-05-07 by A.J. Riggs

function DMtransInd = falco_discretize_FPM_surf(DMsurf,t_nm_vec, dt_nm)

%--Convert surface profiles from meters to nanometers
DMsurf = 1e9*DMsurf;

%--Stay within the thickness range since material properties are not defined outside it
DMsurf(DMsurf<min(t_nm_vec)) = min(t_nm_vec);
DMsurf(DMsurf>max(t_nm_vec)) = max(t_nm_vec);

%--Discretize to find the index the the complex transmission array
DMtransInd = 1 + round(1/dt_nm*(DMsurf - min(t_nm_vec)));

end %--END OF FUNCTION