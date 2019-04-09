% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Vout = falco_fit_dm_surf(dm,Vin)
% 
%--Function to fit a surface to a deformable mirror (DM) commands using PROPER.
%
%--INPUTS
% dm: structure of DM parameters
% Vin: 2-D array of desired DM surface heights at each actuator
%
%--OUTPUT
% Vout: 2-D array of output DM voltage commands to give the desired surface
%       heights in Vin
% 
%--VERSION HISTORY
% Created on 2019-01-23 by A.J. Riggs.

function Vout = falco_fit_dm_surf(dm,Vin)

%--Starting influence function
inf1 = dm.inf0;
N1 = length(inf1);
actres1 = dm.dm_spacing/dm.dx_inf0;
x = (-(N1-1)/2:(N1-1)/2)/actres1;
[X,Y] = meshgrid(x);

%--Influence function resampled to actuator map resolution
actres2 = 1; %--pixels per actuator width
N2 = ceil_even(N1*actres2/actres1)+1; %--Make odd to have peak of 1
xq = (-(N2-1)/2:(N2-1)/2)/actres2;
[Xq,Yq] = meshgrid(xq);
inf2 = interp2(X,Y,inf1,Xq,Yq,'cubic',0);

%--Perform the fit
[Vout, ~] = prop_fit_dm(Vin, inf2);

end %--END OF FUNCTION