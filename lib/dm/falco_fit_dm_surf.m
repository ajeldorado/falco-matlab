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
% surfaceToFit: 2-D array of desired DM surface heights at each actuator
%
%--OUTPUT
% Vout: 2-D array of output DM voltage commands to give the desired surface
%       heights in Vin
% 
%--VERSION HISTORY
% Created on 2019-01-23 by A.J. Riggs.

function Vout = falco_fit_dm_surf(dm,surfaceToFit)

%--Set the order of operations
orderOfOps = 'XYZ';
if(isfield(dm,'flagZYX'))
    if(dm.flagZYX)
        orderOfOps = 'ZYX'; 
    end
end

%--Error checks
[mSurface, nSurface] = size(surfaceToFit);
if(mSurface ~= nSurface); error('surfaceToFit must be a square matrix.'); end
if(mod(nSurface,2) ~= 0); error('surfaceToFit must have an even number of points.'); end

%--Original influence function
infFuncAtOrigRes = dm.inf0;
nPixAcrossOrigInfFunc = length(infFuncAtOrigRes);
origPixPerAct = dm.dm_spacing/dm.dx_inf0;
xOrig = (-(nPixAcrossOrigInfFunc-1)/2:(nPixAcrossOrigInfFunc-1)/2)/origPixPerAct;
[Xorig,Yorig] = meshgrid(xOrig);

%--Influence function resampled to actuator map resolution
newPixPerAct = 1; %--pixels per actuator width
nPixAcrossNewInfFunc = ceil_even(nPixAcrossOrigInfFunc*newPixPerAct/origPixPerAct)+1; %--Make odd to have peak of 1
xNew = (-(nPixAcrossNewInfFunc-1)/2:(nPixAcrossNewInfFunc-1)/2)/newPixPerAct;
[Xnew,Ynew] = meshgrid(xNew);
infFuncAtActRes = interp2(Xorig,Yorig,infFuncAtOrigRes,Xnew,Ynew,'cubic',0);


%--Perform the fit
if(nSurface == dm.Nact)
    gridDerotAtActRes = surfaceToFit;
    
elseif(nSurface > dm.Nact)
    %--Adjust the centering of the output DM surface. The shift needs to be in
    %units of actuators, not meters, for prop_dm.m.
    wArray = nSurface*dm.dx;
    switch dm.centering % 0 shift for pixel-centered pupil, or -Darray/2/Narray shift for inter-pixel centering
        case {'interpixel'}
            cshift = -wArray/2/nSurface/dm.dm_spacing; 
        case {'pixel'}
            cshift = 0;
        otherwise
            error('falco_gen_dm_surf: centering variable must be either pixel or interpixel')
    end

    gridDerotAtActRes = propcustom_derotate_resize_dm_surface(surfaceToFit, ...
        dm.dx, dm.Nact, dm.xc-cshift, dm.yc-cshift, dm.dm_spacing,...
        'XTILT',dm.xtilt,'YTILT',dm.ytilt,'ZTILT',dm.zrot,orderOfOps,...
        'inf_sign',dm.inf_sign, 'inf_fn', dm.inf_fn);
    
elseif(nSurface < dm.Nact)
    error('surfaceToFit cannot be smaller than [Nact x Nact].')
    
end

[Vout, ~] = prop_fit_dm(gridDerotAtActRes, infFuncAtActRes);


end %--END OF FUNCTION