% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function DMsurf = falco_gen_dm_surf(dm,dx,N)
% 
%--Function to generate a deformable mirror (DM) surface using PROPER.
%
%--INPUTS
% dm: structure of DM parameters
%
%--OUTPUT
% DMsurf: DM surface in meters
% 
%--VERSION HISTORY
% - Modified on 2019-02-14 by G. Ruane and D. Echeverri to allow for custom
% infuence functions and quantizing the DM commands. 
% - Modified on 2018-02-06 by A.J. Riggs to have dx and N be inputs as well 
% in order to avoid having to define a few variables (such as the DM 
% commands) twice for the compact and full models.
% - Created falco_gen_dm_poke_cube_PROPER.m on 2017-11-17 by A.J. Riggs.
%
% Need to fix:
% - Don't allow both dm.dx and dx to be used.
%
function DMsurf = falco_gen_dm_surf(dm,dx,N)

% fprintf('falco_gen_dm_surf: dx = %.8g   dm.dx = %.8g\n',dx,dm.dx); %--DEBUGGING

%--Set the order of operations
orderOfOps = 'XYZ';
if(isfield(dm,'flagZYX'))
    if(dm.flagZYX)
        orderOfOps = 'ZYX'; 
    end
end

%--Adjust the centering of the output DM surface. The shift needs to be in
%units of actuators, not meters, for prop_dm.m.
Darray = dm.NdmPad*dm.dx;
Narray = dm.NdmPad;
switch dm.centering % 0 shift for pixel-centered pupil, or -Darray/2/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -Darray/2/Narray/dm.dm_spacing; 
    case {'pixel'}
        cshift = 0;
    otherwise
        error('falco_gen_dm_surf: centering variable must be either pixel or interpixel')
end

pupil_ratio = 1; % beam diameter fraction
wl_dummy = 1e-6; %--dummy value needed to initialize wavelength in PROPER (meters)

bm  = prop_begin(N*dx, wl_dummy, N, pupil_ratio);
%[~,DMsurf] = prop_dm(bm, dm.VtoH.*dm.V, dm.xc-cshift, dm.yc-cshift, dm.dm_spacing,'XTILT',dm.xtilt,'YTILT',dm.ytilt,'ZTILT',dm.zrot,orderOfOps);

%--Quantization of DM actuation steps based on least significant bit of the
% DAC (digital-analog converter). In height, so called HminStep
% If HminStep (minimum step in H) is defined, then quantize the DM voltages
if(isfield(dm,'HminStep') && ~any(isnan(dm.HminStep(:))))

    % If desired method is not defined, set it to the default. 
    if(~isfield(dm,'HminStepMethod'))
        dm.HminStepMethod = 'round';
    end

    % Discretize/Quantize the DM voltages (creates dm.Vquantized)
	dm = falco_discretize_dm_surf(dm, dm.HminStepMethod);
    H = dm.VtoH.*dm.Vquantized;
else

    % Quantization not desired; send raw, continuous voltages
    H = dm.VtoH.*dm.V;
end

[~,DMsurf] = propcustom_dm(bm, H, dm.xc-cshift, dm.yc-cshift, dm.dm_spacing,'XTILT',dm.xtilt,'YTILT',dm.ytilt,'ZTILT',dm.zrot,orderOfOps,...
    'inf_sign',dm.inf_sign, 'inf_fn', dm.inf_fn);


end %--END OF FUNCTION








