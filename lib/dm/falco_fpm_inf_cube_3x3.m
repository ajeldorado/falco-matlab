% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute the influence functions for a DM and save them as a
%  datacube. Uses interpolation to resize, rotate, and translate the
%  influence functions to the actuator centers.
%
% --Some Key Outputs: 
% dm.inf_datacube = datacube of influence functions
% dm.NdmPad = number of points to pad the beam to to include all the actuators.
% dm.r_cent_act = distances of each actuator from the beam center (for zeroing out ones
%   outside the pupil later.
% dm.xy_cent_act_box = positions of the influence functions in the padded
%   beam array. This is essential to adding the influence functions
%   together later via superposition. Units of pixels assuming dm.Ndmpad
%   points across.
%
%--VERSION CHANGE HISTORY
% -Modified by A.J. Riggs on September 23, 2017 to allow for pixel centering.
%  Also added model_params (aka, mp or mp) as an input. Now computes NboxAS and 
%  dmPad based on the angular spectrum padding necessary for inter-DM propagation. This 
%  is to prevent the indexing outside the array if the Fresnel number becomes too small 
%  for a hard-coded value of the padding of the full DM array.
% -Written by A.J. Riggs on September 1, 2016.
%
%--SUGGESTED CHANGES
% -Clean up uses of dm.dx_dm vs dm.dx_dm to avoid any possible differences
% since they need to be the same value.

%-------------------------------------------------------------------
function dm = falco_fpm_inf_cube_3x3(dm)

fprintf('Computing datacube of FPM influence functions... ');

%--Compute sampling of the pupil. Assume that it is square.
dm.dx_dm = dm.dxi;

%--Default to being centered on a pixel (FFT-style) if no centering is specified
if(~isfield(dm,'centering'))
    dm.centering = 'pixel'; %--Centered on a pixel (default if not specified)
end

%--Compute coordinates of original influence function
Ninf0 = length(dm.inf0);
x_inf0 = (-(Ninf0-1)/2:(Ninf0-1)/2)*dm.dx_inf0; % True for even- or odd-sized influence function maps as long as they are centered on the array.
[Xinf0,Yinf0] = meshgrid(x_inf0);

%--Compute list of initial actuator center coordinates (in actuator widths).
%--Square grid
[dm.Xact,dm.Yact] = meshgrid((0:dm.Nact-1)-dm.xcent_dm,(0:dm.Nact-1)-dm.ycent_dm); % in actuator widths
x_vec = dm.Xact(:);
y_vec = dm.Yact(:);

dm.NactTotal = length(x_vec); %--Total number of actuators in the 2-D array

dm.infMaster = dm.inf0;
Nbox = length(dm.inf0);
dm.Nbox = Nbox;
fprintf('FPM influence function size =\t%dx%d ',Nbox,Nbox);

dm.xy_cent_act  = [x_vec,y_vec].'; % in actuator widths

%% Pad the pupil to at least the size of the DM(s) surface(s) to allow all actuators to be located outside the pupil.
% (Same for both DMs)

%--Find actuator farthest from center:
dm.r_cent_act = sqrt(dm.xy_cent_act(1,:).^2 + dm.xy_cent_act(2,:).^2);
dm.rmax = max(abs(dm.r_cent_act));
dm.absxymax = max( abs(dm.xy_cent_act(:)) );
NpixPerActWidth = dm.dm_spacing/dm.dx_dm;

dm.NdmPad = 0 + ceil_even( (dm.Nact+2)*NpixPerActWidth ); % prevent indexing outside the array 

%--Compute coordinates (in meters) of the full DM array
if(strcmpi(dm.centering,'pixel')  ) 
    dm.x_pupPad = (-dm.NdmPad/2:(dm.NdmPad/2 - 1))*dm.dx_dm; % meters, coords for the full DM arrays. Origin is centered on a pixel
else
    dm.x_pupPad = (-(dm.NdmPad-1)/2:(dm.NdmPad-1)/2)*dm.dx_dm; % meters, coords for the full DM arrays. Origin is centered between pixels for an even-sized array
end
dm.y_pupPad = dm.x_pupPad;

%% DM: (use NboxPad-sized postage stamps,

%--Find the locations of the postage stamps arrays in the larger pupilPad array
dm.xy_cent_act_inPix = dm.xy_cent_act*(dm.dm_spacing/dm.dx_dm); % Convert units to pupil-file pixels
if(~strcmpi(dm.centering,'pixel')  ) 
    error('falco_fpm_inf_cube_3x3.m: Not adapted for non-pixel centering.');
end
dm.xy_cent_act_box = round(dm.xy_cent_act_inPix); % Center locations of the postage stamps (in between pixels), in actuator widths
dm.xy_cent_act_box_inM = dm.xy_cent_act_box*dm.dx_dm; % now in meters 
if(mod(Nbox,2)==0)
    dm.xy_box_lowerLeft = dm.xy_cent_act_box + (dm.NdmPad-Nbox)/2 + 1; % indices of pixel in lower left of the postage stamp within the whole pupilPad array
else
    dm.xy_box_lowerLeft = dm.xy_cent_act_box + (dm.NdmPad)/2 -floor(Nbox/2) + 1; % indices of pixel in lower left of the postage stamp within the whole pupilPad array
end
%--Starting coordinates (in actuator widths) for updated influence function. This is
% interpixel centered, so do not translate!
dm.x_box0 = (-(Nbox-1)/2:(Nbox-1)/2)*dm.dx_dm;
[dm.Xbox0,dm.Ybox0] = meshgrid(dm.x_box0); %--meters, interpixel-centered coordinates for the master influence function

%--Limit the actuators used to those within 1 actuator width of the pupil
r_cent_act_box_inM = sqrt(dm.xy_cent_act_box_inM(1,:).^2 + dm.xy_cent_act_box_inM(2,:).^2);
%--Compute and store all the influence functions:
dm.inf_datacube = zeros(Nbox,Nbox,dm.NactTotal); %--initialize array of influence function "postage stamps"
dm.act_ele = []; % Indices of nonzero-ed actuators
for iact=1:dm.NactTotal
   dm.act_ele = [dm.act_ele; iact]; % Add actuator index to the keeper list
   dm.inf_datacube(:,:,iact) = dm.inf0;
end
fprintf('done.\n');

end %--END OF FUNCTION