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
%
%--Inputs
%- distance from pupil that the actuators aren't computed. --> Move this
%outside the function and let that be the user's choice which ones to zero.
%Instead, just output the 

% %--Load the beam diameter and size
% Npup = length(pupil);
% Dpup = 46.3e-3; % meters, distance across pupil file
% 
% %--DM parameters
% dm.Nact = 48; % number of actuators across DM1
% dm.VtoH = 1*1e-9*ones(dm.Nact); % Gains: volts to meters in surface height;
% 
% dm.xtilt = 0;
% dm.ytilt = 0;
% dm.zrot = 0;
% dm.xcent_dm = dm.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% dm.ycent_dm = dm.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% 
% %--DM Actuator characteristics
% dm.dx_inf0 = 1e-4; % meters, sampling of the influence function;
% dm.dm_spacing = 0.9906e-3; % meters, pitch of DM actuators
% dm.inf0 = -1*fitsread('influence_dm5v2.fits');s
%
%
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
function dm = falco_fpm_inf_cube(dm)

fprintf('Computing datacube of FPM influence functions... ');

% if(isfield(dm,'flag_hex_array')==false) %--Define this flag if it doesn't exist (in older code for square actuator arrays only)
%     dm.flag_hex_array=false;
% end
    
%--Compute sampling of the pupil. Assume that it is square.
dm.dx_dm = dm.dxi;
% dm.dx_dm = D/N;
% dm.dx_dm = dx_dm;
% dy_pup = dx_pup;

%--Default to being centered on a pixel (FFT-style) if no centering is specified
if(isfield(dm,'centering'))
    %dm.centering = dm.centering;
else
    dm.centering = 'pixel'; %--Centered on a pixel (default if not specified)
end

%--Compute coordinates of original influence function
Ninf0 = length(dm.inf0);
% if(mod(Ninf0,2)==0 && ( strcmpi(dm.centering,'pixel') ) ) %--if influence function is in an even-sized array and pixel-centered
%     x_inf0 = (-Ninf0/2:(Ninf0/2 - 1))*dm.dx_inf0; 
% else
x_inf0 = (-(Ninf0-1)/2:(Ninf0-1)/2)*dm.dx_inf0; % True for even- or odd-sized influence function maps as long as they are centered on the array.
% end
[Xinf0,Yinf0] = meshgrid(x_inf0);


%--Compute list of initial actuator center coordinates (in actutor widths).
% if(dm.flag_hex_array) %--Hexagonal, hex-packed grid
%     Nrings = dm.Nrings;
%     x_vec = [];
%     y_vec = [];
%     % row number (rowNum) is 1 for the center row and 2 is above it, etc.
%     % Nacross is the total number of segments across that row
%     for rowNum = 1:Nrings
%         Nacross = 2*Nrings - rowNum; % Number of actuators across at that row (for hex tiling in a hex shape)
%         yval = sqrt(3)/2*(rowNum-1);
%         bx = Nrings - (rowNum+1)/2; % x offset from origin
% 
%         xs = (0:Nacross-1).' - bx; % x values are 1 apart
%         ys = yval*ones(Nacross,1); % same y-value for the entire row
% 
%         if(rowNum==1)
%             x_vec = [x_vec;xs];
%             y_vec = [y_vec;ys]; 
%         else
%             x_vec = [x_vec;xs;xs];
%             y_vec = [y_vec;ys;-ys]; % rows +/-n have +/- y coordinates
%         end
%     end    
% else %--Square grid
    [dm.Xact,dm.Yact] = meshgrid((0:dm.Nact-1)-dm.xcent_dm,(0:dm.Nact-1)-dm.ycent_dm); % in actuator widths
    x_vec = dm.Xact(:);
    y_vec = dm.Yact(:);
% end
dm.NactTotal = length(x_vec); %--Total number of actuators in the 2-D array

dm.xy_cent_act  = [x_vec*cosd(dm.ytilt),y_vec*cosd(dm.xtilt)].'; % in actuator widths, Perform x/y projections because of foreshortening
% dm.xy_cent_act = [dm.Xact(:)*cosd(dm.xtilt),dm.Yact(:)*cosd(dm.ytilt)].'; % in actuator widths, Perform x/y projections because of foreshortening
rotMat = [cosd(-dm.zrot),-sind(-dm.zrot); ...
          sind(-dm.zrot),cosd(-dm.zrot)]; % rotation matrix for z-axis
for iact=1:dm.NactTotal %dm.Nact^2
    dm.xy_cent_act(:,iact) = rotMat*dm.xy_cent_act(:,iact); % Rotate the coordinates around z-axis
end

%--Apply x- and y-projections and then z-rotation to the original influence
%    function to make a master influence function.
dm.Xrot = Xinf0/cosd(dm.ytilt); % Divide coords by projection factor to squeeze the influence function
dm.Yrot = Yinf0/cosd(dm.xtilt);
dm.infMaster = interp2(Xinf0,Yinf0,dm.inf0,dm.Xrot,dm.Yrot,'cubic',0);
if(dm.zrot~=0)
    dm.infMaster = imrotate(dm.infMaster,dm.zrot,'bicubic','crop');
end
%--Compute the size of the postage stamps.
Nbox = ceil_even(((Ninf0*dm.dx_inf0)/dm.dx_dm)); % Number of points across the influence function in the pupil file's spacing. Want as even
dm.Nbox = Nbox;
fprintf('FPM influence function size =\t%dx%d ',Nbox,Nbox);



%% Pad the pupil to at least the size of the DM(s) surface(s) to allow all actuators to be located outside the pupil.
% (Same for both DMs)

%--Find actuator farthest from center:
dm.r_cent_act = sqrt(dm.xy_cent_act(1,:).^2 + dm.xy_cent_act(2,:).^2);
dm.rmax = max(abs(dm.r_cent_act));
dm.absxymax = max( abs(dm.xy_cent_act(:)) );
NpixPerAct = dm.dm_spacing/dm.dx_dm;
% if(dm.flag_hex_array)
%     %dm.NdmPad = 2*ceil(1/2*Nbox*2) + 2*ceil((1/2*2*(dm.rmax)*dm.dx_inf0_act)*Nbox); %2*ceil((dm.rmax+3)*dm.dm_spacing/Dpup*Npup);
%     dm.NdmPad = ceil_even((2*(dm.rmax+2))*NpixPerAct + 1); % padded 2 actuators past the last actuator center to avoid trying to index outside the array 
% else
% %     %dm.NdmPad = ceil_even( ( 2*(dm.rmax*NpixPerAct + 1)) ); % padded 1/2 an actuator past the farthest actuator center (on each side) to prevent indexing outside the array 
% %     %dm.NdmPad = ceil_even( ( 2*(dm.rmax*NpixPerAct) + 1) ); % padded 1 pixel past the farthest actuator center (on each side) to prevent indexing outside the array 
% %     dm.NdmPad = ceil_even( ( 2*(dm.absxymax*NpixPerAct) + Nbox+1) ); % padded 1/2 an actuator past the farthest actuator center (on each side) to prevent indexing outside the array 
%     dm.NdmPad = 2 + ceil_even( (dm.Nact+10)*NpixPerAct ); % prevent indexing outside the array 
% end
% % dm.NxiFPM = ceil_even( dm.Nact*NpixPerAct ); % Cropped down number of points for plotting

dm.NdmPad = 2 + ceil_even( (dm.Nact+10)*NpixPerAct ); % prevent indexing outside the array 

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
if(strcmpi(dm.centering,'pixel')  ) 
   dm.xy_cent_act_inPix = dm.xy_cent_act_inPix + 0.5; %--For the half-pixel offset if pixel centered. 
end
dm.xy_cent_act_box = round(dm.xy_cent_act_inPix); % Center locations of the postage stamps (in between pixels), in actuator widths
dm.xy_cent_act_box_inM = dm.xy_cent_act_box*dm.dx_dm; % now in meters 
dm.xy_box_lowerLeft = dm.xy_cent_act_box + (dm.NdmPad-Nbox)/2 + 1; % indices of pixel in lower left of the postage stamp within the whole pupilPad array

%--Starting coordinates (in actuator widths) for updated influence function. This is
% interpixel centered, so do not translate!
dm.x_box0 = (-(Nbox-1)/2:(Nbox-1)/2)*dm.dx_dm;
[dm.Xbox0,dm.Ybox0] = meshgrid(dm.x_box0); %--meters, interpixel-centered coordinates for the master influence function

%--Limit the actuators used to those within 1 actuator width of the pupil
r_cent_act_box_inM = sqrt(dm.xy_cent_act_box_inM(1,:).^2 + dm.xy_cent_act_box_inM(2,:).^2);
%--Compute and store all the influence functions:
dm.inf_datacube = zeros(Nbox,Nbox,dm.NactTotal);%dm.Nact^2); %--initialize array of influence function "postage stamps"
dm.act_ele = []; % Indices of nonzero-ed actuators
for iact=1:dm.NactTotal %dm.Nact^2
%     if(r_cent_act_box_inM(iact) < D/2 + dm.edgeBuffer*Nbox*dm.dx_dm) %--Don't use actuators too far outside the beam
        dm.act_ele = [dm.act_ele; iact]; % Add actuator index to the keeper list
        dm.Xbox = dm.Xbox0 - (dm.xy_cent_act_inPix(1,iact)-dm.xy_cent_act_box(1,iact))*dm.dx_dm; % X = X0 -(x_true_center-x_box_center)
        dm.Ybox = dm.Ybox0 - (dm.xy_cent_act_inPix(2,iact)-dm.xy_cent_act_box(2,iact))*dm.dx_dm; % Y = Y0 -(y_true_center-y_box_center)
        dm.inf_datacube(:,:,iact) = interp2(Xinf0,Yinf0,dm.infMaster,dm.Xbox,dm.Ybox,'cubic',0);
%     end
end
fprintf('done.\n');


end %--END OF FUNCTION































