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
% dm.xc = dm.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% dm.yc = dm.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
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
% -Clean up uses of dx_dm vs dx_dm to avoid any possible differences
% since they need to be the same value.



%-------------------------------------------------------------------
function dm = falco_gen_dm_poke_cube(dm,mp,dx_dm,varargin)




%--Enable the ability to turn off the data cube calculation with a keyword.
flagGenCube = true;
icav = 0;             % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'nocube'}
        flagGenCube  = false; %--Don't compute the data cube of influence functions. (For the full model.)
      otherwise
        error('falco_gen_dm_poke_cube: Unknown keyword: %s\n', varargin{icav});
    end
end

%--Define this flag if it doesn't exist in the older code for square actuator arrays only
if(isfield(dm,'flag_hex_array')==false) 
    dm.flag_hex_array=false;
end

%--Set the order of operations
orderOfOps = 'XYZ';
zyx = false;
if(isfield(dm,'flagZYX')); 
    if(dm.flagZYX)
        orderOfOps = 'ZYX'; 
        zyx = true;
    end
end

    
%--Compute sampling of the pupil. Assume that it is square.
dm.dx_dm = dx_dm;
dm.dx = dx_dm;

%--Default to being centered on a pixel (FFT-style) if no centering is specified
if(isfield(dm,'centering'))
    %dm.centering = dm.centering;
else
    dm.centering = 'pixel'; %--Centered on a pixel (default if not specified)
end

%--Compute coordinates of original influence function
Ninf0 = length(dm.inf0); %--Number of points across the influence function at its native resolution
x_inf0 = (-(Ninf0-1)/2:(Ninf0-1)/2)*dm.dx_inf0; % True for even- or odd-sized influence function maps as long as they are centered on the array.

Ndm0 = ceil_even( Ninf0 + (dm.Nact - 1)*(dm.dm_spacing/dm.dx_inf0) ); %--Number of points across the DM surface at native influence function resolution
dm.NdmMin = ceil_even( Ndm0*(dm.dx_inf0/dm.dx))+2; %--Number of points across the (un-rotated) DM surface at new, desired resolution.
dm.Ndm = ceil_even( max(abs([sqrt(2)*cosd(45-dm.zrot),sqrt(2)*sind(45-dm.zrot)]))*Ndm0*(dm.dx_inf0/dm.dx))+2; %--Number of points across the array to fully contain the DM surface at new, desired resolution and z-rotation angle.

[Xinf0,Yinf0] = meshgrid(x_inf0);



%--Compute list of initial actuator center coordinates (in actutor widths).
if(dm.flag_hex_array) %--Hexagonal, hex-packed grid
    Nrings = dm.Nrings;
    x_vec = [];
    y_vec = [];
    % row number (rowNum) is 1 for the center row and 2 is above it, etc.
    % Nacross is the total number of segments across that row
    for rowNum = 1:Nrings
        Nacross = 2*Nrings - rowNum; % Number of actuators across at that row (for hex tiling in a hex shape)
        yval = sqrt(3)/2*(rowNum-1);
        bx = Nrings - (rowNum+1)/2; % x offset from origin

        xs = (0:Nacross-1).' - bx; % x values are 1 apart
        ys = yval*ones(Nacross,1); % same y-value for the entire row

        if(rowNum==1)
            x_vec = [x_vec;xs];
            y_vec = [y_vec;ys]; 
        else
            x_vec = [x_vec;xs;xs];
            y_vec = [y_vec;ys;-ys]; % rows +/-n have +/- y coordinates
        end
    end    
else %--Square grid
    [dm.Xact,dm.Yact] = meshgrid((0:dm.Nact-1)-dm.xc,(0:dm.Nact-1)-dm.yc); % in actuator widths
    x_vec = dm.Xact(:);
    y_vec = dm.Yact(:);
end
dm.NactTotal = length(x_vec); %--Total number of actuators in the 2-D array



tlt  = zeros(1, 3);
tlt(1) = dm.xtilt;
tlt(2) = dm.ytilt;
tlt(3) = -dm.zrot;

sa   = sind(tlt(1));
ca   = cosd(tlt(1));
sb   = sind(tlt(2));
cb   = cosd(tlt(2));
sg   = sind(tlt(3));
cg   = cosd(tlt(3));

if zyx == true
    Mrot = [               cb * cg,               -cb * sg,       sb, 0.0; ...
            ca * sg + sa * sb * cg, ca * cg - sa * sb * sg, -sa * cb, 0.0; ...
            sa * sg - ca * sb * cg, sa * cg + ca * sb * sg,  ca * cb, 0.0; ...
                               0.0,                    0.0,      0.0, 1.0];
else
    Mrot = [ cb * cg, sa * sb * cg - ca * sg, ca * sb * cg + sa * sg, 0.0; ...
             cb * sg, sa * sb * sg + ca * cg, ca * sb * sg - sa * cg, 0.0; ...
            -sb,      sa * cb,                ca * cb,                0.0; ...
                 0.0,                    0.0,                    0.0, 1.0];
end

for iact=1:dm.NactTotal
    xyzVals = [x_vec(iact); y_vec(iact); 0; 1];
    xyzValsRot = Mrot*xyzVals;
    dm.xy_cent_act(:,iact) = xyzValsRot(1:2);
end


  
N0 = max(size(dm.inf0));
Npad = ceil_odd( sqrt(2)*max(size(dm.inf0)) );
inf0pad = zeros(Npad,Npad);

inf0pad( ceil(Npad/2)-floor(N0/2):ceil(Npad/2)+floor(N0/2), ceil(Npad/2)-floor(N0/2):ceil(Npad/2)+floor(N0/2) ) = dm.inf0;

[ydim,xdim] = size(inf0pad);

xd2  = fix(xdim / 2) + 1;
yd2  = fix(ydim / 2) + 1;
cx   = ([1 : xdim] - xd2) ;
cy   = ([1 : ydim] - yd2) ;
[Xs0, Ys0] = meshgrid(cx, cy);

xsNew = 0*Xs0;
ysNew = 0*Ys0;
 
for ii=1:numel(Xs0)
    xyzVals = [Xs0(ii); Ys0(ii); 0; 1];
    xyzValsRot = Mrot*xyzVals;
    xsNew(ii) = xyzValsRot(1);
    ysNew(ii) = xyzValsRot(2);
end

% Calculate the interpolated DM grid (set extrapolated values to 0.0)
dm.infMaster = griddata(xsNew,ysNew,inf0pad,Xs0,Ys0,'cubic');%,'cubic',0);
dm.infMaster(isnan(dm.infMaster)) = 0;
 
% x_inf0 = (-(Npad-1)/2:(Npad-1)/2)*dm.dx_inf0; % True for even- or odd-sized influence function maps as long as they are centered on the array.
% [Xinf0,Yinf0] = meshgrid(x_inf0);
 



%--Crop down the influence function until it has no zero padding left
infSum = sum(dm.infMaster(:));
infDiff = 0; counter = 0;
while( abs(infDiff) <= 1e-7)
    counter = counter + 2;
    %Ninf0pad = length(dm.infMaster)-counter; %--Number of points across the rotated, cropped-down influence function at the original resolution
    infDiff = infSum - sum(sum( abs(dm.infMaster(1+counter/2:end-counter/2,1+counter/2:end-counter/2)) )); %--Subtract an extra 2 to negate the extra step that overshoots.
end
counter = counter - 2;
Ninf0pad = length(dm.infMaster)-counter; %Ninf0pad = Ninf0pad+2;
infMaster2 = dm.infMaster(1+counter/2:end-counter/2,1+counter/2:end-counter/2); % padOrCropEven(dm.infMaster,Ncrop); %--The cropped-down Lyot stop for the compact model       
% figure; imagesc(log10(abs(infMaster2)));


dm.infMaster = infMaster2;
Npad = Ninf0pad;

x_inf0 = (-(Npad-1)/2:(Npad-1)/2)*dm.dx_inf0; % True for even- or odd-sized influence function maps as long as they are centered on the array.
[Xinf0,Yinf0] = meshgrid(x_inf0);



%%%%%%%%%%%%%%%%%%%%%%%


% %--Apply x- and y-projections and then z-rotation to the original influence
% %    function to make a master influence function.
% dm.Xrot = Xinf0/cosd(dm.ytilt); % Divide coords by projection factor to squeeze the influence function
% dm.Yrot = Yinf0/cosd(dm.xtilt);
% dm.infMaster = interp2(Xinf0,Yinf0,dm.inf0,dm.Xrot,dm.Yrot,'cubic',0);
% dm.infMaster = imrotate(dm.infMaster,dm.zrot,'bicubic','crop');

%--Compute the size of the postage stamps.
Nbox = ceil_even(((Ninf0pad*dm.dx_inf0)/dx_dm)); % Number of points across the influence function array at the DM plane's resolution. Want as even
dm.Nbox = Nbox;
%--Also compute their padded sizes for the angular spectrum (AS) propagation between P2 and DM1 or between DM1 and DM2
Nmin = ceil_even( max(mp.sbp_centers)*max(abs([mp.d_P2_dm1, mp.d_dm1_dm2,(mp.d_P2_dm1+mp.d_dm1_dm2)]))/dx_dm^2 ); %--Minimum number of points across for accurate angular spectrum propagation
% dm.NboxAS = 2^(nextpow2(max([Nbox,Nmin])));  %--Zero-pad for FFTs in A.S. propagation. Use a larger array if the max sampling criterion for angular spectrum propagation is violated
dm.NboxAS = max([Nbox,Nmin]);  %--Use a larger array if the max sampling criterion for angular spectrum propagation is violated

% dm.NdmPad = ceil_even( dm.Ndm + (dm.NboxAS-dm.Nbox) ); %--Number of points across the DM surface (with padding for angular spectrum propagation) at new, desired resolution.

% if( Nbox < Nmin ) %--Use a larger array if the max sampling criterion for angular spectrum propagation is violated
%     dm.NboxAS = 2^(nextpow2(Nmin)); %2*ceil(1/2*min(mp.sbp_centers)*mp.d_dm1_dm2/dx_dm^2);
% else
%     dm.NboxAS = 2^(nextpow2(Nbox));
% end

%% Pad the pupil to at least the size of the DM(s) surface(s) to allow all actuators to be located outside the pupil.
% (Same for both DMs)

%--Find actuator farthest from center:
dm.r_cent_act = sqrt(dm.xy_cent_act(1,:).^2 + dm.xy_cent_act(2,:).^2);
dm.rmax = max(abs(dm.r_cent_act));
NpixPerAct = dm.dm_spacing/dx_dm;
if(dm.flag_hex_array)
    %dm.NdmPad = 2*ceil(1/2*Nbox*2) + 2*ceil((1/2*2*(dm.rmax)*dm.dx_inf0_act)*Nbox); %2*ceil((dm.rmax+3)*dm.dm_spacing/Dpup*Npup);
    dm.NdmPad = ceil_even((2*(dm.rmax+2))*NpixPerAct + 1); % padded 2 actuators past the last actuator center to avoid trying to index outside the array 
else
    dm.NdmPad = ceil_even( ( dm.NboxAS + 2*(1+ (max(max(abs(dm.xy_cent_act)))+0.5)*NpixPerAct)) ); % DM surface array padded by the width of the padded influence function to prevent indexing outside the array. The 1/2 term is because the farthest actuator center is still half an actuator away from the nominal array edge. 
end

%--Compute coordinates (in meters) of the full DM array
if(strcmpi(dm.centering,'pixel')  ) 
    dm.x_pupPad = (-dm.NdmPad/2:(dm.NdmPad/2 - 1))*dx_dm; % meters, coords for the full DM arrays. Origin is centered on a pixel
else
    dm.x_pupPad = (-(dm.NdmPad-1)/2:(dm.NdmPad-1)/2)*dx_dm; % meters, coords for the full DM arrays. Origin is centered between pixels for an even-sized array
end
dm.y_pupPad = dm.x_pupPad;



%% DM: (use NboxPad-sized postage stamps)

if(flagGenCube)
    if(dm.flag_hex_array==false)
        fprintf('  Influence function padded from %d to %d points for A.S. propagation.\n',Nbox,dm.NboxAS);
        %fprintf('  Influence function padded to 2^nextpow2(%d) = %d for A.S. propagation.\n',2*ceil(1/2*max([Nbox,Nmin])),dm.NboxAS);
    end
    tic
    fprintf('Computing datacube of DM influence functions... ');

    %--Find the locations of the postage stamps arrays in the larger pupilPad array
    dm.xy_cent_act_inPix = dm.xy_cent_act*(dm.dm_spacing/dx_dm); % Convert units to pupil-file pixels
%     if(strcmpi(dm.centering,'pixel')  ) 
       dm.xy_cent_act_inPix = dm.xy_cent_act_inPix + 0.5; %--For the half-pixel offset if pixel centered. 
%     end
    dm.xy_cent_act_box = round(dm.xy_cent_act_inPix); % Center locations of the postage stamps (in between pixels), in actuator widths
    dm.xy_cent_act_box_inM = dm.xy_cent_act_box*dx_dm; % now in meters 
    dm.xy_box_lowerLeft = dm.xy_cent_act_box + (dm.NdmPad-Nbox)/2 + 1; % indices of pixel in lower left of the postage stamp within the whole pupilPad array

    %--Starting coordinates (in actuator widths) for updated influence function. This is
    % interpixel centered, so do not translate!
    dm.x_box0 = (-(Nbox-1)/2:(Nbox-1)/2)*dx_dm;
    [dm.Xbox0,dm.Ybox0] = meshgrid(dm.x_box0); %--meters, interpixel-centered coordinates for the master influence function

    %--Limit the actuators used to those within 1 actuator width of the pupil
    r_cent_act_box_inM = sqrt(dm.xy_cent_act_box_inM(1,:).^2 + dm.xy_cent_act_box_inM(2,:).^2);
    %--Compute and store all the influence functions:
    dm.inf_datacube = zeros(Nbox,Nbox,dm.NactTotal);%dm.Nact^2); %--initialize array of influence function "postage stamps"
    dm.act_ele = []; % Indices of nonzero-ed actuators
    for iact=1:dm.NactTotal %dm.Nact^2
%         if(r_cent_act_box_inM(iact) < D/2 + dm.edgeBuffer*Nbox*dx_dm) %--Don't use actuators too far outside the beam
            dm.act_ele = [dm.act_ele; iact]; % Add actuator index to the keeper list
            dm.Xbox = dm.Xbox0 - (dm.xy_cent_act_inPix(1,iact)-dm.xy_cent_act_box(1,iact))*dx_dm; % X = X0 -(x_true_center-x_box_center)
            dm.Ybox = dm.Ybox0 - (dm.xy_cent_act_inPix(2,iact)-dm.xy_cent_act_box(2,iact))*dx_dm; % Y = Y0 -(y_true_center-y_box_center)
            dm.inf_datacube(:,:,iact) = interp2(Xinf0,Yinf0,dm.infMaster,dm.Xbox,dm.Ybox,'spline',0);
%         end
    end
    
    fprintf('done.  Time = %.1fs\n',toc);

else
    dm.act_ele = (1:dm.NactTotal).';    
end

end %--END OF FUNCTION































