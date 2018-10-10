% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to set variables to default values if they are not already defined.
%
%
% New or changed variables in going from HLC to SPHLC:
%  mp.F3.full.res --> mp.F3.full.in.res, mp.F3.full.out.res
%  mp.F3.compact.res --> mp.F3.compact.in.res, mp.F3.compact.out.res
%
%
% REVISION HISTORY:
% ----------------
% Modified on 2018-05-29 by A.J. Riggs from the HLC defaults to the SPHLC
%   defaults.
% Modified on 2018-03-22 by A.J. Riggs to have default values that can be
%   overwritten if the variable is already defined.
% Created on 2017-10-31 by A.J. Riggs.
% Modified on 2018-01-08 by A.J. Riggs to be vortex specific and move key
%   parameters outside this function.
% Created on 2017-10-31 by A.J. Riggs.

function mp = falco_config_defaults_SPHLC(mp)

%--Initialize some structures if they don't already exist
mp.P1.full.dummy = 1;
mp.P1.compact.dummy = 1;

mp.P2.dummy = 1;

mp.dm1.dummy = 1;
mp.dm2.dummy = 1;
mp.dm1.dummy = 1;
mp.dm2.dummy = 1;

mp.P3.dummy = 1;

mp.F3.full.in.dummy = 1;
mp.F3.full.out.dummy = 1;
mp.F3.compact.in.dummy = 1;
mp.F3.compact.out.dummy = 1;

mp.P4.full.dummy = 1;
mp.P4.compact.dummy = 1;

mp.F4.compact.dummy = 1;
mp.F4.full.dummy = 1;
mp.F4.corr.dummy = 1;
mp.F4.score.dummy = 1;


%%
%%--Record Keeping
if(isfield(mp,'SeriesNum')==false); mp.SeriesNum = 867; end %--Use the same Series # for sets of similar trials.
if(isfield(mp,'TrialNum')==false); mp.TrialNum = 5309; end%--Always use a diffrent Trial # for different calls of FALCO.

%%--WFSC. Iterations and Control Matrix Relinearization
if(isfield(mp,'controller')==false); mp.controller = 'EFC'; end %  end% Controller options: 'EFC' for EFC as an empirical grid search over tuning parametrs, 'conEFC' for constrained EFC using CVX.
if(isfield(mp,'Nitr')==false); mp.Nitr = 10; end%--Number of estimation+control iterations to perform
if(isfield(mp,'relinItrVec')==false); mp.relinItrVec = 1:mp.Nitr; end %--Which correction iterations at which to re-compute the control Jacobian

%%--Special Computational Settings
if(isfield(mp,'flagParfor')==false); mp.flagParfor = false; end %--whether to use parfor for Jacobian calculation
if(isfield(mp,'useGPU')==false); mp.useGPU = false; end %--whether to use GPUs for Jacobian calculation

%%--Coronagraph and Pupil Type
if(isfield(mp,'coro')==false); mp.coro = 'HLC';  end %--Tested Options: 'HLC','SPLC','SPHLC'
if(isfield(mp,'flagApod')==false); mp.flagApod = false;  end %--Can be any mask, even just an iris 
if(isfield(mp,'whichPupil')==false); mp.whichPupil = 'LUVOIR_B_offaxis'; end %'Simple';

%%--General:
if(isfield(mp,'centering')==false); mp.centering = 'pixel'; end %--Centering on the arrays at each plane: pixel or interpixel
if(isfield(mp,'fl')==false); mp.fl = 1;  end % meters. Arbitrary value chose for this design configuration. Keep as 1. Used for all focal lengths to keep magnification at each pupil at 1x. 

%%--Pupil Plane and DM Plane Properties
if(isfield(mp.dm1,'Nact')==false); mp.dm1.Nact = 48; end  % number of actuators across DM1
if(isfield(mp.dm2,'Nact')==false); mp.dm2.Nact = 48; end  % number of actuators across DM2
if(isfield(mp,'d_P2_dm1')==false); mp.d_P2_dm1 = 0; end  % distance (along +z axis) from P2 pupil to DM1 (meters)
if(isfield(mp,'d_dm1_dm2')==false); mp.d_dm1_dm2 = 1; end  % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
if(isfield(mp,'lambda0')==false); mp.lambda0 = 500e-9; end  % central wavelength of bandpass (meters)
if(isfield(mp,'fracBW')==false); mp.fracBW = 0.01; end   % fractional bandwidth of correction (Delta lambda / lambda)
if(isfield(mp,'Nsbp')==false); mp.Nsbp = 1; end  % number of sub-bandpasses across correction band 
if(isfield(mp,'Nwpsbp')==false); mp.Nwpsbp = 1; end % number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.


%% Controller Settings

%--Voltage range restrictions
if(isfield(mp.dm1,'maxAbsV')==false); mp.dm1.maxAbsV = 250.; end 
if(isfield(mp.dm2,'maxAbsV')==false); mp.dm2.maxAbsV = 250.; end 

%% Make new PSD maps for full model (NOT USED FOR DESIGN-ONLY CODE)
if(isfield(mp,'flagNewPSD')==false); mp.flagNewPSD = false; end  %--to make new PSD maps

%% Whether to include planet in the images
if(isfield(mp,'planetFlag')==false); mp.planetFlag = false; end 

%% Throughput calculations
if(isfield(mp,'thput_radius')==false); mp.thput_radius = 0.7; end  %--lambda_c/D, MAX radius to score core throughput
if(isfield(mp,'thput_eval_x')==false); mp.thput_eval_x = 6; end  %--lambda_c/D, x-axis value at which to evaluate throughput
if(isfield(mp,'thput_eval_y')==false); mp.thput_eval_y = 0; end  %--lambda_c/D, x-axis value at which to evaluate throughput

%% Coronagraphic Mask Properties:

%%--Pupil Masks
mp = falco_config_load_pupil_defaults(mp);

%%--Apodizer (Shaped Pupil) Properties (Plane P3). Also contains FPM and
%   Lyot stop specifications for that apodized coronagraph design. 
mp = falco_config_load_apodizer_defaults(mp);

%--FPM parameters  
if(isfield(mp.F3.full.in,'res')==false);    mp.F3.full.in.res = 30; end % sampling of inner FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.full.out,'res')==false);   mp.F3.full.out.res = 4; end % sampling of outer FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.compact.in,'res')==false);  mp.F3.compact.in.res = 30; end % sampling of inner FPM for compact model, in pixels per lambda0/D
if(isfield(mp.F3.compact.out,'res')==false); mp.F3.compact.out.res = 4; end % sampling of outer FPM for compact model, in pixels per lambda0/D
% if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 30; end % sampling of FPM for full model, in pixels per lambda0/D
% if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 30; end % sampling of FPM for compact model, in pixels per lambda0/D

%--FPM shape (if not already defined in falco_config_load_apodizer_defaults)
if(isfield(mp.F3,'Rin')==false); mp.F3.Rin = 2.8; end % inner hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = 20; end % outer hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees

%%--Final Focal Plane (F4) Properties
mp = falco_config_load_F4_defaults(mp,DM);

%% Controller weights and thresholds

%--Spatial pixel weighting of the Control Jacobian
if(isfield(mp,'WspatialDef')==false); mp.WspatialDef = []; end %[mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1]; end   %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--Threshold level for culling actuators from the Jacobian
if(isfield(mp,'logGmin')==false); mp.logGmin = -6; end  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators


%% Deformable Mirror (DM) Parameters

if(isfield(mp,'dm_ind')==false); mp.dm_ind = [1 2]; end% vector of which DMs to use for control.
if(isfield(mp,'dm_weights')==false); mp.dm_weights = ones(9,1);  end % vector of relative weighting of DMs for EFC

%--DM1 parameters
% if(isfield(mp.dm1,'Nact')==false); mp.dm1.Nact = 48; end  % number of actuators across DM1
if(isfield(mp.dm1,'VtoH')==false); mp.dm1.VtoH = 1*1e-9*ones(mp.dm1.Nact); end  % Gains: volts to meters in surface height;
if(isfield(mp.dm1,'xtilt')==false); mp.dm1.xtilt = 0; end 
if(isfield(mp.dm1,'ytilt')==false); mp.dm1.ytilt = 0; end 
if(isfield(mp.dm1,'zrot')==false); mp.dm1.zrot = 0; end  %--clocking angle (degrees)
if(isfield(mp.dm1,'xc')==false); mp.dm1.xc = (mp.dm1.Nact/2 - 1/2); end  % x-center of DM in actuator widths
if(isfield(mp.dm1,'yc')==false); mp.dm1.yc = (mp.dm1.Nact/2 - 1/2); end  % x-center of DM in actuator widths
if(isfield(mp.dm1,'edgeBuffer')==false); mp.dm1.edgeBuffer = 1; end  % Radius (in actuator spacings) outside of pupil to compute influence functions for.

%--DM2 parameters
% if(isfield(mp.dm2,'Nact')==false); mp.dm2.Nact = 48; end  % number of actuators across DM1
if(isfield(mp.dm2,'VtoH')==false); mp.dm2.VtoH = 1*1e-9*ones(mp.dm2.Nact); end  % Gains: volts to meters in surface height;
if(isfield(mp.dm2,'xtilt')==false); mp.dm2.xtilt = 0; end 
if(isfield(mp.dm2,'ytilt')==false); mp.dm2.ytilt = 0; end 
if(isfield(mp.dm2,'zrot')==false); mp.dm2.zrot = 0; end  %--clocking angle (degrees)
if(isfield(mp.dm2,'xc')==false); mp.dm2.xc = mp.dm2.Nact/2 - 1/2; end  % x-center of DM in actuator widths
if(isfield(mp.dm2,'yc')==false); mp.dm2.yc = mp.dm2.Nact/2 - 1/2; end  % x-center of DM in actuator widths
if(isfield(mp.dm2,'edgeBuffer')==false); mp.dm2.edgeBuffer = 1; end  % Radius (in actuator spacings) outside of pupil to compute influence functions for.

%--DM Actuator characteristics
if(isfield(mp.dm1,'dx_inf0')==false); mp.dm1.dx_inf0 = 1e-4; end  % meters, sampling of the influence function;
if(isfield(mp.dm1,'dm_spacing')==false); mp.dm1.dm_spacing = 1e-3; end  % meters, pitch of DM actuators
if(isfield(mp.dm1,'inf0')==false); mp.dm1.inf0 = 1*fitsread('influence_dm5v2.fits'); end     %  -1*fitsread('inf64.3.fits');                              
if(isfield(mp.dm2,'dx_inf0')==false); mp.dm2.dx_inf0 = 1e-4; end  % meters, sampling of the influence function;
if(isfield(mp.dm2,'dm_spacing')==false); mp.dm2.dm_spacing = 1e-3; end %0.9906e-3; % meters, pitch of DM actuators
if(isfield(mp.dm2,'inf0')==false); mp.dm2.inf0 = 1*fitsread('influence_dm5v2.fits'); end     

%--Aperture stops at DMs
if(isfield(mp,'flagDM1stop')==false); mp.flagDM1stop = false; end  %--logical flag whether to include the stop at DM1 or not
if(isfield(mp,'flagDM2stop')==false); mp.flagDM2stop = false; end  %--logical flag whether to include the stop at DM2 or not
if(isfield(mp.dm1,'Dstop')==false); mp.dm1.Dstop = mp.dm1.Nact*1e-3; end  %--diameter of circular stop at DM1 and centered on the beam
if(isfield(mp.dm2,'Dstop')==false); mp.dm2.Dstop = mp.dm2.Nact*1e-3; end  %--diameter of circular stop at DM2 and centered on the beam


%--DM9: Phase control at the FPM
% % mp.dm9.Nrings = 25; % Number of radial hexagonal rings (including the center) for the FPM's phase "DM"
% mp.dm8.VtoHavg = 1e-9; %--meters, gain of each "actuator" for the FPM nickel.
% mp.dm9.VtoHavg = 1e-9; %--meters, gain of each "actuator" for the FPM dielectric.

%--DM9 parameters
if(isfield(mp.dm9,'flag_hex_array')==false); mp.dm9.flag_hex_array = false; end %--true->hex grid. false->square grid
if(isfield(mp.dm9,'actres')==false);  mp.dm9.actres = 10; end % number of "actuators" per lambda0/D in the FPM's focal plane. Only used for square grid

if(isfield(mp.dm8,'Vmin')==false);  mp.dm8.Vmin = 0; end % minimum thickness of FPM metal layer (nm)
if(isfield(mp.dm8,'Vmax')==false);  mp.dm8.Vmax = 300; end % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)
if(isfield(mp.dm9,'Vmin')==false);  mp.dm9.Vmin = 0; end % minimum thickness of FPM dielectric layer (nm)
if(isfield(mp.dm9,'Vmax')==false);  mp.dm9.Vmax = 500; end % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

mp.dm9.Nact = ceil_even(2*mp.F3.Rin*mp.dm9.actres); % number of actuators across DM9 (if not in a hex grid)
mp.dm9.xtilt = 0;
mp.dm9.ytilt = 0;
mp.dm9.zrot = 0; %--clocking angle (degrees)
% mp.dm9.xc = mp.dm9.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% mp.dm9.yc = mp.dm9.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
mp.dm9.edgeBuffer = 1;%0; %1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.

%--DM9 Actuator characteristics
mp.dm9.dx_inf0_act = 1/10; % units of actuator spacings, sampling of the influence function (for the Xinetics influence function, always 1/10 of the actuator pitch).
mp.dm9.inf0 = 1*fitsread('influence_dm5v2.fits');    %  -1*fitsread('inf64.3.fits');                              
%--Cropped influence function for FPM phase
%mp.dm9.Ncrop_inf0 = 45;%20;%45;  %--Half the amount to extract from the center of the 91x91 pixel DM influence function. >=45 gives the uncropped influence functionNhalf = ceil(length(mp.dm9.inf0)/2);
%mp.dm9.FPMbuffer = 0.2;%0.25; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)

%%--Pupil Plane Properties
if(isfield(mp.P2,'D')==false);  mp.P2.D = (mp.dm1.Nact-2)*1e-3; end %46.3e-3; end % beam diameter at pupil closest to the DMs  (meters)
if(isfield(mp.P3,'D')==false);  mp.P3.D = mp.P2.D; end  % beam diameter at apodizer plane pupil (meters).
if(isfield(mp.P4,'D')==false);  mp.P4.D = mp.P2.D; end  % beam diameter at Lyot plane pupil (meters).


%% Final Focal Plane (F4) Properties


%--Specs for Correction (Corr) region and the Scoring (Score) region.
% if(isfield(mp.F4.corr,'Rin')==false); mp.F4.corr.Rin  = mp.F3.Rin; end  %--lambda0/D, inner radius of correction region
if(isfield(mp.F4.corr,'Rin')==false)
    if(isfield(mp.F3,'RinA')) 
        mp.F4.corr.Rin  = mp.F3.RinA;
    else
        mp.F4.corr.Rin  = mp.F3.Rin;
    end
end  %--lambda0/D, inner radius of correction region
if(isfield(mp.F4.score,'Rin')==false); mp.F4.score.Rin = mp.F4.corr.Rin; end  %--Needs to be >= that of Correction mask
if(isfield(mp.F4.corr,'Rout')==false); mp.F4.corr.Rout  = min( [floor(mp.dm1.Nact/2*(1-mp.fracBW/2)), mp.F3.Rout ] ); end  %--lambda0/D, outer radius of correction region
if(isfield(mp.F4.score,'Rout')==false); mp.F4.score.Rout = mp.F4.corr.Rout; end  %--Needs to be <= that of Correction mask

if(isfield(mp.F4.corr,'ang')==false); mp.F4.corr.ang  = 180; end  %--degrees per side
if(isfield(mp.F4.score,'ang')==false); mp.F4.score.ang = 180; end  %--degrees per side
if(isfield(mp.F4,'sides')==false); mp.F4.sides = 'both'; end  %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%%--Final Focal Plane (F4) Properties
if(isfield(mp.F4.compact,'res')==false); mp.F4.res = 3; end  %--Pixels per lambda_c/D
if(isfield(mp.F4.full,'res')==false); mp.F4.full.res = 6; end  %--Pixels per lambda_c/D
if(isfield(mp.F4,'FOV')==false); mp.F4.FOV = 1 + mp.F4.corr.Rout; end  % minimum desired field of view (along both axes) in lambda0/D

end %--END OF FUNCTION



