% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Function to set variables to default values if they are not already defined.
%
%
% Modified from falco_config_defaults_LC on 2018-11-20 by G. Ruane
% Modified on 2018-03-22 by A.J. Riggs to have default values that can be
%   overwritten if the variable is already defined.
% Modified on 2018-01-08 by A.J. Riggs to be vortex specific and move key
%   parameters outside this function.
% Created on 2017-10-31 by A.J. Riggs.
%




function mp = falco_config_defaults_Roddier(mp)


%% Initialize some structures if they don't already exist
mp.P1.full.dummy = 1;
mp.P1.compact.dummy = 1;

mp.P2.dummy = 1;

mp.dm1.dummy = 1;
mp.dm2.dummy = 1;
mp.dm1.dummy = 1;
mp.dm2.dummy = 1;

mp.P3.dummy = 1;

mp.F3.full.dummy = 1;
mp.F3.compact.dummy = 1;

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
if(isfield(mp,'coro')==false); mp.coro = 'Roddier';  end %--Tested Options: 'HLC','SPLC','SPHLC'
if(isfield(mp,'flagApod')==false); mp.flagApod = false;  end %--Can be any mask, even just an iris 
if(isfield(mp,'whichPupil')==false); mp.whichPupil = 'Simple'; end %'Simple';

%%--General:
if(isfield(mp,'centering')==false); mp.centering = 'pixel'; end %--Centering on the arrays at each plane: pixel or interpixel


%%--Pupil Plane and DM Plane Properties
if(isfield(mp,'d_P2_dm1')==false); mp.d_P2_dm1 = 0; end  % distance (along +z axis) from P2 pupil to DM1 (meters)
if(isfield(mp,'d_dm1_dm2')==false); mp.d_dm1_dm2 = 1; end  % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
if(isfield(mp,'lambda0')==false); mp.lambda0 = 500e-9; end  % central wavelength of bandpass (meters)
if(isfield(mp,'fracBW')==false); mp.fracBW = 0.01; end   % fractional bandwidth of correction (Delta lambda / lambda)
if(isfield(mp,'Nsbp')==false); mp.Nsbp = 1; end  % number of sub-bandpasses across correction band 
if(isfield(mp,'Nwpsbp')==false); mp.Nwpsbp = 1; end % number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.


%% Need TO DOUBLE CHECK THAT THIS CELL DOESN'T OVERWRITE/OVERRULE SOMETHING ACCIDENTALLY
%% Coronagraphic Mask Properties:

%%--Pupil Masks
mp = falco_config_load_pupil_defaults(mp);

%%--Apodizer (Shaped Pupil) Properties (Plane P3). Also contains FPM and
%   Lyot stop specifications for that apodized coronagraph design. 
mp = falco_config_load_apodizer_defaults(mp);

%--FPM parameters  
% if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 8; end % sampling of FPM for full model, in pixels per lambda0/D
% if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 8; end % sampling of FPM for compact model, in pixels per lambda0/D

% %--FPM (F3) values (if not already defined in falco_config_load_apodizer_defaults)
% if(isfield(mp.F3,'Rin')==false); mp.F3.Rin = 1.5; end % inner hard-edge radius of the focal plane mask, in lambda0/D
% if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = Inf; end % outer hard-edge radius of the focal plane mask, in lambda0/D
% if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees

%%--Final Focal Plane (F4) Properties
mp = falco_config_load_F4_defaults(mp);


%%

%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        
        if(isfield(mp.P1,'D')==false); mp.P1.D = 4; end  %--meters, diameter of telescope (This is like HabEx A)
        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 250; end  
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 250; end 
        
        if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        if(isfield(mp.P1,'ODnorm')==false); mp.P1.ODnorm = 1; end % Outer diameter (fraction of Nbeam) 
        
        if(isfield(mp.P4,'IDnorm')==false); mp.P4.IDnorm = 0; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        if(isfield(mp.P4,'ODnorm')==false); mp.P4.ODnorm = 0.95; end % Outer diameter (fraction of Nbeam) 
        
        if(isfield(mp.P1,'num_strut')==false); mp.P1.num_strut = 0; end % Number of struts 
        if(isfield(mp.P1,'strut_angs')==false); mp.P1.strut_angs = []; end %Array of angles of the radial struts (deg)
        if(isfield(mp.P1,'strut_width')==false); mp.P1.strut_width = []; end  % Width of the struts (fraction of pupil diam.)
        
        if(isfield(mp.P4,'num_strut')==false); mp.P4.num_strut = 0; end % Number of struts 
        if(isfield(mp.P4,'strut_angs')==false); mp.P4.strut_angs = []; end %Array of angles of the radial struts (deg)
        if(isfield(mp.P4,'strut_width')==false); mp.P4.strut_width = []; end  % Width of the struts (fraction of pupil diam.)
                
    case{'WFIRST180718''WFIRST20180103','WFIRST_onaxis'}
        
        if(isfield(mp.P1,'D')==false); mp.P1.D = 2.3631; end  %--meters, diameter of telescope (used only for mas to lambda/D conversion)
        if(isfield(mp,'pup_strut_width')==false); mp.pup_strut_width = 3.22/100.; end  %--3.22% is the new value as of 2018-01-03
        if(isfield(mp,'LS_strut_width')==false); mp.LS_strut_width = 3.8/100.; end  

        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 324; end  %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 324; end 
        
        if(isfield(mp.P4,'IDnorm')==false); mp.P4.IDnorm = 0.50; end  % value for original pupil was 2.46/8 = 0.3075
        if(isfield(mp.P4,'ODnorm')==false); mp.P4.ODnorm = 0.80; end 
        
    case{'LUVOIRA5'}
        if(isfield(mp.P1,'D')==false); mp.P1.D = 15.2;  end %--meters, circumscribing diameter of telescope (used only for mas-to-lambda/D conversion)
        if(isfield(mp.P1,'Dfac')==false); mp.P1.Dfac = 15.2/13.7; end  %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        
        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 1000; end %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 500; end 
        
        if(isfield(mp,'LS_strut_width')==false); mp.LS_strut_width = 1.4/100.; end  

        
    case 'LUVOIR_B_offaxis'         % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        if(isfield(mp.P1,'D')==false); mp.P1.D = 7.989; end  %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.
        
        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 1000; end  %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering)

        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 500; end 
        
        if(isfield(mp.P4,'IDnorm')==false); mp.P4.IDnorm = 0; end 
        if(isfield(mp.P4,'ODnorm')==false); mp.P4.ODnorm = 0.80; end 
        
        
end

mp.P4.full.Nbeam = mp.P1.full.Nbeam;  % P4 must be the same as P1 for Vortex and HLC. 
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;  % P4 must be the same as P1 for Vortex and HLC. 
        

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



%%--FPM parameters
if(isfield(mp.F3,'Rin')==false); mp.F3.Rin = 1.5/2; end % inner hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = Inf; end % outer hard-edge radius of the focal plane mask, in lambda0/D

if(isfield(mp,'FPMampFac')==false); mp.FPMampFac = 1.0; end %--amplitude transmission value of the spot (achromatic)

if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 50; end % sampling of FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 30; end % sampling of FPM for compact model, in pixels per lambda0/D

if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees

if(isfield(mp,'fl')==false); mp.fl = 1;  end % meters. Arbitrary value chose for this design configuration. Keep as 1. Used for all focal lengths to keep magnification at each pupil at 1x. 


if(isfield(mp,'FPMmaterial')==false); mp.FPMmaterial = 'FS';  end % Use fused silica by default
if(isfield(mp.F3,'t')==false); mp.F3.t = mp.lambda0/2/(1.45-1); end % outer hard-edge radius of the focal plane mask, in lambda0/D



%% Controller weights and thresholds

%%--Tip/Tilt, Spatial, and Chromatic Weighting of the Control Jacobian
if(isfield(mp,'Ntt')==false); mp.Ntt = 1; end  %--Number of tip/tilt offsets, including none (so always set >=1). 1, 4, or 5
if(mp.Ntt == 1) %--mp.Ntt=1 means no tip/tilt offsets used
    mp.NlamForTT = 0;
end
if(isfield(mp,'NlamForTT')==false); mp.NlamForTT = 3; end  %--Number of wavelengths to compute tip/tilt at. 0,1, 2, 3, or inf (for all)
if(isfield(mp,'TToffset')==false); mp.TToffset = 1; end  %--tip/tilt offset in mas

%--Spatial pixel weighting of the Control Jacobian
if(isfield(mp,'WspatialDef')==false); mp.WspatialDef = []; end %[mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1]; end   %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--Chromatic weighting

%--Threshold level for culling actuators from the Jacobian
if(isfield(mp,'logGmin')==false); mp.logGmin = -7; end  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators


%% Deformable Mirror (DM) Parameters

if(isfield(mp,'dm_ind')==false); mp.dm_ind = [1 2]; end% vector of which DMs to use for control.
if(isfield(mp,'dm_weights')==false); mp.dm_weights = ones(9,1);  end % vector of relative weighting of DMs for EFC

%--DM1 parameters
if(isfield(mp.dm1,'Nact')==false); mp.dm1.Nact = 48; end  % number of actuators across DM1
if(isfield(mp.dm1,'VtoH')==false); mp.dm1.VtoH = 1*1e-9*ones(mp.dm1.Nact); end  % Gains: volts to meters in surface height;
if(isfield(mp.dm1,'xtilt')==false); mp.dm1.xtilt = 0; end 
if(isfield(mp.dm1,'ytilt')==false); mp.dm1.ytilt = 0; end 
if(isfield(mp.dm1,'zrot')==false); mp.dm1.zrot = 0; end  %--clocking angle (degrees)
if(isfield(mp.dm1,'xc')==false); mp.dm1.xc = (mp.dm1.Nact/2 - 1/2); end  % x-center of DM in actuator widths
if(isfield(mp.dm1,'yc')==false); mp.dm1.yc = (mp.dm1.Nact/2 - 1/2); end  % x-center of DM in actuator widths
if(isfield(mp.dm1,'edgeBuffer')==false); mp.dm1.edgeBuffer = 1; end  % Radius (in actuator spacings) outside of pupil to compute influence functions for.

%--DM2 parameters
if(isfield(mp.dm2,'Nact')==false); mp.dm2.Nact = 48; end  % number of actuators across DM1
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
% mp.dm9.Nrings = 25; % Number of radial hexagonal rings (including the center) for the FPM's phase "DM"
mp.dm9.VtoHavg = 1e-9; %--meters, gain of each "actuator" for the FPM phase.

%--DM9 parameters
mp.dm9.flag_hex_array = false; %--true->hex grid. false->square grid
mp.dm9.actres = 10; % number of "actuators" per lambda0/D in the FPM's focal plane. Only used for square grid
mp.dm9.Nact = ceil_even(2*mp.F3.Rin*mp.dm9.actres); % number of actuators across DM9 (if not in a hex grid)
mp.dm9.xtilt = 0;
mp.dm9.ytilt = 0;
mp.dm9.zrot = 0; %--clocking angle (degrees)
% mp.dm9.xc = mp.dm9.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% mp.dm9.yc = mp.dm9.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
mp.dm9.edgeBuffer = 0; %1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.

%--DM9 Actuator characteristics
mp.dm9.dx_inf0_act = 1/10; % units of actuator spacings, sampling of the influence function (for the Xinetics influence function, always 1/10 of the actuator pitch).
mp.dm9.inf0 = 1*fitsread('influence_dm5v2.fits');    %  -1*fitsread('inf64.3.fits');                              
%--Cropped influence function for FPM phase
mp.dm9.Ncrop_inf0 = 45;%20;%45;  %--Half the amount to extract from the center of the 91x91 pixel DM influence function. >=45 gives the uncropped influence functionNhalf = ceil(length(mp.dm9.inf0)/2);
mp.dm9.FPMbuffer = 0.2;%0.25; %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)

%%--Pupil Plane Properties
if(isfield(mp.P2,'D')==false);  mp.P2.D = (mp.dm1.Nact-2)*1e-3; end %46.3e-3; end % beam diameter at pupil closest to the DMs  (meters)
if(isfield(mp.P3,'D')==false);  mp.P3.D = mp.P2.D; end  % beam diameter at apodizer plane pupil (meters).
if(isfield(mp.P4,'D')==false);  mp.P4.D = mp.P2.D; end  % beam diameter at Lyot plane pupil (meters).



%% Apodizer (Shaped Pupil) Properties (Plane P3)


   



%% Lyot Stop Properties

%--Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle
mp.P4.full.Nbeam = mp.P1.full.Nbeam; %--Number of pixels across the re-imaged pupil at the Lyot plane (independent of beam centering)      
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)


switch mp.whichPupil
	case{'Simple'}
        %--Make or read in Lyot stop (LS) for the 'full' model
        if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0; end % Inner radius of Lyot stop
        if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.99; end % Outer radius of Lyot stop 
        if(isfield(mp.P4,'LS_num_strut')==false); mp.LS_num_strut = 4; end  % Number of struts in Lyot stop 
        if(isfield(mp.P4,'LS_strut_angs')==false); mp.LS_strut_angs = [0 90 180 270]; end %Angles of the struts 
        if(isfield(mp.P4,'LS_strut_width')==false); mp.LS_strut_width = 0.01; end % Size of Lyot stop spiders 
        
    case{'WFIRST_onaxis','WFIRST20180103'} % WFIRST is case specific 

    case{'LUVOIRA5'} % LUVOIR

end

%% Final Focal Plane (F4) Properties


%--Specs for Correction (Corr) region and the Scoring (Score) region.
if(isfield(mp.F4.corr,'Rin')==false); mp.F4.corr.Rin  = mp.F3.Rin; end  %--lambda0/D, inner radius of correction region
if(isfield(mp.F4.score,'Rin')==false); mp.F4.score.Rin = mp.F4.corr.Rin; end  %--Needs to be >= that of Correction mask
if(isfield(mp.F4.corr,'Rout')==false); mp.F4.corr.Rout  = floor(mp.dm1.Nact/2/(1+mp.fracBW/2)); end  %--lambda0/D, outer radius of correction region
if(isfield(mp.F4.score,'Rout')==false); mp.F4.score.Rout = mp.F4.corr.Rout; end  %--Needs to be <= that of Correction mask

if(isfield(mp.F4.corr,'ang')==false); mp.F4.corr.ang  = 180; end  %--degrees per side
if(isfield(mp.F4.score,'ang')==false); mp.F4.score.ang = 180; end  %--degrees per side
if(isfield(mp.F4,'sides')==false); mp.F4.sides = 'both'; end  %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%%--Final Focal Plane (F4) Properties
if(isfield(mp.F4.compact,'res')==false); mp.F4.res = 4; end  %--Pixels per lambda_c/D
if(isfield(mp.F4.full,'res')==false); mp.F4.full.res = 6; end  %--Pixels per lambda_c/D
if(isfield(mp.F4,'FOV')==false); mp.F4.FOV = 1 + mp.F4.corr.Rout; end  % minimum desired field of view (along both axes) in lambda0/D


end %--END OF FUNCTION



