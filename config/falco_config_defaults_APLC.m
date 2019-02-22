% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Function to set variables to default values if they are not already defined.
%
% Modified on 2018-05-19 by G. Ruane to be for the APLC.
% Modified on 2018-03-27 by A.J. Riggs to be for the SPLC.
% Modified on 2018-03-22 by A.J. Riggs to have default values that can be
%   overwritten if the variable is already defined.
% Created on 2017-10-31 by A.J. Riggs.
% Modified on 2018-01-08 by A.J. Riggs to be vortex specific and move key
%   parameters outside this function.
% Created on 2017-10-31 by A.J. Riggs.
%




function mp = falco_config_defaults_APLC(mp)

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
if(isfield(mp,'coro')==false); mp.coro = 'SPLC';  end %--Tested Options: 'HLC','SPLC','Vortex'
if(isfield(mp,'flagApod')==false); mp.flagApod = true;  end %--Can be any mask, even just an iris 
if(isfield(mp,'whichPupil')==false); mp.whichPupil = 'LUVOIRA0'; end %'Simple';

%%--General:
if(isfield(mp,'centering')==false); mp.centering = 'pixel'; end %--Centering on the arrays at each plane: pixel or interpixel


%%--Pupil Plane and DM Plane Properties
if(isfield(mp,'d_P2_dm1')==false); mp.d_P2_dm1 = 0; end  % distance (along +z axis) from P2 pupil to DM1 (meters)
if(isfield(mp,'d_dm1_dm2')==false); mp.d_dm1_dm2 = 3; end  % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
if(isfield(mp,'lambda0')==false); mp.lambda0 = 500e-9; end  % central wavelength of bandpass (meters)
if(isfield(mp,'fracBW')==false); mp.fracBW = 0.01; end   % fractional bandwidth of correction (Delta lambda / lambda)
if(isfield(mp,'Nsbp')==false); mp.Nsbp = 1; end  % number of sub-bandpasses across correction band 
if(isfield(mp,'Nwpsbp')==false); mp.Nwpsbp = 1; end % number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.

%% Coronagraphic Mask Properties:

%%--Pupil Masks
mp = falco_config_load_pupil_defaults(mp);

%% Controller Settings

%--Voltage range restrictions
if(isfield(mp.dm1,'maxAbsV')==false); mp.dm1.maxAbsV = 250./2.; end 
if(isfield(mp.dm2,'maxAbsV')==false); mp.dm2.maxAbsV = 250./2.; end 

%% Make new PSD maps for full model (NOT USED FOR DESIGN-ONLY CODE)
if(isfield(mp,'flagNewPSD')==false); mp.flagNewPSD = false; end  %--to make new PSD maps

%% Whether to include planet in the images
if(isfield(mp,'planetFlag')==false); mp.planetFlag = false; end 

%% Throughput calculations
if(isfield(mp,'thput_radius')==false); mp.thput_radius = 0.7; end  %--lambda_c/D, MAX radius to score core throughput
if(isfield(mp,'thput_eval_x')==false); mp.thput_eval_x = 6; end  %--lambda_c/D, x-axis value at which to evaluate throughput
if(isfield(mp,'thput_eval_y')==false); mp.thput_eval_y = 0; end  %--lambda_c/D, x-axis value at which to evaluate throughput

%% Controller weights and thresholds

%--Spatial pixel weighting of the Control Jacobian
if(isfield(mp,'WspatialDef')==false); mp.WspatialDef = []; end %[mp.F4.corr.Rin, mp.F4.corr.Rin+2, 1]; end   %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--Chromatic weighting

%--Threshold level for culling actuators from the Jacobian
if(isfield(mp,'logGmin')==false); mp.logGmin = -6; end  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators


%% Deformable Mirror (DM) Parameters

mp = falco_config_load_dm_defaults(mp);

%--Aperture stops at DMs
if(isfield(mp,'flagDM1stop')==false); mp.flagDM1stop = false; end  %--logical flag whether to include the stop at DM1 or not
if(isfield(mp,'flagDM2stop')==false); mp.flagDM2stop = false; end  %--logical flag whether to include the stop at DM2 or not
if(isfield(mp.dm1,'Dstop')==false); mp.dm1.Dstop = mp.dm1.Nact*mp.dm1.dm_spacing; end  %--diameter of circular stop at DM1 and centered on the beam
if(isfield(mp.dm2,'Dstop')==false); mp.dm2.Dstop = mp.dm2.Nact*mp.dm2.dm_spacing; end  %--diameter of circular stop at DM2 and centered on the beam

%%--Pupil Plane Properties
if(isfield(mp.P2,'D')==false);  mp.P2.D = (mp.dm1.Nact-2)*mp.dm1.dm_spacing; end %46.3e-3; end % beam diameter at pupil closest to the DMs  (meters)
if(isfield(mp.P3,'D')==false);  mp.P3.D = mp.P2.D; end  % beam diameter at apodizer plane pupil (meters).
if(isfield(mp.P4,'D')==false);  mp.P4.D = mp.P2.D; end  % beam diameter at Lyot plane pupil (meters).


%% Apodizer (Shaped Pupil Properties (Plane P3)

mp = falco_config_load_apodizer_defaults(mp);

%% Coronagraphic Mask Properties:

%%--FPM parameters
if(isfield(mp.F3,'Rin')==false); mp.F3.Rin = 3; end % inner hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = Inf; end % outer hard-edge radius of the focal plane mask, in lambda0/D

if(isfield(mp,'FPMampFac')==false); mp.FPMampFac = 0; end %--amplitude transmission value of the spot (achromatic)

if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 50; end % sampling of FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 30; end % sampling of FPM for compact model, in pixels per lambda0/D

if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees

if(isfield(mp,'fl')==false); mp.fl = 1;  end % meters. Arbitrary value chose for this design configuration. Keep as 1. Used for all focal lengths to keep magnification at each pupil at 1x. 


%% Lyot Stop Properties

% Now handled in falco_config_load_apodizer_defaults

%% Final Focal Plane (F4) Properties

%%--Final Focal Plane (F4) Properties
mp = falco_config_load_F4_defaults(mp);


end %--END OF FUNCTION



