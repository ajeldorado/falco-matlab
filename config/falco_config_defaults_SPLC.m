% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Function to set variables to default values if they are not already defined.
%
% Modified on 2018-03-27 by A.J. Riggs to be for the SPLC.
% Modified on 2018-03-22 by A.J. Riggs to have default values that can be
%   overwritten if the variable is already defined.
% Created on 2017-10-31 by A.J. Riggs.
% Modified on 2018-01-08 by A.J. Riggs to be vortex specific and move key
%   parameters outside this function.
% Created on 2017-10-31 by A.J. Riggs.
%




function mp = falco_config_defaults_SPLC(mp)

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

%%--Apodizer (Shaped Pupil) Properties (Plane P3). Also contains FPM and
%   Lyot stop specifications for that apodized coronagraph design. 
mp = falco_config_load_apodizer_defaults(mp);

%--FPM parameters  
% if(isfield(mp.F3.full.in,'res')==false);    mp.F3.full.in.res = 30; end % sampling of inner FPM for full model, in pixels per lambda0/D
% if(isfield(mp.F3.full.out,'res')==false);   mp.F3.full.out.res = 4; end % sampling of outer FPM for full model, in pixels per lambda0/D
% if(isfield(mp.F3.compact.in,'res')==false);  mp.F3.compact.in.res = 30; end % sampling of inner FPM for compact model, in pixels per lambda0/D
% if(isfield(mp.F3.compact.out,'res')==false); mp.F3.compact.out.res = 4; end % sampling of outer FPM for compact model, in pixels per lambda0/D
if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 10; end % sampling of FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 6; end % sampling of FPM for compact model, in pixels per lambda0/D

%--FPM shape (if not already defined in falco_config_load_apodizer_defaults)
if(isfield(mp.F3,'Rin')==false); mp.F3.Rin = 2.8; end % inner hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = 20; end % outer hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees

%%--Final Focal Plane (F4) Properties
mp = falco_config_load_F4_defaults(mp);


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

%% Coronagraphic Mask Properties:

if(isfield(mp,'fl')==false); mp.fl = 1;  end % meters. Arbitrary value chose for this design configuration. Keep as 1. Used for all focal lengths to keep magnification at each pupil at 1x. 

%%--Pupil Plane Properties
if(isfield(mp.P2,'D')==false);  mp.P2.D = 46.3e-3; end % beam diameter at pupil closest to the DMs  (meters)
if(isfield(mp.P3,'D')==false);  mp.P3.D = mp.P2.D; end  % beam diameter at apodizer plane pupil (meters).
if(isfield(mp.P4,'D')==false);  mp.P4.D = mp.P2.D; end  % beam diameter at Lyot plane pupil (meters).


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
if(isfield(mp,'logGmin')==false); mp.logGmin = -6; end  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators


%% Deformable Mirror (DM) Parameters

if(isfield(mp,'dm_ind')==false); mp.dm_ind = [1 2]; end% vector of which DMs to use for control. dm9 is the FPM phase
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

% %% Apodizer (Shaped Pupil Properties (Plane P3)
% if(isfield(mp,'SPname')==false); mp.SPname = 'luvoirA5bw10'; end
% 
% switch mp.whichPupil
%     case{'WFIRST_onaxis','WFIRST20180103'}
%         switch mp.SPname
%             case '20170714' %-- 18% Bandwidth
%                 %--Dummy values
%                 mp.rEdgesLeft = 1;
%                 mp.rEdgesRight = 1;
%                 
%             case '32WA194' %-- 18% Bandwidth, concentric rings
%                 mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot100_3300Dpup9850_32WA194_45LS79_8FPres4_BW18N9_c90_Nring9_D20mm_10um_left_WS.dat');
%                 mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot100_3300Dpup9850_32WA194_45LS79_8FPres4_BW18N9_c90_Nring9_D20mm_10um_right_WS.dat');  
%             case '31WA220' %-- 10% Bandwidth,  concentric rings
%                 mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_3300Dpup9850_31WA220_43LS82_8FPres4_BW10N9_c90_Nring10_D20mm_10um_left_WS.dat');
%                 mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_3300Dpup9850_31WA220_43LS82_8FPres4_BW10N9_c90_Nring10_D20mm_10um_right_WS.dat');  
%         end
% 
%     case{'LUVOIRA0'}
%         if(strcmpi(mp.SPname,'luvoirA0bw10')) %-- 10% Bandwidth
%             mp.rEdgesLeft  = load('out_RSPLC1D_A_maxTrPH_2848Dpup9960_33WA228_38LS83_10FPres5_BW10N11_c100_Nring10_left.dat');
%             mp.rEdgesRight = load('out_RSPLC1D_A_maxTrPH_2848Dpup9960_33WA228_38LS83_10FPres5_BW10N11_c100_Nring10_right.dat');  
%         end  
% 
%         
%     case{'LUVOIRA5'}
%         if(strcmpi(mp.SPname,'luvoirA5bw10')) %-- 10% Bandwidth
%             mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_1123Dpup9900_34WA264_24LS78_10FPres4_BW10N6_c100_Nring14_D20mm_10um_left_WS.dat');
%             mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_1123Dpup9900_34WA264_24LS78_10FPres4_BW10N6_c100_Nring14_D20mm_10um_right_WS.dat');  
%         end
%         
%     case{'LUVOIR_B_offaxis'}
%         if(strcmpi(mp.SPname,'luvoirBbw20')) %-- 10% Bandwidth
%             mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot400_0Dpup9900_25WA350_6LS74_10FPres4_BW20N9_c103_Nring19_D20mm_10um_left_WS.dat');
%             mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot400_0Dpup9900_25WA350_6LS74_10FPres4_BW20N9_c103_Nring19_D20mm_10um_right_WS.dat');  
%         end
% end
% 
% %--Scale the ring diameters from the inscribed diameter values to the
% %  circumscribed diameter values.
% mp.rEdgesLeft = mp.rEdgesLeft/mp.P1.Dfac;
% mp.rEdgesRight = mp.rEdgesRight/mp.P1.Dfac;
% 
% 
% if(strcmpi(mp.SPname,'20170714'))
%     mp.F3.Rin = 2.6; % inner hard-edge radius of the focal plane mask, in lambda0/D
%     mp.F3.Rout = 9.0; % outer hard-edge radius of the focal plane mask, in lambda0/D
%     mp.F3.ang = 65; % opening angle of horizontal bowtie (degrees)
% elseif(strcmpi(mp.SPname,'32WA194'))
%     mp.F3.Rin = 3.2; % inner hard-edge radius of the focal plane mask, in lambda0/D
%     mp.F3.Rout = 19.4; % outer hard-edge radius of the focal plane mask, in lambda0/D
% elseif(strcmpi(mp.SPname,'31WA220'))
%     mp.F3.Rin = 3.1; % inner hard-edge radius of the focal plane mask, in lambda0/D
%     mp.F3.Rout = 22.0; % outer hard-edge radius of the focal plane mask, in lambda0/D
% elseif(strcmpi(mp.SPname,'luvoirA0bw10'))
%     mp.F3.Rin = 3.3*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
%     mp.F3.Rout = 22.8*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
% elseif(strcmpi(mp.SPname,'luvoirA5bw10'))
% %     if( exist('mp.P1.Dfac','var')==false ); mp.P1.Dfac = 15.2/13.7; end %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
%     mp.F3.Rin = 3.367*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
%     mp.F3.Rout = 26.4*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
% elseif(strcmpi(mp.SPname,'luvoirBbw20'))
% %     if( exist('mp.P1.Dfac','var')==false ); mp.P1.Dfac = 15.2/13.7; end %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
%     mp.F3.Rin = 2.5*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
%     mp.F3.Rout = 35.0*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
% end

%% Lyot Stop Properties
% switch mp.whichPupil
% 	case{'Simple'}
%         %--Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle
%         mp.P4.full.Nbeam = mp.P1.full.Nbeam; %--Number of pixels across the re-imaged pupil at the Lyot plane (independent of beam centering)
%         mp.P4.full.Nbeam
%         mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
%         %--Make or read in Lyot stop (LS) for the 'full' model
%         if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.3; end % Inner radius of Lyot stop
%         if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.85; end % Outer radius of Lyot stop 
%         mp.LS_num_strut = 4; % Number of struts in Lyot stop 
%         mp.LS_strut_angs = [0 90 180 270];%Angles of the struts 
%         mp.LS_strut_width = 0.01;% Size of Lyot stop spiders 
%         
%     case{'WFIRST_onaxis','WFIRST20180103'} % WFIRST is case specific 
% 
%         %--Make or read in Lyot stop (LS) for the 'full' model
%         if(strcmpi(mp.SPname,'20170714'))
%             if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.25; mp.P4.ODnorm = 0.82; mp.P4.ang = 115;  end
%             mp.LS_strut_width = 0; %--Dummy value
%         elseif(strcmpi(mp.SPname,'32WA194'))
%             if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.45; mp.P4.ODnorm = 0.79; end
%         elseif(strcmpi(mp.SPname,'31WA220'))
%             if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.43; mp.P4.ODnorm = 0.82; end
%         end
% 
%     case{'LUVOIRA0'} % LUVOIR with big central obscuration
%         if(strcmpi(mp.SPname,'luvoirA0bw10'))
%             if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.38/mp.P1.Dfac; end
%             if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.83/mp.P1.Dfac; end
%         end
%         
%     case{'LUVOIRA5'} % LUVOIR
%         if(strcmpi(mp.SPname,'luvoirA5bw10'))
%             if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.24/mp.P1.Dfac; end
%             if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.78/mp.P1.Dfac; end
%         end
%         
%         
%     case{'LUVOIR_B_offaxis'} % LUVOIR
%         if(strcmpi(mp.SPname,'luvoirBbw20'))
%             if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.06/mp.P1.Dfac; end
%             if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.74/mp.P1.Dfac; end
%         end
% end

%% Last Two Focal Plane (F3 and F4) Properties

%%--FPM parameters
if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 10; end % sampling of FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 5; end % sampling of FPM for compact model, in pixels per lambda0/D

if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees

switch mp.coro
    case {'SPLC'} %--Occulting spot and diaphragm FPM (opaque spot and opaque outer diaphragm)
        if(isfield(mp,'FPMampFac')==false); mp.FPMampFac = 0; end %--amplitude transmission value of the spot (achromatic)
    case {'SPHLC'} %--Complex occulting spot and opaque diaphragm FPM
%         mp.FPMampFac = 0.025; %--amplitude transmission value of the spot (achromatic)
end



%--Specs for Correction (Corr) region and the Scoring (Score) region.
% if(isfield(mp.F4.corr,'Rin')==false); mp.F4.corr.Rin  = mp.F3.Rin; end  %--lambda0/Dcircumscribed, inner radius of correction region
if(isfield(mp.F4.corr,'Rin')==false)
    if(isfield(mp.F3,'RinA')) 
        mp.F4.corr.Rin  = mp.F3.RinA;
    else
        mp.F4.corr.Rin  = mp.F3.Rin;
    end
end  %--lambda0/D, inner radius of correction region
if(isfield(mp.F4.score,'Rin')==false); mp.F4.score.Rin = mp.F4.corr.Rin; end  %--Needs to be >= that of Correction mask
if(isfield(mp.F4.corr,'Rout')==false); mp.F4.corr.Rout  = min( [floor(mp.dm1.Nact/2*(1-mp.fracBW/2)), mp.F3.Rout ] ); end %--lambda0/Dcircumscribed, outer radius of correction region
if(isfield(mp.F4.score,'Rout')==false); mp.F4.score.Rout = mp.F4.corr.Rout; end  %--Needs to be <= that of Correction mask
if(isfield(mp.F4.corr,'ang')==false); mp.F4.corr.ang  = 180; end  %--degrees per side
if(isfield(mp.F4.score,'ang')==false); mp.F4.score.ang = 180; end  %--degrees per side

if(isfield(mp.F4,'sides')==false); mp.F4.sides = 'both'; end  %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%%--Final Focal Plane (F4) Properties
if(isfield(mp.F4.compact,'res')==false); mp.F4.res = 2.5; end  %--Pixels per lambda_c/D
if(isfield(mp.F4.full,'res')==false); mp.F4.full.res = 4; end  %--Pixels per lambda_c/D
if(isfield(mp.F4,'FOV')==false); mp.F4.FOV = 1 + mp.F4.corr.Rout; end  % minimum desired field of view (along both axes) in lambda0/D


end %--END OF FUNCTION



