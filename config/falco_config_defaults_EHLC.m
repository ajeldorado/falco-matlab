% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Function to set variables to default values if they are not already defined.
%
% Modified on 2018-03-22 by A.J. Riggs to have default values that can be
%   overwritten if the variable is already defined.
% Created on 2017-10-31 by A.J. Riggs.
% Modified on 2018-01-08 by A.J. Riggs to be vortex specific and move key
%   parameters outside this function.
% Created on 2017-10-31 by A.J. Riggs.
%




function mp = falco_config_defaults_EHLC(mp)


%% Initialize some structures if they don't already exist
mp.ctrl.dummy = 1;
mp.est.dummy = 1;

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
if(isfield(mp,'coro')==false); mp.coro = 'EHLC';  end %--Tested Options: 'HLC','SPLC','SPHLC'
if(isfield(mp,'flagApod')==false); mp.flagApod = false;  end %--Can be any mask, even just an iris 
if(isfield(mp,'whichPupil')==false); mp.whichPupil = 'WFIRST180718'; end %'Simple';

%%--General:
if(isfield(mp,'centering')==false); mp.centering = 'pixel'; end %--Centering on the arrays at each plane: pixel or interpixel


%%--Pupil Plane and DM Plane Properties
if(isfield(mp,'d_P2_dm1')==false); mp.d_P2_dm1 = 0; end  % distance (along +z axis) from P2 pupil to DM1 (meters)
if(isfield(mp,'d_dm1_dm2')==false); mp.d_dm1_dm2 = 1; end  % distance between DM1 and DM2 (meters)

%%--Bandwidth and Wavelength Specs
if(isfield(mp,'lambda0')==false); mp.lambda0 = 575e-9; end  % central wavelength of bandpass (meters)
if(isfield(mp,'fracBW')==false); mp.fracBW = 0.01; end   % fractional bandwidth of correction (Delta lambda / lambda)
if(isfield(mp,'Nsbp')==false); mp.Nsbp = 1; end  % number of sub-bandpasses across correction band 
if(isfield(mp,'Nwpsbp')==false); mp.Nwpsbp = 1; end % number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Can be odd or even value.

%% Focal Lengths

if(isfield(mp,'fl')==false); mp.fl = 1;  end % meters. Arbitrary value chose for this design configuration. Keep as 1. Used for all focal lengths to keep magnification at each pupil at 1x. 

%% Input Pupil Mask Properties:
mp = falco_config_load_pupil_defaults(mp);


%% FPM (focal plane F3) 

mp.FPM.d0fac = 4; %--For HLC occulter: Reference plane's distance from the substrate surface [waves]

%--FPM resolution  
if(isfield(mp.F3.full,'res')==false); mp.F3.full.res = 20; end % sampling of FPM for full model, in pixels per lambda0/D
if(isfield(mp.F3.compact,'res')==false); mp.F3.compact.res = 20; end % sampling of FPM for compact model, in pixels per lambda0/D

%--FPM (F3) values (if not already defined in falco_config_load_apodizer_defaults)
if(isfield(mp.F3,'Rin')==false); mp.F3.Rin = 2.7; end % inner hard-edge radius of the focal plane mask, in lambda0/D
if(isfield(mp.F3,'RinMaxMetal')==false); mp.F3.RinMaxMetal = 6; end % max extent of actuatable metal part of FPM, in lambda0/D #NEW4EHLC
if(isfield(mp.F3,'RinMaxDiel')==false); mp.F3.RinMaxDiel = 10; end % max extent of actuatable dielectric part of FPM, in lambda0/D #NEW4EHLC
if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = 30; end % outer hard-edge radius of the focal plane mask, in lambda0/D

if(isfield(mp.F3,'ang')==false); mp.F3.ang = 180; end% angular opening on each side of the focal plane mask, in degrees


%% Controller Settings

%--Voltage range restrictions
if(isfield(mp.dm1,'maxAbsV')==false); mp.dm1.maxAbsV = 250./2.; end 
if(isfield(mp.dm2,'maxAbsV')==false); mp.dm2.maxAbsV = 250./2.; end 

if(isfield(mp.dm8,'dVmax')==false); mp.dm8.dVmax = 20; end 
 

%% Make new PSD maps for full model (NOT USED FOR DESIGN-ONLY CODE)
if(isfield(mp,'flagNewPSD')==false); mp.flagNewPSD = false; end  %--to make new PSD maps

%% Whether to include planet in the images
if(isfield(mp,'planetFlag')==false); mp.planetFlag = false; end 

%% Throughput calculations
if(isfield(mp,'thput_radius')==false); mp.thput_radius = 0.7; end  %--lambda_c/D, MAX radius to score core throughput
if(isfield(mp,'thput_eval_x')==false); mp.thput_eval_x = 6; end  %--lambda_c/D, x-axis value at which to evaluate throughput
if(isfield(mp,'thput_eval_y')==false); mp.thput_eval_y = 0; end  %--lambda_c/D, x-axis value at which to evaluate throughput

%% Controller weights and thresholds

%--DM9 regularization weight variation.
if(isfield(mp.ctrl,'dm9regfacVec')==false); mp.ctrl.dm9regfacVec = 10.^(-2:1:4); end

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
% % mp.dm9.Nrings = 25; % Number of radial hexagonal rings (including the center) for the FPM's phase "DM"
% mp.dm8.VtoHavg = 1e-9; %--meters, gain of each "actuator" for the FPM nickel.
% mp.dm9.VtoHavg = 1e-9; %--meters, gain of each "actuator" for the FPM dielectric.

%--DM9 parameters
% if(isfield(mp.dm9,'flag_hex_array')==false); mp.dm9.flag_hex_array = false; end %--true->hex grid. false->square grid
if(isfield(mp.dm9,'actres')==false);  mp.dm9.actres = 10; end % number of "actuators" per lambda0/D in the FPM's focal plane for DM9. Only used for square grid
if(isfield(mp.dm8,'actres')==false);  mp.dm8.actres = mp.dm9.actres; end % number of "actuators" per lambda0/D in the FPM's focal plane for DM8. Only used for square grid

% % % if(isfield(mp.dm8,'Vmin')==false);  mp.dm8.Vmin = 0; end % minimum thickness of FPM metal layer (nm)
% % % if(isfield(mp.dm8,'Vmax')==false);  mp.dm8.Vmax = 300; end % maximum thickness (from one actuator, not of the facesheet) of FPM metal layer (nm)
% % % if(isfield(mp.dm9,'Vmin')==false);  mp.dm9.Vmin = -mp.t_diel_bias_nm; end % minimum thickness of FPM dielectric layer (nm)
% % % if(isfield(mp.dm9,'Vmax')==false);  mp.dm9.Vmax = 400; end % maximum thickness (from one actuator, not of the facesheet) of FPM dielectric layer (nm)

if(isfield(mp.F3,'RinMaxDiel')==false);  mp.F3.RinMaxDiel = 10; end % max extent of actuatable dielectric part of FPM, in lambda0/D #NEW4EHLC
if(isfield(mp.F3,'RinMaxMetal')==false);  mp.F3.RinMaxMetal = 6; end % max extent of actuatable metal part of FPM, in lambda0/D #NEW4EHLC


mp.dm8.Nact = ceil_even(2*mp.F3.RinMaxMetal*mp.dm8.actres); % number of actuators across DM8 (assuming a square grid)
mp.dm9.Nact = ceil_even(2*mp.F3.RinMaxDiel*mp.dm9.actres); % number of actuators across DM9 (assuming a square grid)
% mp.dm9.Nact = ceil_even(2*mp.F3.Rin*mp.dm9.actres); % number of actuators across DM9 (assuming a square grid)
mp.dm9.xtilt = 0;
mp.dm9.ytilt = 0;
mp.dm9.zrot = 0; %--clocking angle (degrees)
% mp.dm9.xc = mp.dm9.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% mp.dm9.yc = mp.dm9.Nact/2 - 1/2; % x-center of DM in mm, in actuator widths
% mp.dm9.edgeBuffer = 1;%0; %1; % Radius (in actuator spacings) outside of pupil to compute influence functions for.


if(isfield(mp,'dm_weights')==false);   mp.dm_weights = ones(9,1); end % Jacobian weights

if(isfield(mp.dm8,'act_sens')==false);  mp.dm8.act_sens = 1;  end %--Change in oomph (E-field sensitivity) of DM8 actuators. Chosen empirically based on how much DM8 actuates during a control step.
if(isfield(mp.dm8,'stepFac')==false);  mp.dm8.stepFac = 5; end%--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.
if(isfield(mp.dm9,'act_sens')==false);  mp.dm9.act_sens = 1;  end %--Change in oomph (E-field sensitivity) of DM8 actuators. Chosen empirically based on how much DM8 actuates during a control step.
if(isfield(mp.dm9,'stepFac')==false);  mp.dm9.stepFac = 200; end%--Adjust the step size in the Jacobian, then divide back out. Used for helping counteract effect of discretization.



%--DM9 Actuator characteristics

%--Set the DM9 influence function
if(isfield(mp.dm9,'inf0name'))
    switch mp.dm9.inf0name
        case '3x3'
            mp.dm9.inf0 = 1/4.*[1, 2, 1; 2, 4, 2; 1, 2, 1];  % influence function
            mp.dm9.dx_inf0_act = 0.5;  % number of inter-actuator widths per pixel 
            %--FPM resolution (pixels per lambda0/D) in the compact and full models.
            mp.F3.compact.res = mp.dm9.actres/mp.dm9.dx_inf0_act;
            mp.F3.full.res = mp.dm9.actres/mp.dm9.dx_inf0_act;

        case 'Lanczos3'
            N = 91;
            xs = (-(N-1)/2:(N-1)/2)/N*10*0.906;
            %xsPlot = (-(N-1)/2:(N-1)/2);

            a = 3;
            Lx0 = a*sin(pi*xs).*sin(pi*xs/a)./(pi*xs).^2;
            Lx0(xs==0) = 1;
            Lx = Lx0;
            Lx(xs>=a) = 0;
            Lx(xs<=-a) = 0;
            Lxy = Lx.'*Lx; %--The 2-D Lanczos kernel
            Nhalf = ceil(N/2);
            Nrad = 30;
            inds = Nhalf-Nrad:Nhalf+Nrad;
            LxyCrop = Lxy(inds,inds);
            
            mp.dm9.inf0 = LxyCrop; % influence function
            mp.dm9.dx_inf0_act = 1/10;  % number of inter-actuator widths per pixel 
%             figure(11); imagesc(xsPlot,xsPlot,Lxy); axis xy equal tight; colorbar; set(gca,'Fontsize',20); title('Lanczos3 Influence Function','Fontsize',20); drawnow;
%             figure(12); imagesc(mp.dm9.inf0); axis xy equal tight; colorbar; set(gca,'Fontsize',20); drawnow; %title('FALCO','Fontsize',20);

        case 'Xinetics'
            mp.dm9.inf0 = 1*fitsread('influence_dm5v2.fits');
            mp.dm9.dx_inf0_act = 1/10;  % number of inter-actuator widths per pixel 
    end
    
    
else

    if(isfield(mp.dm9,'dx_inf0_act')==false); mp.dm9.dx_inf0_act = 1/10; end % units of actuator spacings, sampling of the influence function (for the Xinetics influence function, always 1/10 of the actuator pitch).
    if(isfield(mp.dm9,'inf0')==false); mp.dm9.inf0 = 1*fitsread('influence_dm5v2.fits');  end
end

%--Cropped influence function for FPM phase
%mp.dm9.Ncrop_inf0 = 45;%20;%45;  %--Half the amount to extract from the center of the 91x91 pixel DM influence function. >=45 gives the uncropped influence functionNhalf = ceil(length(mp.dm9.inf0)/2);
if(isfield(mp.dm9,'FPMbuffer')==false); mp.dm9.FPMbuffer = 0;  end %0.2; end %--Zero out DM9 actuators too close to the outer edge (within mp.dm9.FPMbuffer lambda0/D of edge)

%%--Pupil Plane Properties
if(isfield(mp.P2,'D')==false);  mp.P2.D = (mp.dm1.Nact-2)*1e-3; end %46.3e-3; end % beam diameter at pupil closest to the DMs  (meters)
if(isfield(mp.P3,'D')==false);  mp.P3.D = mp.P2.D; end  % beam diameter at apodizer plane pupil (meters).
if(isfield(mp.P4,'D')==false);  mp.P4.D = mp.P2.D; end  % beam diameter at Lyot plane pupil (meters).



%% Apodizer (Shaped Pupil) Properties (Plane P3)

% % %%--Apodizer (Shaped Pupil) Properties (Plane P3). Also contains FPM and
% % %   Lyot stop specifications for that apodized coronagraph design. 
% % mp = falco_config_load_apodizer_defaults(mp);

   

%% Lyot Stop Properties

%--Lyot stop pupil size
if(isfield(mp.P4.full,'Nbeam')==false); mp.P4.full.Nbeam = 200;  end
if(isfield(mp.P4.compact,'Nbeam')==false); mp.P4.compact.Nbeam = 200;  end


%% Final Focal Plane (F4) Properties

%--Resolution and FOV
if(isfield(mp.F4,'res')==false); mp.F4.res = 2.5; end  %--Pixels per lambda_c/D
if(isfield(mp.F4,'FOV')==false); mp.F4.FOV = 1.0 + mp.F4.corr.Rout; end  % minimum desired field of view (along both axes) in lambda0/D


%--Specs for Correction (Corr) region and the Scoring (Score) region.
if(isfield(mp.F4.corr,'Rin')==false);  mp.F4.corr.Rin  = mp.F3.Rin;  end  %--lambda0/D, inner radius of correction region
if(isfield(mp.F4.score,'Rin')==false); mp.F4.score.Rin = mp.F4.corr.Rin; end  %--Needs to be >= that of Correction mask
if(isfield(mp.F4.corr,'Rout')==false); mp.F4.corr.Rout  = min( [floor(mp.dm1.Nact/2*(1-mp.fracBW/2)), mp.F3.Rout ] ); end %--lambda0/Dcircumscribed, outer radius of correction region
if(isfield(mp.F4.score,'Rout')==false); mp.F4.score.Rout = mp.F4.corr.Rout; end  %--Needs to be <= that of Correction mask
if(isfield(mp.F4.corr,'ang')==false); mp.F4.corr.ang  = 180; end  %--degrees per side
if(isfield(mp.F4.score,'ang')==false); mp.F4.score.ang = 180; end  %--degrees per side
if(isfield(mp.F4,'sides')==false); mp.F4.sides = 'both'; end  %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


end %--END OF FUNCTION



