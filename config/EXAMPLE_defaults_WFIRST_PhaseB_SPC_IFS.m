
% %--Initialize some structures if they don't already exist

%% Misc

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 34;

%--Special Computational Settings
mp.flagParfor = true;
mp.useGPU = false;
mp.flagPlot = false;

%--General
mp.centering = 'pixel';

%--Whether to include planet in the images
mp.planetFlag = false;

%--Method of computing core throughput:
% - 'HMI' for energy within half-max isophote divided by energy at telescope pupil
% - 'EE' for encircled energy within a radius (mp.thput_radius) divided by energy at telescope pupil
mp.thput_metric = 'HMI'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 6; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%% Bandwidth and Wavelength Specs

mp.lambda0 = 730e-9;   %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.15;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 7;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp' for pairwise probing with batch process estimation
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
% - 'pwp-iekf' for pairwise probing with iterated extended Kalman filter  [NOT AVAILABLE YET]
mp.estimator = 'perfect';

%--New variables for pairwise probing estimation:
mp.est.probe.Npairs = 3;%2;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;%20;    % Max x/y extent of probed region [actuators].
mp.est.probe.offsetX = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.offsetY = 14;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%--New variables for pairwise probing with a Kalman filter
%  mp.est.ItrStartKF =  %Which correction iteration to start recursive estimate
%  mp.est.tExp =
%  mp.est.num_im =
%  mp.readNoiseStd =
%  mp.peakCountsPerPixPerSec =
%  mp.est.Qcoef =
%  mp.est.Rcoef =

%% Wavefront Control: General

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*[1]; %1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is not used; it is a placeholder for correct indexing of other modes.)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = 2:6; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [3., 4.;...
                4., 8.;
                8., 9.];  

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
% % mp.ctrl.dm9regfacVec = 1;        %--Additional regularization factor applied to DM9
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1.;
mp.dm2.weight = 1.;

%--Voltage range restrictions
mp.dm1.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.dm2.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.maxAbsdV = 1000;     %--Max +/- delta voltage step for each actuator for DMs 1 and 2 [volts] %--NOT ENFORCED YET

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 5; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use

% % % PLANNED SEARCH EFC DEFAULTS     
% mp.dm_ind = [1 2 ]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
% mp.ctrl.dmfacVec = 1;
% %--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
%     % Column 1: # of iterations, 
%     % Column 2: log10(regularization), 
%     % Column 3: which DMs to use (12, 128, 129, or 1289) for control
%     % Column 4: flag (0 = false, 1 = true), whether to re-linearize
%     %   at that iteration.
%     % Column 5: flag (0 = false, 1 = true), whether to perform an
%     %   EFC parameter grid search to find the set giving the best
%     %   contrast .
%     % The imaginary part of the log10(regularization) in column 2 is
%     %  replaced for that iteration with the optimal log10(regularization)
%     % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time
% 
% mp.ctrl.sched_mat = [...
%     repmat([1,1j,12,1,1],[4,1]);...
%     repmat([1,1j-1,12,1,1],[25,1]);...
%     repmat([1,1j,12,1,1],[1,1]);...
%     ];
% [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


%% Deformable Mirrors: Influence Functions
%--Influence Function Options:
% - 'influence_dm5v2.fits' for one type of Xinetics DM
% - 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM
% - 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM

mp.dm1.inf_fn = 'influence_dm5v2.fits';
mp.dm2.inf_fn = 'influence_dm5v2.fits';

mp.dm1.dm_spacing = 0.9906e-3;%1e-3; %--User defined actuator pitch
mp.dm2.dm_spacing = 0.9906e-3;%1e-3; %--User defined actuator pitch

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Deformable Mirrors: Optical Layout Parameters

%--DM1 parameters
mp.dm1.Nact = 48;               % # of actuators across DM array
mp.dm1.VtoH = 1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 5.83;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 48;               % # of actuators across DM array
mp.dm2.VtoH = 1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 5.55;%8;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;              % clocking of DM surface [degrees]
mp.dm2.xc = (48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = 100e-3;  %--Diameter of iris [meters]
mp.flagDM2stop = true;  %--Whether to apply an iris or not
mp.dm2.Dstop = 50e-3;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
mp.d_dm1_dm2 = 1.000;   % distance between DM1 and DM2 [meters]


%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'Fourier';  %--Which optical layout to use
mp.coro = 'SPLC';
mp.flagApod = true;    %--Whether to use an apodizer or not

%--Final Focal Plane Properties
mp.Fend.res = 6;%(730/660)*2.; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 11.; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 2.6;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 9;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 65;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 3;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 9;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 65;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'both'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1.; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 46.2987e-3;
mp.P3.D = 46.2987e-3;
mp.P4.D = 46.2987e-3;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 386;
mp.P2.compact.Nbeam = 386;
mp.P3.compact.Nbeam = 386;
mp.P4.compact.Nbeam = 60;

%--Shaped Pupil Mask: Load and downsample.
mp.SPname = 'SPC-20190130';
SP0 = fitsread('SPM_SPC-20190130.fits');
SP0 = SP0(2:end,2:end);
% figure(1); imagesc(SP0); axis xy equal tight; colormap jet; colorbar;
% figure(11); imagesc(SP0-fliplr(SP0)); axis xy equal tight; colormap jet; colorbar;
dx0 = 1/1000;
dx1 = 1/mp.P3.compact.Nbeam;
N0 = size(SP0,1);
switch lower(mp.centering)
    case{'pixel'}
        N1 = ceil_odd(N0*dx0/dx1);
    case{'interpixel'}
        N1 = ceil_even(N0*dx0/dx1);
end
x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
[X0,Y0] = meshgrid(x0);
R0 = sqrt(X0.^2+Y0.^2);
Window = 0*R0;
Window(R0<=dx1) = 1; Window = Window/sum(sum(Window));
% figure(10); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
SP0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(SP0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
SP0 = circshift(SP0,[1 1]); %--Undo a centering shift
x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
[X1,Y1] = meshgrid(x1);
SP1 = interp2(X0,Y0,SP0,X1,Y1,'cubic',0); %--Downsample by interpolation

switch lower(mp.centering)
    case{'pixel'}
        mp.P3.compact.mask = zeros(N1+1,N1+1);
        mp.P3.compact.mask(2:end,2:end) = SP1;
    otherwise
        mp.P3.compact.mask = SP1;
end
figure(2); imagesc(SP0); axis xy equal tight; colormap jet; colorbar;
figure(3); imagesc(SP1); axis xy equal tight; colormap jet; colorbar;
% figure(12); imagesc(SP0-fliplr(SP0)); axis xy equal tight; colormap jet; colorbar;
% figure(13); imagesc(SP1-fliplr(SP1)); axis xy equal tight; colormap jet; colorbar;

%%
%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;

%--FPM resolution
mp.F3.compact.res = 6;    % sampling of FPM for compact model [pixels per lambda0/D]

%--Load and downsample the FPM. To get good grayscale edges, convolve with the correct window before downsampling. 
FPM0 = fitsread('FPM_res100_SPC-20190130.fits'); %--Resolution of 100 pixels per lambda0/D
FPM0 = padOrCropOdd(FPM0,1821);
% figure(1); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% figure(11); imagesc(FPM0-rot90(FPM0,2)); axis xy equal tight; colormap jet; colorbar;
dx0 = 1/100;
dx1 = 1/mp.F3.compact.res;
N0 = size(FPM0,1);
switch lower(mp.centering)
    case{'pixel'}
        N1 = ceil_odd(N0*dx0/dx1);
    case{'interpixel'}
        N1 = ceil_even(N0*dx0/dx1);
end
x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
[X0,Y0] = meshgrid(x0);
R0 = sqrt(X0.^2+Y0.^2);
Window = 0*R0;
Window(R0<=dx1/2) = 1; Window = Window/sum(sum(Window));
% figure(10); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
FPM0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(FPM0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
FPM0 = circshift(FPM0,[1 1]); %--Undo a centering shift
x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
[X1,Y1] = meshgrid(x1);
FPM1 = interp2(X0,Y0,FPM0,X1,Y1,'cubic',0); %--Downsample by interpolation
switch lower(mp.centering)
    case{'pixel'}
        mp.F3.compact.mask.amp = zeros(N1+1,N1+1);
        mp.F3.compact.mask.amp(2:end,2:end) = FPM1;
    otherwise
        mp.F3.compact.mask.amp = FPM1;
end
% figure(2); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% figure(3); imagesc(FPM1); axis xy equal tight; colormap jet; colorbar;
% figure(12); imagesc(FPM0-rot90(FPM0,2)); axis xy equal tight; colormap jet; colorbar;
% figure(13); imagesc(FPM1-rot90(FPM1,2)); axis xy equal tight; colormap jet; colorbar;

%% Optical Layout: Full Model 

%--Focal Lengths
% mp.fl = 1; 

%--Pupil Plane Resolutions
mp.P1.full.Nbeam = 1000;
mp.P2.full.Nbeam = 1000;
mp.P3.full.Nbeam = 1000;
mp.P4.full.Nbeam = 200;

%--Shaped Pupil Mask.
mp.P3.full.mask = fitsread('SPM_SPC-20190130.fits');
mp.SPname = 'SPC-20190130';

%--FPM resolution
mp.F3.full.res = 20;    % sampling of FPM for full model [pixels per lambda0/D]

%--Load and downsample the FPM. To get good grayscale edges, convolve with the correct window before downsampling. 
FPM0 = fitsread('FPM_res100_SPC-20190130.fits'); %--Resolution of 100 pixels per lambda0/D
FPM0 = padOrCropOdd(FPM0,1821);
% figure(1); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
dx0 = 1/100;
dx1 = 1/mp.F3.full.res;
N0 = size(FPM0,1);
switch lower(mp.centering)
    case{'pixel'}
        N1 = ceil_odd(N0*dx0/dx1);
    case{'interpixel'}
        N1 = ceil_even(N0*dx0/dx1);
end
x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
[X0,Y0] = meshgrid(x0);
R0 = sqrt(X0.^2+Y0.^2);
Window = 0*R0;
Window(R0<=dx1/2) = 1; Window = Window/sum(sum(Window));
% figure(10); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
FPM0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(FPM0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
FPM0 = circshift(FPM0,[1 1]); %--Undo a centering shift
x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
[X1,Y1] = meshgrid(x1);
FPM1 = interp2(X0,Y0,FPM0,X1,Y1,'cubic',0); %--Downsample by interpolation
switch lower(mp.centering)
    case{'pixel'}
        mp.F3.full.mask.amp = zeros(N1+1,N1+1);
        mp.F3.full.mask.amp(2:end,2:end) = FPM1;
    otherwise
        mp.F3.full.mask.amp = FPM1;
end

% figure(2); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% figure(3); imagesc(FPM1); axis xy equal tight; colormap jet; colorbar;


%% Mask Definitions

%--Pupil definition
mp.whichPupil = 'WFIRST180718';
mp.P1.IDnorm = 0.303; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.D = 2.3631; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.

%--Lyot stop shape
mp.LSshape = 'bowtie';
mp.P4.IDnorm = 0.38; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.92; %--Lyot stop OD [Dtelescope]
mp.P4.ang = 90;      %--Lyot stop opening angle [degrees]
mp.P4.wStrut = 0;    %--Lyot stop strut width [pupil diameters]

%--FPM size
mp.F3.Rin = 2.6;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin 
mp.F3.Rout = 9;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 65;    % on each side, opening angle [degrees]


%% LC-Specific Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

