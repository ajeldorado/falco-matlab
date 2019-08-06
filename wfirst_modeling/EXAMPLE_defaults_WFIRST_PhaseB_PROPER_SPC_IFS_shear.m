
% %--Initialize some structures if they don't already exist

%% Misc

%--Record Keeping
mp.SeriesNum = 44;
mp.TrialNum = 0;

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
mp.thput_eval_x = 7; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%--Where to shift the source to compute the intensity normalization value.
mp.source_x_offset_norm = 7;  % x location [lambda_c/D] in dark hole at which to compute intensity normalization
mp.source_y_offset_norm = 0;  % y location [lambda_c/D] in dark hole at which to compute intensity normalization

%% Bandwidth and Wavelength Specs

mp.lambda0 = 730e-9;   %--Central wavelength of the whole spectral bandpass [meters]
% mp.fracBW = 0.15;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 5;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 3;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

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

mp.ctrl.flagUseModel = true; %--Whether to perform a model-based (vs empirical) grid search for the controller

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*[1]; %1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is not used; it is a placeholder for correct indexing of other modes.)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = [];%2:3;%2:6; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [];%...
%                 [3., 4.;...
%                 4., 8.;
%                 8., 9.];  

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

% % % % GRID SEARCH EFC DEFAULTS     
% %--WFSC Iterations and Control Matrix Relinearization
% mp.controller = 'gridsearchEFC';
% mp.Nitr = 5; %--Number of estimation+control iterations to perform
% mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
% mp.dm_ind = [1 2]; %--Which DMs to use

% % PLANNED SEARCH EFC DEFAULTS 
mp.controller = 'plannedEFC';
mp.dm_ind = [1 2 ]; % vector of DMs used in controller at ANY time (not necessarily all at once or all the time). 
mp.ctrl.dmfacVec = 1;
%--CONTROL SCHEDULE. Columns of mp.ctrl.sched_mat are: 
    % Column 1: # of iterations, 
    % Column 2: log10(regularization), 
    % Column 3: which DMs to use (12, 128, 129, or 1289) for control
    % Column 4: flag (0 = false, 1 = true), whether to re-linearize
    %   at that iteration.
    % Column 5: flag (0 = false, 1 = true), whether to perform an
    %   EFC parameter grid search to find the set giving the best
    %   contrast .
    % The imaginary part of the log10(regularization) in column 2 is
    %  replaced for that iteration with the optimal log10(regularization)
    % A row starting with [0, 0, 0, 1...] is for relinearizing only at that time

% mp.ctrl.sched_mat = [...
%     repmat([1,1j,12,1,1],[4,1]);...
%     repmat([1,1j-1,12,1,1],[25,1]);...
%     repmat([1,1j,12,1,1],[1,1]);...
%     ];

mp.ctrl.sched_mat = [...
    [0,0,0,1,0];
    repmat([1,1j,12,0,1],[5,1]);...
    [1,-5,12,0,0];...
    repmat([1,1j,12,0,1],[9,1]);...
    ];
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);


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
mp.dm1.ytilt = 5.7;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = 23.5;%(48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = 23.5;%(48/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 48;               % # of actuators across DM array
mp.dm2.VtoH = 1e-9*ones(48);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 5.7;%               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;              % clocking of DM surface [degrees]
mp.dm2.xc = 23.5;%(48/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = 23.5;%(48/2 - 1/2);       % y-center location of DM surface [actuator widths]
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
mp.layout = 'wfirst_phaseb_proper';  %--Which optical layout to use
mp.coro = 'SPLC';
mp.flagApod = true;    %--Whether to use an apodizer or not
mp.flagDMwfe = false;  %--Whether to use BMC DM quilting maps

%--Final Focal Plane Properties
mp.Fend.res = 3;%(730/660)*2.; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 11.; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin = 2.6;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 9;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 65;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin = 3;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = 9;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang = 65;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'both'; %--Which side(s) for correction: 'both', 'left', 'right', 'top', 'bottom'

%% Define custom SPC values.

%--Full model mask files and resolutions
file_dir = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/spc_ifs_custom/';

% mp.full.pupil_mask_file = [file_dir, 'SPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW15Nlam6.fits'];        mp.fracBW = 0.15; mp.Nsbp = 5;%--SPM file name
mp.full.pupil_mask_file = [file_dir, 'minpadSPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW15Nlam6.fits'];  mp.fracBW = 0.15; mp.Nsbp = 5;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'SPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW2Nlam6.fits'];         mp.fracBW = 0.02; mp.Nsbp = 1;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'minpadSPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW2Nlam6.fits'];   mp.fracBW = 0.02; mp.Nsbp = 1;%--SPM file name

% mp.full.pupil_file = [file_dir, 'pupil_SPC-20190130_rotated.fits'];
mp.full.pupil_file = [file_dir, 'unpaddedpupil_full_symm_N1000_rotated.fits'];
mp.full.pupil_diam_pix = 1000;

%--Add pupil shear
shear_x_pix = 1;%10; % [pixels]
shear_y_pix = 0;  % [pixels]
pupil0 = fitsread(mp.full.pupil_file);
Npad = size(pupil0,1)+ceil_even(max([shear_x_pix,shear_y_pix])*3);
pupil = padOrCropEven(pupil0,Npad);
pupil = circshift(pupil,[shear_y_pix,shear_x_pix]);
mp.full.pupil_file = [mp.full.pupil_file(1:end-5) '_sheared.fits'];
fitswrite(pupil,mp.full.pupil_file); %--Write the sheared pupil to a file
figure(110); imagesc(pupil-padOrCropEven(pupil0,Npad)); axis xy equal tight; drawnow;

mp.full.fpm_file = [file_dir, 'fpm_sharp_26WA90_65deg_res20.fits'];
mp.full.fpm_sampling_lam0 = 0.05; 	% sampling in lambda0/D of FPM mask
mp.full.lyot_stop_file = [file_dir,'LS_sharp_20D96_95deg_N1000.fits'];
mp.full.lambda0_m = mp.lambda0;       % FPM scaled for this central wavelength

%--Compact model FPM
mp.compact.flagGenFPM = true;
mp.F3.Rin = 2.6;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin 
mp.F3.Rout = 9;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 65;    % on each side, opening angle [degrees]

%--Lyot stop shape
mp.compact.flagGenLS = true;
mp.LSshape = 'bowtie';
mp.P4.IDnorm = 0.20; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.96; %--Lyot stop OD [Dtelescope]
mp.P4.ang = 95;      %--Lyot stop opening angle [degrees]
mp.P4.wStrut = 0;    %--Lyot stop strut width [pupil diameters]


%% Mask Definitions (and generation flags)

%--Pupil definition
mp.compact.flagGenPupil = false;
mp.whichPupil = 'WFIRST180718';
mp.P1.IDnorm = 0.303; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
mp.P1.D = 2.3631; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
mp.P1.Dfac = 1; %--Factor scaling inscribed OD to circumscribed OD for the telescope pupil.

%--Lyot stop shape
mp.compact.flagGenLS = true;
mp.LSshape = 'bowtie';
mp.P4.IDnorm = 0.20; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.96; %--Lyot stop OD [Dtelescope]
mp.P4.ang = 95;      %--Lyot stop opening angle [degrees]
mp.P4.wStrut = 0;    %--Lyot stop strut width [pupil diameters]

%--FPM size
mp.compact.flagGenFPM = true;
mp.F3.Rin = 2.6;   % inner hard-edge radius of the focal plane mask [lambda0/D]. Needs to be <= mp.F3.Rin 
mp.F3.Rout = 9;   % radius of outer opaque edge of FPM [lambda0/D]
mp.F3.ang = 65;    % on each side, opening angle [degrees]

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1.; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 46.3e-3; %46.2987e-3;
mp.P3.D = 46.3e-3; %46.2987e-3;
mp.P4.D = 46.3e-3; %46.2987e-3;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 386;
mp.P2.compact.Nbeam = 386;
mp.P3.compact.Nbeam = 386;
mp.P4.compact.Nbeam = 60;

%%
%--Pupil 
% mp.P1.compact.mask = fitsread(mp.full.pupil_file);
% mp.P1.compact.mask(2:end,2:end) = rot90(mp.P1.compact.mask(2:end,2:end),2);

%--Pupil Mask: Load and downsample.
pupil0 = fitsread(mp.full.pupil_file);
pupil0(2:end,2:end) = rot90(pupil0(2:end,2:end),2);

if(mp.P1.compact.Nbeam==1000)
    mp.P1.compact.mask = pupil0;
else
    pupil0 = pupil0(2:end,2:end);
    % figure(1); imagesc(pupil0); axis xy equal tight; colormap jet; colorbar;
    % figure(11); imagesc(pupil0-fliplr(pupil0)); axis xy equal tight; colormap jet; colorbar;
    dx0 = 1/1000;
    dx1 = 1/mp.P1.compact.Nbeam;
    N0 = size(pupil0,1);
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
    pupil0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(pupil0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
    pupil0 = circshift(pupil0,[1 1]); %--Undo a centering shift
    x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
    [X1,Y1] = meshgrid(x1);
    pupil1 = interp2(X0,Y0,pupil0,X1,Y1,'cubic',0); %--Downsample by interpolation

    switch lower(mp.centering)
        case{'pixel'}
            mp.P1.compact.mask = zeros(N1+1,N1+1);
            mp.P1.compact.mask(2:end,2:end) = pupil1;
        otherwise
            mp.P1.compact.mask = pupil1;
    end
    
%     figure(2); imagesc(pupil0); axis xy equal tight; colormap jet; colorbar; drawnow;
%     figure(3); imagesc(pupil1); axis xy equal tight; colormap jet; colorbar; drawnow;
%     figure(12); imagesc(pupil0-fliplr(pupil0)); axis xy equal tight; colormap jet; colorbar;
%     figure(13); imagesc(pupil1-fliplr(pupil1)); axis xy equal tight; colormap jet; colorbar;
end


%%
%--Shaped Pupil Mask: Load and downsample.
mp.SPname = 'SPC-custom';
SP0 = fitsread(mp.full.pupil_mask_file);
SP0 = padOrCropEven(SP0,Npad);
% SP0(2:end,2:end) = rot90(SP0(2:end,2:end),2);

if(mp.P1.compact.Nbeam==1000)
    mp.P3.compact.mask = SP0;
else
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
    
    figure(2); imagesc(SP0); axis xy equal tight; colormap jet; colorbar; drawnow;
    figure(3); imagesc(SP1); axis xy equal tight; colormap jet; colorbar; drawnow;
    % figure(12); imagesc(SP0-fliplr(SP0)); axis xy equal tight; colormap jet; colorbar;
    % figure(13); imagesc(SP1-fliplr(SP1)); axis xy equal tight; colormap jet; colorbar;

end


%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 1; %--How many times to rotate the final image by 180 degrees

%--FPM resolution
mp.F3.compact.res = 6;    % sampling of FPM for compact model [pixels per lambda0/D]

% %--Load and downsample the FPM. To get good grayscale edges, convolve with the correct window before downsampling. 
% FPM0 = fitsread('FPM_res100_SPC-20190130.fits'); %--Resolution of 100 pixels per lambda0/D
% FPM0 = padOrCropOdd(FPM0,1821);
% % figure(1); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% % figure(11); imagesc(FPM0-rot90(FPM0,2)); axis xy equal tight; colormap jet; colorbar;
% dx0 = 1/100;
% dx1 = 1/mp.F3.compact.res;
% N0 = size(FPM0,1);
% switch lower(mp.centering)
%     case{'pixel'}
%         N1 = ceil_odd(N0*dx0/dx1);
%     case{'interpixel'}
%         N1 = ceil_even(N0*dx0/dx1);
% end
% x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
% [X0,Y0] = meshgrid(x0);
% R0 = sqrt(X0.^2+Y0.^2);
% Window = 0*R0;
% Window(R0<=dx1/2) = 1; Window = Window/sum(sum(Window));
% % figure(10); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
% FPM0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(FPM0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
% FPM0 = circshift(FPM0,[1 1]); %--Undo a centering shift
% x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
% [X1,Y1] = meshgrid(x1);
% FPM1 = interp2(X0,Y0,FPM0,X1,Y1,'cubic',0); %--Downsample by interpolation
% switch lower(mp.centering)
%     case{'pixel'}
%         mp.F3.compact.mask.amp = zeros(N1+1,N1+1);
%         mp.F3.compact.mask.amp(2:end,2:end) = FPM1;
%     otherwise
%         mp.F3.compact.mask.amp = FPM1;
% end
% figure(2); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% figure(3); imagesc(FPM1); axis xy equal tight; colormap jet; colorbar;
% figure(12); imagesc(FPM0-rot90(FPM0,2)); axis xy equal tight; colormap jet; colorbar;
% figure(13); imagesc(FPM1-rot90(FPM1,2)); axis xy equal tight; colormap jet; colorbar;

%% Optical Layout: Full Model 

mp.full.data_dir = '/Users/ajriggs/Repos/proper-models/wfirst_cgi/data_phaseb/'; % mask design data path
mp.full.cor_type = 'spc_ifs_custom'; %   'hlc', 'spc', or 'none' (none = clear aperture, no coronagraph)

mp.full.flagGenFPM = false;
mp.full.flagPROPER = true; %--Whether the full model is a PROPER prescription

% %--Pupil Plane Resolutions
mp.P1.full.Nbeam = 1000;
mp.P1.full.Narr = 1002;

mp.full.output_dim = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D

mp.full.pol_conds = 10;% [-2,-1,1,2]; %--Which polarization states to use when creating an image.
mp.full.polaxis = 10;                %   polarization condition (only used with input_field_rootname)
mp.full.use_errors = 1;

mp.full.zindex = 4;
mp.full.zval_m = 0.19e-9;
mp.full.use_hlc_dm_patterns = false; % whether to use design WFE maps for HLC
mp.full.lambda0_m = mp.lambda0;
mp.full.input_field_rootname = '';	%   rootname of files containing aberrated pupil
% mp.full.use_hlc_dm_patterns = 1;	%   use Dwight-generated HLC default DM wavefront patterns? 1 or 0
mp.full.use_dm1 = 0;                %   use DM1? 1 or 0
mp.full.use_dm2 = 0;                %   use DM2? 1 or 0
mp.full.dm_sampling_m = 0.9906e-3;     %   actuator spacing in meters; default is 1 mm
mp.full.dm1_xc_act = 23.5;          %   for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
mp.full.dm1_yc_act = 23.5;
mp.full.dm1_xtilt_deg = 0;   		%   tilt around X axis
mp.full.dm1_ytilt_deg = 5.7;		%   effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
mp.full.dm1_ztilt_deg = 0;
mp.full.dm2_xc_act = 23.5;		
mp.full.dm2_yc_act = 23.5;
mp.full.dm2_xtilt_deg = 0;   
mp.full.dm2_ytilt_deg = 5.7;
mp.full.dm2_ztilt_deg = 0;
mp.full.use_fpm  = 1;
mp.full.fpm_axis = 'p';             %   HLC FPM axis: '', 's', 'p'

mp.full.dm1.flatmap = fitsread('errors_polaxis10_dm.fits');
mp.full.dm2.flatmap = 0;



% %--Pupil Plane Resolutions
% mp.P1.full.Nbeam = 1000;
% mp.P2.full.Nbeam = 1000;
% mp.P3.full.Nbeam = 1000;
% mp.P4.full.Nbeam = 200;

% %--Shaped Pupil Mask.
% mp.P3.full.mask = fitsread('SPM_SPC-20190130.fits');
% mp.SPname = 'SPC-20190130';

% %--FPM resolution
% mp.F3.full.res = 20;    % sampling of FPM for full model [pixels per lambda0/D]
% %--Load and downsample the FPM. To get good grayscale edges, convolve with the correct window before downsampling. 
% FPM0 = fitsread('FPM_res100_SPC-20190130.fits'); %--Resolution of 100 pixels per lambda0/D
% FPM0 = padOrCropOdd(FPM0,1821);
% % figure(1); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% dx0 = 1/100;
% dx1 = 1/mp.F3.full.res;
% N0 = size(FPM0,1);
% switch lower(mp.centering)
%     case{'pixel'}
%         N1 = ceil_odd(N0*dx0/dx1);
%     case{'interpixel'}
%         N1 = ceil_even(N0*dx0/dx1);
% end
% x0 = (-(N0-1)/2:(N0-1)/2)*dx0;
% [X0,Y0] = meshgrid(x0);
% R0 = sqrt(X0.^2+Y0.^2);
% Window = 0*R0;
% Window(R0<=dx1/2) = 1; Window = Window/sum(sum(Window));
% % figure(10); imagesc(Window); axis xy equal tight; colormap jet; colorbar;
% FPM0 = ifftshift(  ifft2( fft2(fftshift(Window)).*fft2(fftshift(FPM0)) )); %--To get good grayscale edges, convolve with the correct window before downsampling.
% FPM0 = circshift(FPM0,[1 1]); %--Undo a centering shift
% x1 = (-(N1-1)/2:(N1-1)/2)*dx1;
% [X1,Y1] = meshgrid(x1);
% FPM1 = interp2(X0,Y0,FPM0,X1,Y1,'cubic',0); %--Downsample by interpolation
% switch lower(mp.centering)
%     case{'pixel'}
%         mp.F3.full.mask.amp = zeros(N1+1,N1+1);
%         mp.F3.full.mask.amp(2:end,2:end) = FPM1;
%     otherwise
%         mp.F3.full.mask.amp = FPM1;
% end
% % figure(2); imagesc(FPM0); axis xy equal tight; colormap jet; colorbar;
% % figure(3); imagesc(FPM1); axis xy equal tight; colormap jet; colorbar;




%% LC-Specific Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

