% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to perform wavefront control with the WFIRST CGI Phase B SP(L)C-IFS design.
%  1) Load the default model parameters for an SPLC.
%  2) Specify the values to overwrite.
%  3) Run a single trial of WFC using FALCO.
%
% REVISION HISTORY:
% --------------
% Created on 2019-04-17 by A.J. Riggs.
% ---------------

clear all;


%% Step 1: Define Necessary Paths on Your Computer System

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


%% Step 2: Load default model parameters

EXAMPLE_defaults_WFIRST_PhaseB_SPC_IFS_simple_shear

%% Step 3: Overwrite default values as desired

% %%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;
mp.propMethodPTP = 'mft';

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_SeriesX_TrialY.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING:
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass
% % mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

% % % GRID-SEARCH EFC     
mp.Nitr = 10; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1;%1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.ctrl.flagUseModel = true; %--Use the compact model for the grid search (much faster)

%% Define Pupil and Shaped Pupil (full model). Change this cell

mp.full.flagGenPupil = false;
mp.full.flagGenApod = false;

file_dir = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/spc_ifs_custom/';

%--FULL MODEL: PUPIL
% mp.full.pupil_file = [file_dir, 'pupil_SPC-20190130_rotated.fits'];
pupil_file = [file_dir, 'unpaddedpupil_full_symm_N1000_rotated.fits'];
% mp.P1.full.mask = fitsread(pupil_file);

%--Add shear to pupil mask
shear_x_pix = 1;%10; % [pixels] 
shear_y_pix = 0;  % [pixels]
Npad = 1050; %--Make sure is large enough to include the sheared pupil
pupil0 = fitsread(pupil_file);
pupil0(2:end,2:end) = rot90(pupil0(2:end,2:end),2); %--Derotate pupil
pupil = padOrCropEven(pupil0,Npad);
pupil = circshift(pupil,[shear_y_pix,shear_x_pix]);
mp.P1.full.mask = pupil;


pupil_file = [pupil_file(1:end-5) '_sheared.fits'];
fitswrite(pupil,pupil_file); %--Write the sheared pupil to a file
if(mp.flagPlot);  figure(110); imagesc(pupil-padOrCropEven(pupil0,Npad)); axis xy equal tight; drawnow;  end

%--FULL MODEL: SPM
% mp.full.pupil_mask_file = [file_dir, 'SPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW15Nlam6.fits'];        mp.fracBW = 0.15; mp.Nsbp = 5;%--SPM file name
sp_file = [file_dir, 'minpadSPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW15Nlam6.fits'];  mp.fracBW = 0.15; mp.Nsbp = 5;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'SPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW2Nlam6.fits'];         mp.fracBW = 0.02; mp.Nsbp = 1;%--SPM file name
% mp.full.pupil_mask_file = [file_dir, 'minpadSPM_jg36_79c81_PH40_65deg_26WA90_20LS96_RoC1_LS95deg_BW2Nlam6.fits'];   mp.fracBW = 0.02; mp.Nsbp = 1;%--SPM file name
% mp.P3.full.mask = fitsread('SPM_SPC-20190130.fits');
 mp.P3.full.mask = fitsread(sp_file);
 mp.P3.full.mask = padOrCropEven(mp.P3.full.mask,Npad);

%% Downsample the Pupil and SPM

mp.compact.flagGenPupil = false;
mp.compact.flagGenApod = false;

%--DOWNSAMPLE PUPIL
pupil0 = pupil;
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
end
    
    

%--DOWNSAMPLE SPM
SP0 = mp.P3.full.mask;
if(mp.P1.compact.Nbeam==1000)
    mp.P3.compact.mask = SP0;
else
    SP0 = SP0(2:end,2:end);
    % figure(1); imagesc(SP0); axis xy equal tight; colormap jet; colorbar;
    % figure(11); imagesc(SP0-fliplr(SP0)); axis xy equal tight; colormap jet; colorbar;
    dx0 = 1/1000;
    dx1 = 1/mp.P3.compact.Nbeam;
    %--Use N0 from earlier %N0 = size(SP0,1);
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

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

out = falco_wfsc_loop(mp);
