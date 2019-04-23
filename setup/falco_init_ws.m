% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to provide input parameters as structures
% Setup parameters for a mock SPLC at the HCIT.
%
% Created by A.J. Riggs on 2017-10-31.

function [mp,out] = falco_init_ws(fn_config)

%% Read inputs as structures from a .mat config file

load(fn_config,'mp');

mainPath = mp.path.falco;

disp(['DM 1-to-2 Fresnel number (using radius) = ',num2str((mp.P2.D/2)^2/(mp.d_dm1_dm2*mp.lambda0))]);

%% Intializations of structures (if they don't exist yet)
mp.jac.dummy = 1;

%% Optional/Hidden flags
if(isfield(mp,'flagSaveWS')==false);  mp.flagSaveWS = false;  end  %--Whehter to save otu the entire workspace at the end of the trial. Can take up lots of space.
if(isfield(mp,'flagSaveEachItr')==false);  mp.flagSaveEachItr = false;  end  %--Whether to save out the performance at each iteration. Useful for long trials in case it crashes or is stopped early.
if(isfield(mp,'flagSVD')==false);  mp.flagSVD = false;  end    %--Whether to compute and save the singular mode spectrum of the control Jacobian (each iteration)
if(isfield(mp,'flagFiber')==false);  mp.flagFiber = false;  end  %--Whether to couple the final image through lenslets and a single mode fiber.
if(isfield(mp,'flagDMwfe')==false);  mp.flagDMwfe = false;  end  %--Temporary for BMC quilting study
if(isfield(mp,'flagTrainModel')==false);  mp.flagTrainModel = false;  end  %--Whether to call the Expectation-Maximization (E-M) algorithm to improve the linearized model. 
if(isfield(mp,'flagUseLearnedJac')==false);  mp.flagUseLearnedJac = false;  end  %--Whether to load and use an improved Jacobian from the Expectation-Maximization (E-M) algorithm 
if(isfield(mp.est,'flagUseJac')==false); mp.est.flagUseJac = false; end   %--Whether to use the Jacobian or not for estimation. (If not using Jacobian, model is called and differenced.)

%% Optional/Hidden variables
if(isfield(mp,'propMethodPTP')==false);  mp.propMethodPTP = 'fft';  end %--Propagation method for postage stamps around the influence functions. 'mft' or 'fft'
if(isfield(mp,'SPname')==false);  mp.SPname = 'none';  end %--Apodizer name default
%--Training Data: mp.NitrTrain = 5;  %--The number of correction iterations to use per round of training data for the adaptive Jacobian (E-M) algorithm.
%--Zernike sensitivities to 1nm RMS: which noll indices in which annuli, given by mp.eval.indsZnoll and mp.eval.Rsens 
%--Tied actuator pair definitions: See Section with variables mp.dmX.tied for X=1:9
%--Quantization of DM actuation steps based on least significant bit of the
% DAC (digital-analog converter). In height, so called HminStep. If HminStep (minimum step in H) is defined, then quantize the DM voltages
% Variables to define if wanted: mp.dm1.HminStep, mp.dm2.HminStep

%% File Paths

%--Storage directories (data in these folders will not be synced via Git
if(isfield(mp.path,'ws')==false); mp.path.ws = [mainPath 'data' filesep 'ws' filesep]; end % Store final workspace data here
if(isfield(mp.path,'maps')==false); mp.path.falcoaps = [mainPath 'maps' filesep]; end % Maps go here
if(isfield(mp.path,'jac')==false); mp.path.jac = [mainPath 'data' filesep 'jac' filesep]; end % Store the control Jacobians here
if(isfield(mp.path,'images')==false); mp.path.images = [mainPath 'data' filesep 'images' filesep]; end % Store all full, reduced images here
if(isfield(mp.path,'dm')==false); mp.path.dm = [mainPath 'data' filesep 'DM' filesep]; end % Store DM command maps here
if(isfield(mp.path,'ws_inprogress')==false); mp.path.ws_inprogress = [mainPath 'data' filesep 'ws_inprogress' filesep]; end % Store in progress workspace data here

%% Loading previous DM commands as the starting point
%--Stash DM8 and DM9 starting commands if they are given in the main script
if(isfield(mp,'dm8'))
    if(isfield(mp.dm8,'V')); mp.DM8V0 = mp.dm8.V; end
    if(isfield(mp.dm9,'V')); mp.DM9V0 = mp.dm9.V; end
end

%% Useful factor
mp.mas2lam0D = 1/(mp.lambda0/mp.P1.D*180/pi*3600*1000); %--Conversion factor: milliarcseconds (mas) to lambda0/D

%% Estimator
if(isfield(mp,'estimator')==false); mp.estimator = 'perfect'; end

%% Bandwidth and Wavelength Specs

if(isfield(mp,'Nwpsbp')==false)
    mp.Nwpsbp = 1;
end
mp.full.Nlam = mp.Nsbp*mp.Nwpsbp; %--Total number of wavelengths in the full model

%--When in simulation and using perfect estimation, use end wavelengths in bandbass, which (currently) requires Nwpsbp=1. 
if(strcmpi(mp.estimator,'perfect') && mp.Nsbp>1)
    if(mp.Nwpsbp>1)
        fprintf('* Forcing mp.Nwpsbp = 1 * \n')
        mp.Nwpsbp = 1; % number of wavelengths per sub-bandpass. To approximate better each finite sub-bandpass in full model with an average of images at these values. Only >1 needed when each sub-bandpass is too large (say >3%).
    end
end

%--Center-ish wavelength indices (ref = reference)
mp.si_ref = ceil(mp.Nsbp/2);
mp.wi_ref = ceil(mp.Nwpsbp/2);

%--Wavelengths used for Compact Model (and Jacobian Model)
mp.sbp_weights = ones(mp.Nsbp,1);
if(strcmpi(mp.estimator,'perfect') && mp.Nsbp>1) %--For design or modeling without estimation: Choose ctrl wvls evenly between endpoints (inclusive) of the total bandpass
    mp.fracBWsbp = mp.fracBW/(mp.Nsbp-1);
    mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
else %--For cases with estimation: Choose est/ctrl wavelengths to be at subbandpass centers.
    mp.fracBWsbp = mp.fracBW/mp.Nsbp;
    mp.fracBWcent2cent = mp.fracBW*(1-1/mp.Nsbp); %--Bandwidth between centers of endpoint subbandpasses.
    mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBWcent2cent/2,1+mp.fracBWcent2cent/2,mp.Nsbp); %--Space evenly at the centers of the subbandpasses.
end
mp.sbp_weights = mp.sbp_weights/sum(mp.sbp_weights); %--Normalize the sum of the weights

fprintf(' Using %d discrete wavelength(s) in each of %d sub-bandpasses over a %.1f%% total bandpass \n', mp.Nwpsbp, mp.Nsbp,100*mp.fracBW);
fprintf('Sub-bandpasses are centered at wavelengths [nm]:\t '); fprintf('%.2f  ',1e9*mp.sbp_centers); fprintf('\n\n');

%--Wavelength factors/weights within sub-bandpasses in the full model
mp.full.lambda_weights = ones(mp.Nwpsbp,1); %--Initialize as all ones. Weights within a single sub-bandpass
if(mp.Nsbp==1)
    mp.full.sbp_facs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nwpsbp);
    if(mp.Nwpsbp>2) %--Include end wavelengths with half weights
        mp.full.lambda_weights(1) = 1/2;
        mp.full.lambda_weights(end) = 1/2;
    end
else %--For cases with estimation (est/ctrl wavelengths at subbandpass centers). Full model only
    mp.full.sbp_facs = linspace(1-(mp.fracBWsbp/2)*(1-1/mp.Nwpsbp),...
                           1+(mp.fracBWsbp/2)*(1-1/mp.Nwpsbp), mp.Nwpsbp);
end
if(mp.Nwpsbp==1);  mp.full.sbp_facs = 1;  end %--Set factor to 1 if only 1 value.

mp.full.lambda_weights = mp.full.lambda_weights/sum(mp.full.lambda_weights); %--Normalize sum of the weights

%--Make vector of all wavelengths and weights used in the full model
mp.full.lambdas = zeros(mp.Nsbp*mp.Nwpsbp,1);
mp.full.weights = zeros(mp.Nsbp*mp.Nwpsbp,1);
counter = 1;
for si=1:mp.Nsbp
    for wi=1:mp.Nwpsbp
        mp.full.lambdas(counter) = mp.sbp_centers(si)*mp.full.sbp_facs(wi);
        mp.full.all_weights = mp.sbp_weights(si)*mp.full.lambda_weights(wi);
        counter = counter+1;
    end
end

%% Zernike and Chromatic Weighting of the Control Jacobian
if(isfield(mp.jac,'zerns')==false);  mp.jac.zerns = 1; end %--Which Zernike modes to include in Jacobian [Noll index]. Always include 1 for piston term.
if(isfield(mp.jac,'Zcoef')==false);  mp.jac.Zcoef = 1e-9*ones(10,1); end %--meters RMS of Zernike aberrations. (piston value is reset to 1 later for correct normalization)

mp = falco_config_jac_weights(mp);

%% Pupil Masks
mp = falco_config_gen_chosen_pupil(mp); %--input pupil mask
mp = falco_config_gen_chosen_apodizer(mp); %--apodizer mask
mp = falco_config_gen_chosen_LS(mp); %--Lyot stop

%% Plot the pupil and Lyot stop on top of each other to make sure they are aligned correctly
%--Only for coronagraphs using Babinet's principle, for which the input
%pupil and Lyot plane have the same resolution.
switch upper(mp.coro)
    case{'FOHLC','HLC','LC','APLC','VC','AVC'}
        if(mp.flagPlot)
            P4mask = padOrCropEven(mp.P4.compact.mask,mp.P1.compact.Narr);
            P4mask = rot90(P4mask,2);
            if(strcmpi(mp.centering,'pixel'))
               P4mask = circshift(P4mask,[1 1]); 
            end
            P1andP4 = mp.P1.compact.mask + P4mask;
            figure(301); imagesc(P1andP4); axis xy equal tight; colorbar; set(gca,'Fontsize',20); title('Pupil and LS Superimposed','Fontsize',16');
            
            if(mp.flagApod)
                P1andP3 = mp.P1.compact.mask + padOrCropEven(mp.P3.compact.mask,length(mp.P1.compact.mask));
                figure(302); imagesc(P1andP3); axis xy equal tight; colorbar; set(gca,'Fontsize',20); title('Pupil and Apod Superimposed','Fontsize',16');
            end
        end
end

%% DM Initialization

%--Initialize the number of actuators (NactTotal) and actuators used (Nele).
mp.dm1.NactTotal=0; mp.dm2.NactTotal=0; mp.dm3.NactTotal=0; mp.dm4.NactTotal=0; mp.dm5.NactTotal=0; mp.dm6.NactTotal=0; mp.dm7.NactTotal=0; mp.dm8.NactTotal=0; mp.dm9.NactTotal=0; %--Initialize for bookkeeping later.
mp.dm1.Nele=0; mp.dm2.Nele=0; mp.dm3.Nele=0; mp.dm4.Nele=0; mp.dm5.Nele=0; mp.dm6.Nele=0; mp.dm7.Nele=0; mp.dm8.Nele=0; mp.dm9.Nele=0; %--Initialize for Jacobian calculations later. 

%% HLC and EHLC FPM: Initialization and Generation

switch lower(mp.layout)
    case 'fourier'
        switch upper(mp.coro)
            case{'HLC'}
                switch mp.dm9.inf0name
                    case '3foldZern'
                        mp = falco_setup_FPM_HLC_3foldZern(mp);
                    otherwise
                        mp = falco_setup_FPM_HLC(mp);
                end
                mp = falco_config_gen_FPM_HLC(mp);
            case{'FOHLC'}
                mp = falco_setup_FPM_FOHLC(mp);
                mp = falco_config_gen_FPM_FOHLC(mp);
                mp.compact.Nfpm = max([mp.dm8.compact.NdmPad,mp.dm9.compact.NdmPad]); %--Width of the FPM array in the compact model.
                mp.full.Nfpm = max([mp.dm8.NdmPad,mp.dm9.NdmPad]); %--Width of the FPM array in the full model.
            case{'EHLC'}
                mp = falco_setup_FPM_EHLC(mp);
                mp = falco_config_gen_FPM_EHLC(mp);
            case 'SPHLC'
                mp = falco_config_gen_FPM_SPHLC(mp);
        end

        %%--Pre-compute the complex transmission of the allowed Ni+PMGI FPMs.
        switch upper(mp.coro)
            case{'EHLC','HLC','SPHLC'}
                [mp.complexTransCompact,mp.complexTransFull] = falco_gen_complex_trans_table(mp);
        end
end

%% Generate FPM

switch upper(mp.coro)
    case {'LC','APLC'} %--Occulting spot FPM (can be HLC-style and partially transmissive)
        mp = falco_config_gen_FPM_LC(mp);
    case{'SPLC','FLC'}
        mp = falco_config_gen_FPM_SPLC(mp);
    case{'RODDIER'}
        mp = falco_config_gen_FPM_Roddier(mp);  
end

%% FPM coordinates, [meters] and [dimensionless]

switch upper(mp.coro)
    case{'VORTEX','VC','AVC'}   %--Nothing needed to run the vortex model
    case 'SPHLC' %--Moved to separate function
    otherwise
        
        switch mp.layout
            case{'wfirst_phaseb_simple','wfirst_phaseb_proper'}
            otherwise
                %--FPM (at F3) Resolution [meters]
                mp.F3.full.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res;
                mp.F3.full.deta = mp.F3.full.dxi;
        end
        mp.F3.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res;
        mp.F3.compact.deta = mp.F3.compact.dxi;
        
        %--Coordinates in FPM plane in the compact model [meters]
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.compact.Nxi,2)==1)
            mp.F3.compact.xis  = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)*mp.F3.compact.dxi;
            mp.F3.compact.etas = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2).'*mp.F3.compact.deta;
        else
            mp.F3.compact.xis  = (-mp.F3.compact.Nxi/2: (mp.F3.compact.Nxi/2-1))*mp.F3.compact.dxi;
            mp.F3.compact.etas = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1)).'*mp.F3.compact.deta;
        end

        switch mp.layout
            case{'wfirst_phaseb_simple','wfirst_phaseb_proper'}
            otherwise
                %--Coordinates (dimensionless [DL]) for the FPMs in the full model
                if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.Nxi,2)==1)
                    mp.F3.full.xisDL  = (-(mp.F3.full.Nxi-1)/2:(mp.F3.full.Nxi-1)/2)/mp.F3.full.res;
                    mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
                else
                    mp.F3.full.xisDL  = (-mp.F3.full.Nxi/2:(mp.F3.full.Nxi/2-1))/mp.F3.full.res;
                    mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
                end
        end
        
        %--Coordinates (dimensionless [DL]) for the FPMs in the compact model
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.compact.Nxi,2)==1)
            mp.F3.compact.xisDL  = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2)/mp.F3.compact.res;
        else
            mp.F3.compact.xisDL  = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1))/mp.F3.compact.res;
        end
end

%% Sampling/Resolution and Scoring/Correction Masks for Final Focal Plane (Fend.

mp.Fend.dxi = (mp.fl*mp.lambda0/mp.P4.D)/mp.Fend.res; % sampling at Fend.[meters]
mp.Fend.deta = mp.Fend.dxi; % sampling at Fend.[meters]    

if(mp.flagFiber)
    mp.Fend.lenslet.D = 2*mp.Fend.res*mp.Fend.lensletWavRad*mp.Fend.dxi;
    mp.Fend.x_lenslet_phys = mp.Fend.dxi*mp.Fend.res*mp.Fend.x_lenslet;
    mp.Fend.y_lenslet_phys = mp.Fend.deta*mp.Fend.res*mp.Fend.y_lenslet;

    mp.F5.dxi = mp.lensletFL*mp.lambda0/mp.Fend.lenslet.D/mp.F5.res;
    mp.F5.deta = mp.F5.dxi;
end
    
%% Software Mask for Correction (corr) and Scoring (score)

%--Set Inputs
maskCorr.pixresFP = mp.Fend.res;
maskCorr.rhoInner = mp.Fend.corr.Rin; %--lambda0/D
maskCorr.rhoOuter = mp.Fend.corr.Rout ; %--lambda0/D
maskCorr.angDeg = mp.Fend.corr.ang; %--degrees
maskCorr.centering = mp.centering;
maskCorr.FOV = mp.Fend.FOV;
maskCorr.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
if(isfield(mp.Fend,'shape'));  maskCorr.shape = mp.Fend.shape;  end

%--Compact Model: Generate Software Mask for Correction 
[mp.Fend.corr.mask, mp.Fend.xisDL, mp.Fend.etasDL] = falco_gen_SW_mask(maskCorr); 
mp.Fend.corr.settings = maskCorr; %--Store values for future reference
%--Size of the output image 
%--Need the sizes to be the same for the correction and scoring masks
mp.Fend.Nxi  = size(mp.Fend.corr.mask,2);
mp.Fend.Neta = size(mp.Fend.corr.mask,1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%--Mask defining the area covered by the lenslet.  Only the immediate area
%around the lenslet is propagated, saving computation time.  This lenslet
%can then be moved around to different positions in Fend.
if(mp.flagFiber)
    maskLenslet.pixresFP = mp.Fend.res;
    maskLenslet.rhoInner = 0;
    maskLenslet.rhoOuter = mp.Fend.lensletWavRad;
    maskLenslet.angDeg = mp.Fend.corr.ang;
    maskLenslet.centering = mp.centering;
    maskLenslet.FOV = mp.Fend.FOV;
    maskLenslet.whichSide = mp.Fend.sides;
    [mp.Fend.lenslet.mask, ~, ~] = falco_gen_SW_mask(maskLenslet);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    %--Dummy mask to calculate the F5 coordinates correctly.
    maskF5.pixresFP = mp.F5.res;
    maskF5.rhoInner = 0;
    maskF5.rhoOuter = 1.22;
    maskF5.angDeg = 180;
    maskF5.centering = mp.centering;
    maskF5.FOV = mp.F5.FOV;
    maskF5.whichSide = mp.Fend.sides;
    [mp.F5.mask, mp.F5.xisDL, mp.F5.etasDL] = falco_gen_SW_mask(maskF5);

    %--Size of the output image in F5
    mp.F5.Nxi = size(mp.F5.mask, 2);
    mp.F5.Neta = size(mp.F5.mask, 1);

    %% Set up the fiber mode in F5

    V = 2*pi/mp.lambda0*mp.fiber.a*mp.fiber.NA; 
    W = 1.1428*V - 0.996;
    U = sqrt(V.^2 - W.^2);

    maskFiberCore.pixresFP = mp.F5.res;
    maskFiberCore.rhoInner = 0;
    maskFiberCore.rhoOuter = mp.fiber.a;
    maskFiberCore.angDeg = 180;
    maskFiberCore.FOV = mp.F5.FOV;
    maskFiberCore.whichSide = mp.Fend.sides;
    [mp.F5.fiberCore.mask, ~, ~] = falco_gen_SW_mask(maskFiberCore);

    maskFiberCladding.pixresFP = mp.F5.res;
    maskFiberCladding.rhoInner = mp.fiber.a;
    maskFiberCladding.rhoOuter = 10;
    maskFiberCladding.angDeg = 180;
    maskFiberCladding.FOV = mp.F5.FOV;
    maskFiberCladding.whichSide = mp.Fend.sides;
    [mp.F5.fiberCladding.mask, ~, ~] = falco_gen_SW_mask(maskFiberCladding);

	[F5XIS, F5ETAS] = meshgrid(mp.F5.xisDL, mp.F5.etasDL);

    mp.F5.RHOS = sqrt((F5XIS - mp.F5.fiberPos(1)).^2 + (F5ETAS - mp.F5.fiberPos(2)).^2);
    mp.F5.fiberCore.mode = mp.F5.fiberCore.mask.*besselj(0, U*mp.F5.RHOS/mp.fiber.a)/besselj(0,U);
    mp.F5.fiberCladding.mode = mp.F5.fiberCladding.mask.*besselk(0, W*mp.F5.RHOS/mp.fiber.a)/besselk(0,W);
    mp.F5.fiberCladding.mode(isnan(mp.F5.fiberCladding.mode)) = 0;
    mp.F5.fiberMode = mp.F5.fiberCore.mode + mp.F5.fiberCladding.mode;
    fiberModeNorm = sqrt(sum(sum(abs(mp.F5.fiberMode).^2)));
    mp.F5.fiberMode = mp.F5.fiberMode./fiberModeNorm;
end
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%--Evaluation Model for Computing Throughput (same as Compact Model but
% with different Fend.resolution)
mp.Fend.eval.dummy = 1; %--Initialize the structure if it doesn't exist.
if(isfield(mp.Fend.eval,'res')==false);  mp.Fend.eval.res = 10;  end 
maskCorr.pixresFP = mp.Fend.eval.res; %--Assign the resolution
[mp.Fend.eval.mask, mp.Fend.eval.xisDL, mp.Fend.eval.etasDL] = falco_gen_SW_mask(maskCorr);  %--Generate the mask
mp.Fend.eval.Nxi  = size(mp.Fend.eval.mask,2);
mp.Fend.eval.Neta = size(mp.Fend.eval.mask,1);
mp.Fend.eval.dxi = (mp.fl*mp.lambda0/mp.P4.D)/mp.Fend.eval.res; % higher sampling at Fend.for evaulation [meters]
mp.Fend.eval.deta = mp.Fend.eval.dxi; % higher sampling at Fend.for evaulation [meters]   

% (x,y) location [lambda_c/D] in dark hole at which to evaluate throughput
[XIS,ETAS] = meshgrid(mp.Fend.eval.xisDL - mp.thput_eval_x, mp.Fend.eval.etasDL - mp.thput_eval_y);
mp.Fend.eval.RHOS = sqrt(XIS.^2 + ETAS.^2);

%--Storage array for throughput at each iteration
mp.thput_vec = zeros(mp.Nitr+1,1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%--Software Mask for Scoring Contrast 
%--Set Inputs
maskScore.rhoInner = mp.Fend.score.Rin; %--lambda0/D
maskScore.rhoOuter = mp.Fend.score.Rout ; %--lambda0/D
maskScore.angDeg = mp.Fend.score.ang; %--degrees
maskScore.centering = mp.centering;
maskScore.FOV = mp.Fend.FOV; %--Determines max dimension length
maskScore.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
if(isfield(mp.Fend,'shape'));  maskScore.shape = mp.Fend.shape;  end
%--Compact Model: Generate Software Mask for Scoring Contrast 
maskScore.Nxi = mp.Fend.Nxi; %--Set min dimension length to be same as for corr 
maskScore.pixresFP = mp.Fend.res;
[mp.Fend.score.mask,~,~] = falco_gen_SW_mask(maskScore); 
mp.Fend.score.settings = maskScore; %--Store values for future reference

%--Number of pixels used in the dark hole
mp.Fend.corr.Npix = sum(sum(mp.Fend.corr.mask));
mp.Fend.score.Npix = sum(sum(mp.Fend.score.mask));

%--Indices of dark hole pixels and logical masks
if(mp.flagFiber)
    mp.Fend.corr.inds = find(sum(mp.Fend.lenslet.mask,3)~=0);
    mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);
else
    mp.Fend.corr.inds = find(mp.Fend.corr.mask~=0);
    mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);
end

mp.Fend.score.inds = find(mp.Fend.score.mask~=0);
mp.Fend.score.maskBool = logical(mp.Fend.score.mask);

%% Spatial weighting of pixel intensity. 
% NOTE: For real instruments and testbeds, only the compact model should be 
% used. The full model spatial weighting is included too if in simulation 
% the full model has a different detector resolution than the compact model.

if(mp.flagFiber)
    mp.WspatialVec = ones(mp.Fend.Nlens,1);
else
    mp = falco_config_spatial_weights(mp);
    %--Extract the vector of weights at the pixel locations of the dark hole pixels.
    mp.WspatialVec = mp.Wspatial(mp.Fend.corr.inds);
end

%% Deformable Mirror (DM) 1 and 2 Parameters

if(isfield(mp,'dm1'))
    % Read the influence function header data from the FITS file
	info = fitsinfo(mp.dm1.inf_fn);
	[~, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1)); % Get index in cell array
	dx1 = info.PrimaryData.Keywords{idef, 2}; % pixel width of the influence function IN THE FILE [meters];
    [~, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
	pitch1 = info.PrimaryData.Keywords{idef, 2}; % actuator spacing x (m)
    
    mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
    mp.dm1.dx_inf0 = mp.dm1.dm_spacing*(dx1/pitch1);
    
    switch lower(mp.dm1.inf_sign(1))
        case{'-','n','m'}
            mp.dm1.inf0 = -1*mp.dm1.inf0;
        otherwise
            %--Leave coefficient as +1
    end
end    

if(isfield(mp,'dm2'))
    % Read the influence function header data from the FITS file
	info = fitsinfo(mp.dm2.inf_fn);
	[~, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1)); % Get index in cell array
	dx2 = info.PrimaryData.Keywords{idef, 2}; % pixel width of the influence function IN THE FILE [meters];
    [~, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
	pitch2 = info.PrimaryData.Keywords{idef, 2}; % actuator spacing x (m)
    
    mp.dm2.inf0 = fitsread(mp.dm2.inf_fn);
    mp.dm2.dx_inf0 = mp.dm2.dm_spacing*(dx2/pitch2);
    
    switch lower(mp.dm2.inf_sign(1))
        case{'-','n','m'}
            mp.dm2.inf0 = -1*mp.dm2.inf0;
        otherwise
            %--Leave coefficient as +1
    end
end  

%--Create influence function datacubes for each DM
if( any(mp.dm_ind==1) )
    mp.dm1.centering = mp.centering;
    mp.dm1.compact = mp.dm1;
    mp.dm1 = falco_gen_dm_poke_cube(mp.dm1, mp, mp.P2.full.dx,'NOCUBE');
    mp.dm1.compact = falco_gen_dm_poke_cube(mp.dm1.compact, mp, mp.P2.compact.dx);
else
    mp.dm1.compact = mp.dm1;
    mp.dm1 = falco_gen_dm_poke_cube(mp.dm1, mp, mp.P2.full.dx,'NOCUBE');
    mp.dm1.compact = falco_gen_dm_poke_cube(mp.dm1.compact, mp, mp.P2.compact.dx,'NOCUBE');
end

if( any(mp.dm_ind==2) ) 
    mp.dm2.centering = mp.centering;
    mp.dm2.compact = mp.dm2;
    mp.dm2.dx = mp.P2.full.dx;
    mp.dm2.compact.dx = mp.P2.compact.dx;
    
    mp.dm2 = falco_gen_dm_poke_cube(mp.dm2, mp, mp.P2.full.dx, 'NOCUBE');
    mp.dm2.compact = falco_gen_dm_poke_cube(mp.dm2.compact, mp, mp.P2.compact.dx);
else
    mp.dm2.compact = mp.dm2;
    mp.dm2.dx = mp.P2.full.dx;
    mp.dm2.compact.dx = mp.P2.compact.dx;
    
    mp.dm2 = falco_gen_dm_poke_cube(mp.dm2, mp, mp.P2.full.dx, 'NOCUBE');
    mp.dm2.compact = falco_gen_dm_poke_cube(mp.dm2.compact, mp, mp.P2.compact.dx,'NOCUBE');
end

%--Initial DM voltages
if(isfield(mp.dm1,'V')==false); mp.dm1.V = zeros(mp.dm1.Nact,mp.dm1.Nact); end
if(isfield(mp.dm2,'V')==false); mp.dm2.V = zeros(mp.dm2.Nact,mp.dm2.Nact); end

%% DM Aperture Masks (moved here because the commands mp.dm2.compact = mp.dm2; and mp.dm1.compact = mp.dm1; otherwise would overwrite the compact model masks)

if(mp.flagDM1stop)
    mp.dm1.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm1.Dstop,mp.centering);
    mp.dm1.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm1.Dstop,mp.centering);
end
if(mp.flagDM2stop)
    mp.dm2.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm2.Dstop,mp.centering);
    mp.dm2.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm2.Dstop,mp.centering);
end

%% %--First delta DM settings are zero (for covariance calculation in Kalman filters or robust controllers)
mp.dm1.dV = zeros(mp.dm1.Nact,mp.dm1.Nact); % delta voltage on DM1;
mp.dm2.dV = zeros(mp.dm2.Nact,mp.dm2.Nact); % delta voltage on DM2;
mp.dm8.dV = zeros(mp.dm8.NactTotal,1); % delta voltage on DM8;
mp.dm9.dV = zeros(mp.dm9.NactTotal,1); % delta voltage on DM9;

%% Array Sizes for Angular Spectrum Propagation with FFTs

%--Compact Model: Set nominal DM plane array sizes as a power of 2 for angular spectrum propagation with FFTs
if(any(mp.dm_ind==1) && any(mp.dm_ind==2))
    NdmPad = 2.^ceil(1 + log2(max([mp.dm1.compact.NdmPad,mp.dm2.compact.NdmPad]))); 
elseif(any(mp.dm_ind==1))
    NdmPad = 2.^ceil(1 + log2(mp.dm1.compact.NdmPad));
elseif(any(mp.dm_ind==2))
    NdmPad = 2.^ceil(1 + log2(mp.dm2.compact.NdmPad));
else
    NdmPad = 2*mp.P1.compact.Nbeam;
end
while((NdmPad < min(mp.sbp_centers)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) || (NdmPad < min(mp.sbp_centers)*abs(mp.d_P2_dm1)/mp.P2.compact.dx^2)) %--Double the zero-padding until the angular spectrum sampling requirement is not violated
    NdmPad = 2*NdmPad; 
end
mp.compact.NdmPad = NdmPad;

%--Full Model: Set nominal DM plane array sizes as a power of 2 for angular spectrum propagation with FFTs
if(any(mp.dm_ind==1) && any(mp.dm_ind==2))
    NdmPad = 2.^ceil(1 + log2(max([mp.dm1.NdmPad,mp.dm2.NdmPad]))); 
elseif(any(mp.dm_ind==1))
    NdmPad = 2.^ceil(1 + log2(mp.dm1.NdmPad));
elseif(any(mp.dm_ind==2))
    NdmPad = 2.^ceil(1 + log2(mp.dm2.NdmPad));
else
    NdmPad = 2*mp.P1.full.Nbeam;    
end
while((NdmPad < min(mp.full.lambdas)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) || (NdmPad < min(mp.full.lambdas)*abs(mp.d_P2_dm1)/mp.P2.full.dx^2)) %--Double the zero-padding until the angular spectrum sampling requirement is not violated
    NdmPad = 2*NdmPad; 
end
mp.full.NdmPad = NdmPad;

%% Initial Electric Fields for Star and Exoplanet

if(isfield(mp.P1.full,'E')==false)
    mp.P1.full.E  = ones(mp.P1.full.Narr,mp.P1.full.Narr,mp.Nwpsbp,mp.Nsbp); % Input E-field at entrance pupil
end

mp.Eplanet = mp.P1.full.E; %--Initialize the input E-field for the planet at the entrance pupil. Will apply the phase ramp later
        
if(isfield(mp.P1.compact,'E')==false)
    mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp);
end
mp.sumPupil = sum(sum(abs(mp.P1.compact.mask.*padOrCropEven(mean(mp.P1.compact.E,3),size(mp.P1.compact.mask,1) )).^2)); %--Throughput is computed with the compact model

%% Off-axis, incoherent point source (exoplanet)

if(~mp.flagFiber)
    mp.c_planet = 1; % contrast of exoplanet
    mp.x_planet = 6; % x position of exoplanet in lambda0/D
    mp.y_planet = 0; % y position of exoplanet in lambda0/D
end

%% Field Stop at Fend.(as a software mask)
mp.Fend.compact.mask = ones(mp.Fend.Neta,mp.Fend.Nxi);

%% Contrast to Normalized Intensity Map Calculation 

%% Get the starlight normalization factor for the compact and full models (to convert images to normalized intensity)
mp = falco_get_PSF_norm_factor(mp);

%--Check that the normalization of the coronagraphic PSF is correct

modvar.ttIndex = 1;
modvar.sbpIndex = mp.si_ref;
modvar.wpsbpIndex = mp.wi_ref;
modvar.whichSource = 'star';     

E0c = model_compact(mp, modvar);
I0c = abs(E0c).^2;
if(mp.flagPlot)
    figure(501); imagesc(log10(I0c)); axis xy equal tight; colorbar;
    title('(Compact Model: Normalization Check Using Starting PSF)'); 
    drawnow;
end
E0f = model_full(mp, modvar);
I0f = abs(E0f).^2;
if(mp.flagPlot)
    figure(502); imagesc(log10(I0f)); axis xy equal tight; colorbar;
    title('(Full Model: Normalization Check Using Starting PSF)'); drawnow;
end

%% Intialize delta DM voltages. Needed for Kalman filters.
%%--Save the delta from the previous command
if(any(mp.dm_ind==1));  mp.dm1.dV = 0;  end
if(any(mp.dm_ind==2));  mp.dm2.dV = 0;  end
if(any(mp.dm_ind==3));  mp.dm3.dV = 0;  end
if(any(mp.dm_ind==4));  mp.dm4.dV = 0;  end
if(any(mp.dm_ind==5));  mp.dm5.dV = 0;  end
if(any(mp.dm_ind==6));  mp.dm6.dV = 0;  end
if(any(mp.dm_ind==7));  mp.dm7.dV = 0;  end
if(any(mp.dm_ind==8));  mp.dm8.dV = 0;  end
if(any(mp.dm_ind==9));  mp.dm9.dV = 0;  end

%% Intialize tied actuator pairs if not already defined. 
% Dimensions of the pair list is [Npairs x 2]
%%--Save the delta from the previous command
if(any(mp.dm_ind==1)); if(isfield(mp.dm1,'tied')==false); mp.dm1.tied = []; end; end
if(any(mp.dm_ind==2)); if(isfield(mp.dm2,'tied')==false); mp.dm2.tied = []; end; end
if(any(mp.dm_ind==3)); if(isfield(mp.dm3,'tied')==false); mp.dm3.tied = []; end; end
if(any(mp.dm_ind==4)); if(isfield(mp.dm4,'tied')==false); mp.dm4.tied = []; end; end
if(any(mp.dm_ind==5)); if(isfield(mp.dm5,'tied')==false); mp.dm5.tied = []; end; end
if(any(mp.dm_ind==6)); if(isfield(mp.dm6,'tied')==false); mp.dm6.tied = []; end; end
if(any(mp.dm_ind==7)); if(isfield(mp.dm7,'tied')==false); mp.dm7.tied = []; end; end
if(any(mp.dm_ind==8)); if(isfield(mp.dm8,'tied')==false); mp.dm8.tied = []; end; end
if(any(mp.dm_ind==9)); if(isfield(mp.dm9,'tied')==false); mp.dm9.tied = []; end; end

%% Storage Arrays for DM Metrics
%--EFC regularization history
out.log10regHist = zeros(mp.Nitr,1);

%--Peak-to-Valley DM voltages
out.dm1.Vpv = zeros(mp.Nitr,1);
out.dm2.Vpv = zeros(mp.Nitr,1);
out.dm8.Vpv = zeros(mp.Nitr,1);
out.dm9.Vpv = zeros(mp.Nitr,1);

%--Peak-to-Valley DM surfaces
out.dm1.Spv = zeros(mp.Nitr,1);
out.dm2.Spv = zeros(mp.Nitr,1);
out.dm8.Spv = zeros(mp.Nitr,1);
out.dm9.Spv = zeros(mp.Nitr,1);

%--RMS DM surfaces
out.dm1.Srms = zeros(mp.Nitr,1);
out.dm2.Srms = zeros(mp.Nitr,1);
out.dm8.Srms = zeros(mp.Nitr,1);
out.dm9.Srms = zeros(mp.Nitr,1);

%--Zernike sensitivities to 1nm RMS
if(isfield(mp.eval,'Rsens')==false);  mp.eval.Rsens = [];   end
if(isfield(mp.eval,'indsZnoll')==false);  mp.eval.indsZnoll = [2,3];   end
Nannuli = size(mp.eval.Rsens,1);
Nzern = length(mp.eval.indsZnoll);
out.Zsens = zeros(Nzern,Nannuli,mp.Nitr);

%--Store the DM commands at each iteration
if(isfield(mp,'dm1')); if(isfield(mp.dm1,'V'));  out.dm1.Vall = zeros(mp.dm1.Nact,mp.dm1.Nact,mp.Nitr+1);  end; end
if(isfield(mp,'dm2')); if(isfield(mp.dm2,'V'));  out.dm2.Vall = zeros(mp.dm2.Nact,mp.dm2.Nact,mp.Nitr+1); end; end
if(isfield(mp,'dm8')); if(isfield(mp.dm8,'V'));  out.dm8.Vall = zeros(mp.dm8.NactTotal,mp.Nitr+1); end; end
if(isfield(mp,'dm9')); if(isfield(mp.dm9,'V'));  out.dm9.Vall = zeros(mp.dm9.NactTotal,mp.Nitr+1); end; end

%% 
fprintf('\nBeginning Trial %d of Series %d.\n',mp.TrialNum,mp.SeriesNum);

end %--END OF FUNCTION