% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to finish initializing the workspace prior to wavefront
% estimation and control.
%
% Created by A.J. Riggs on 2017-10-31.

function [mp,out] = falco_init_ws(fn_config)

%% Read inputs as structures from a .mat config file

load(fn_config,'mp');

mainPath = mp.path.falco;

disp(['DM 1-to-2 Fresnel number (using radius) = ',num2str((mp.P2.D/2)^2/(mp.d_dm1_dm2*mp.lambda0))]);

%% Intializations of structures (if they don't exist yet)
mp.jac.dummy = 1;
mp.est.dummy = 1;
mp.compact.dummy = 1;
mp.full.dummy = 1;
mp.dm1.dummy = 1;
mp.dm2.dummy = 1;

%% Optional/Hidden flags
%--Saving data
if(isfield(mp,'flagSaveWS')==false);  mp.flagSaveWS = false;  end  %--Whehter to save otu the entire workspace at the end of the trial. Can take up lots of space.
if(isfield(mp,'flagSaveEachItr')==false);  mp.flagSaveEachItr = false;  end  %--Whether to save out the performance at each iteration. Useful for long trials in case it crashes or is stopped early.
if(isfield(mp,'flagSVD')==false);  mp.flagSVD = false;  end    %--Whether to compute and save the singular mode spectrum of the control Jacobian (each iteration)
%--Jacobian or controller related
if(isfield(mp,'flagTrainModel')==false);  mp.flagTrainModel = false;  end  %--Whether to call the Expectation-Maximization (E-M) algorithm to improve the linearized model. 
if(isfield(mp,'flagUseLearnedJac')==false);  mp.flagUseLearnedJac = false;  end  %--Whether to load and use an improved Jacobian from the Expectation-Maximization (E-M) algorithm 
if(isfield(mp.est,'flagUseJac')==false);  mp.est.flagUseJac = false;  end   %--Whether to use the Jacobian or not for estimation. (If not using Jacobian, model is called and differenced.)
if(isfield(mp.ctrl,'flagUseModel')==false);  mp.ctrl.flagUseModel = false;  end %--Whether to perform a model-based (instead of empirical) grid search for the controller
%--Deformable mirror actuator constraints or bounds
if(isfield(mp.dm1,'Vmin')==false);  mp.dm1.Vmin = -1000;  end %--Min allowed voltage command
if(isfield(mp.dm1,'Vmax')==false);  mp.dm1.Vmax = 1000;  end %--Max allowed voltage command
if(isfield(mp.dm1,'pinned')==false);  mp.dm1.pinned = [];  end %--Indices of pinned actuators
if(isfield(mp.dm1,'Vpinned')==false);  mp.dm1.Vpinned = [];  end %--(Fixed) voltage commands of pinned actuators
if(isfield(mp.dm1,'tied')==false);  mp.dm1.tied = zeros(0,2);  end %--Indices of paired actuators. Two indices per row
if(isfield(mp.dm1,'flagNbrRule')==false);  mp.dm1.flagNbrRule = false;  end %--Whether to set constraints on neighboring actuator voltage differences. If set to true, need to define mp.dm1.dVnbr
if(isfield(mp.dm2,'Vmin')==false);  mp.dm2.Vmin = -1000;  end %--Min allowed voltage command
if(isfield(mp.dm2,'Vmax')==false);  mp.dm2.Vmax = 1000;  end %--Max allowed voltage command
if(isfield(mp.dm2,'pinned')==false);  mp.dm2.pinned = [];  end %--Indices of pinned actuators
if(isfield(mp.dm2,'Vpinned')==false);  mp.dm2.Vpinned = [];  end %--(Fixed) voltage commands of pinned actuators
if(isfield(mp.dm2,'tied')==false);  mp.dm2.tied = zeros(0,2);  end %--Indices of paired actuators. Two indices per row
if(isfield(mp.dm2,'flagNbrRule')==false);  mp.dm2.flagNbrRule = false;  end %--Whether to set constraints on neighboring actuator voltage differences. If set to true, need to define mp.dm1.dVnbr

%--Model options (Very specialized cases--not for the average user)
if(isfield(mp,'flagFiber')==false);  mp.flagFiber = false;  end  %--Whether to couple the final image through lenslets and a single mode fiber.
if(isfield(mp,'flagLenslet')==false); mp.flagLenslet = false; end %--Flag to propagate through a lenslet array placed in Fend before coupling light into fibers
if(isfield(mp,'flagDMwfe')==false);  mp.flagDMwfe = false;  end  %--Temporary for BMC quilting study. Adds print-through to the DM surface.

%--Whether to generate or load various masks: compact model
if(isfield(mp.compact,'flagGenPupil')==false);  mp.compact.flagGenPupil = true;  end
if(isfield(mp.compact,'flagGenApod')==false);  mp.compact.flagGenApod = false;  end %--Different! Apodizer generation defaults to false.
if(isfield(mp.compact,'flagGenFPM')==false);  mp.compact.flagGenFPM = true;  end
if(isfield(mp.compact,'flagGenLS')==false);  mp.compact.flagGenLS = true;  end
%--Whether to generate or load various masks: full model
if(isfield(mp.full,'flagPROPER')==false);  mp.full.flagPROPER = false;  end %--Whether to use a full model written in PROPER. If true, then load (don't generate) all masks for the full model
if(mp.full.flagPROPER)
    mp.full.flagGenPupil = false;
    mp.full.flagGenApod = false;
    mp.full.flagGenFPM = false;
    mp.full.flagGenLS = false;
end
if(isfield(mp.full,'flagGenPupil')==false);  mp.full.flagGenPupil = true;  end
if(isfield(mp.full,'flagGenApod')==false);  mp.full.flagGenApod = false;  end %--Different! Apodizer generation defaults to false.
if(isfield(mp.full,'flagGenFPM')==false);  mp.full.flagGenFPM = true;  end
if(isfield(mp.full,'flagGenLS')==false);  mp.full.flagGenLS = true;  end

%% Optional/Hidden variables
if(isfield(mp.full,'ZrmsVal')==false);  mp.full.ZrmsVal = 1e-9;  end %--Amount of RMS Zernike mode used to calculate aberration sensitivities [meters]. WFIRST CGI uses 1e-9, and LUVOIR and HabEx use 1e-10. 
if(isfield(mp.full,'pol_conds')==false);  mp.full.pol_conds = 0;  end %--Vector of which polarization state(s) to use when creating images from the full model. Currently only used with PROPER full models from John Krist.
if(isfield(mp,'propMethodPTP')==false);  mp.propMethodPTP = 'fft';  end %--Propagation method for postage stamps around the influence functions. 'mft' or 'fft'
if(isfield(mp,'apodType')==false);  mp.apodType = 'none';  end %--Type of apodizer. Only use this variable when generating the apodizer. Currently only binary-ring or grayscale apodizers can be generated.

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

%% Bandwidth and Wavelength Specs: Compact Model

if(isfield(mp,'Nwpsbp')==false);  mp.Nwpsbp = 1;  end

%--Center(-ish) wavelength indices (ref = reference). (Only the center if
%  an odd number of wavelengths is used.)
mp.si_ref = ceil(mp.Nsbp/2);

%--Wavelengths used for Compact Model (and Jacobian Model)
mp.sbp_weights = ones(mp.Nsbp,1);
if(mp.Nwpsbp==1 && mp.flagSim) %--Set ctrl wavelengths evenly between endpoints (inclusive) of the total bandpass.
    if(mp.Nsbp==1)
        mp.sbp_centers = mp.lambda0;
    else
        mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
        mp.sbp_weights(1) = 1/2; %--Give end sub-bands half weighting
        mp.sbp_weights(end) = 1/2; %--Give end sub-bands half weighting
    end
else %--For cases with multiple sub-bands: Choose wavelengths to be at subbandpass centers since the wavelength samples will span to the full extent of the sub-bands.
    mp.fracBWcent2cent = mp.fracBW*(1-1/mp.Nsbp); %--Bandwidth between centers of endpoint subbandpasses.
    mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBWcent2cent/2,1+mp.fracBWcent2cent/2,mp.Nsbp); %--Space evenly at the centers of the subbandpasses.
end
mp.sbp_weights = mp.sbp_weights/sum(mp.sbp_weights); %--Normalize the sum of the weights

fprintf(' Using %d discrete wavelength(s) in each of %d sub-bandpasses over a %.1f%% total bandpass \n', mp.Nwpsbp, mp.Nsbp,100*mp.fracBW);
fprintf('Sub-bandpasses are centered at wavelengths [nm]:\t '); fprintf('%.2f  ',1e9*mp.sbp_centers); fprintf('\n\n');

%% Bandwidth and Wavelength Specs: Full Model

%--Center(-ish) wavelength indices (ref = reference). (Only the center if an odd number of wavelengths is used.)
mp.wi_ref = ceil(mp.Nwpsbp/2);

%--Wavelength factors/weights within each sub-bandpass. For full model only
mp.full.lambda_weights = ones(mp.Nwpsbp,1); %--Initialize as all ones. Weights within a single sub-bandpass
if(mp.Nwpsbp==1) %--Give equal weighting to all wavelengths
    mp.full.dlam = 0; %--Delta lambda between every wavelength in the sub-band in the full model
else %--Give half weighting to the end wavelengths
    %--Spectral weighting in image
    mp.full.lambda_weights(1) = 1/2; %--Give end wavelengths half weighting
    mp.full.lambda_weights(end) = 1/2; %--Give end wavelengths half weighting
    mp.fracBWsbp = mp.fracBW/mp.Nsbp; %--Bandwidth per sub-bandpass
    %--Indexing of wavelengths in each sub-bandpass
    sbp_facs = linspace(1-mp.fracBWsbp/2,1+mp.fracBWsbp/2,mp.Nwpsbp); %--Factor applied to lambda0 only
    mp.full.dlam = (sbp_facs(2) - sbp_facs(1))*mp.lambda0; %--Delta lambda between every wavelength in the full model 
end
mp.full.lambda_weights = mp.full.lambda_weights/sum(mp.full.lambda_weights); %--Normalize sum of the weights (within the sub-bandpass)

%--Make vector of all wavelengths and weights used in the full model
lambdas = zeros(mp.Nsbp*mp.Nwpsbp,1);
lambda_weights_all = zeros(mp.Nsbp*mp.Nwpsbp,1);
mp.full.lambdasMat = zeros(mp.Nsbp,mp.Nwpsbp);
mp.full.indsLambdaMat = zeros(mp.Nsbp*mp.Nwpsbp,2);
counter = 0;
for si=1:mp.Nsbp
    mp.full.lambdasMat(si,:) = (-(mp.Nwpsbp-1)/2:(mp.Nwpsbp-1)/2)*mp.full.dlam + mp.sbp_centers(si); 
    for wi=1:mp.Nwpsbp
        counter = counter+1;
        lambdas(counter) = mp.full.lambdasMat(si,wi);
        lambda_weights_all(counter) = mp.sbp_weights(si)*mp.full.lambda_weights(wi);
        mp.full.indsLambdaMat(counter,:) = [si,wi];
    end
end

%--Get rid of redundant wavelengths in the complete list, and sum weights for repeated wavelengths
% indices of unique wavelengths
[~, inds_unique] = unique(round(1e12*lambdas)); %--Check equality at the picometer level for wavelength
mp.full.indsLambdaUnique = inds_unique;
% indices of duplicate wavelengths
duplicate_inds = setdiff( 1:length(lambdas) , inds_unique);
% duplicate weight values
duplicate_values = lambda_weights_all(duplicate_inds);

%--Shorten the vectors to contain only unique values. Combine weights for repeated wavelengths.
mp.full.lambdas = lambdas(inds_unique);
mp.full.lambda_weights_all = lambda_weights_all(inds_unique);
for idup=1:length(duplicate_inds)
    lambda = lambdas(duplicate_inds(idup));
    weight = lambda_weights_all(duplicate_inds(idup));
    ind = find(abs(mp.full.lambdas-lambda)<=1e-11);
    mp.full.lambda_weights_all(ind) = mp.full.lambda_weights_all(ind) + weight;
end
mp.full.NlamUnique = length(inds_unique);

%%
% %--Make vector of all wavelengths and weights used in the full model
% mp.full.lambdas = zeros(mp.Nsbp*mp.Nwpsbp,1);
% mp.full.lambda_weights_all = zeros(mp.Nsbp*mp.Nwpsbp,1);
% mp.full.lambdaInds = zeros(mp.Nsbp*mp.Nwpsbp,2); %--Indices as a guide
% counter = 0;
% for si=1:mp.Nsbp
%     for wi=1:mp.Nwpsbp
%         counter = counter+1;
%         mp.full.lambdas(counter) = mp.sbp_centers(si)*mp.full.sbp_facs(wi);
%         mp.full.lambda_weights_all = mp.sbp_weights(si)*mp.full.lambda_weights(wi);
%         mp.full.lambdaInds(counter,:) = [si,wi]; %--Indices
%     end
% end
% mp.full.Nlambdas = mp.Nsbp*mp.Nwpsbp;

%% Zernike and Chromatic Weighting of the Control Jacobian
if(isfield(mp.jac,'zerns')==false);  mp.jac.zerns = 1; end %--Which Zernike modes to include in Jacobian [Noll index]. Always include 1 for piston term.
if(isfield(mp.jac,'Zcoef')==false);  mp.jac.Zcoef = 1e-9*ones(10,1); end %--meters RMS of Zernike aberrations. (piston value is reset to 1 later for correct normalization)

mp = falco_config_jac_weights(mp);

%% Pupil Masks
mp = falco_config_gen_chosen_pupil(mp); %--input pupil mask
mp = falco_config_gen_chosen_apodizer(mp); %--apodizer mask
mp = falco_config_gen_chosen_LS(mp); %--Lyot stop

%--Compare apodizer and pupil mask overlap
if(mp.flagApod)
    if(mp.flagPlot)
        figure(600); imagesc(mp.P3.compact.mask - mp.P1.compact.mask,[-1 1]); axis xy equal tight; colorbar; drawnow;
    end
end

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
    case {'LC','APLC'} %--Occulting spot FPM (can be partially transmissive)
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

if(mp.flagLenslet)
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
if(isfield(mp.Fend,'clockAngDeg'));  maskCorr.clockAngDeg = mp.Fend.clockAngDeg;  end

%--Compact Model: Generate Software Mask for Correction 
[mp.Fend.corr.mask, mp.Fend.xisDL, mp.Fend.etasDL] = falco_gen_SW_mask(maskCorr); 
mp.Fend.corr.settings = maskCorr; %--Store values for future reference
%--Size of the output image 
%--Need the sizes to be the same for the correction and scoring masks
mp.Fend.Nxi  = size(mp.Fend.corr.mask,2);
mp.Fend.Neta = size(mp.Fend.corr.mask,1);

[XIS,ETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
mp.Fend.RHOS = sqrt(XIS.^2 + ETAS.^2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(mp.flagFiber)
    
    V(1,1,:) = 2*pi*mp.fiber.a_phys*mp.fiber.NA./mp.sbp_centers;
    W = 1.1428*V - 0.996;
    U = sqrt(V.^2 - W.^2);
    
    maskFiberCore.rhoInner = 0;
    maskFiberCore.rhoOuter = mp.fiber.a;
    maskFiberCore.angDeg = 180;
    maskFiberCore.whichSide = mp.Fend.sides;
    
    maskFiberCladding.rhoInner = mp.fiber.a;
    maskFiberCladding.rhoOuter = 20;
    maskFiberCladding.angDeg = 180;
    maskFiberCladding.whichSide = mp.Fend.sides;
    
    if(mp.flagLenslet)
        %--Mask defining the area covered by the lenslet.  Only the immediate area
        %around the lenslet is propagated, saving computation time.  This lenslet
        %can then be moved around to different positions in Fend.
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
        
        maskFiberCore.pixresFP = mp.F5.res;
        maskFiberCore.FOV = mp.F5.FOV;
        [mp.F5.fiberCore.mask, ~, ~] = falco_gen_SW_mask(maskFiberCore);
        
        maskFiberCladding.pixresFP = mp.F5.res;
        maskFiberCladding.FOV = mp.F5.FOV;
        [mp.F5.fiberCladding.mask, ~, ~] = falco_gen_SW_mask(maskFiberCladding);
        
        [F5XIS, F5ETAS] = meshgrid(mp.F5.xisDL, mp.F5.etasDL);
        
        mp.F5.RHOS = sqrt((F5XIS - mp.F5.fiberPos(1)).^2 + (F5ETAS - mp.F5.fiberPos(2)).^2);
        mp.F5.fiberCore.mode = mp.F5.fiberCore.mask.*besselj(0, U.*mp.F5.RHOS/mp.fiber.a)./besselj(0,U);
        mp.F5.fiberCladding.mode = mp.F5.fiberCladding.mask.*besselk(0, W.*mp.F5.RHOS/mp.fiber.a)./besselk(0,W);
        mp.F5.fiberCladding.mode(isnan(mp.F5.fiberCladding.mode)) = 0;
        mp.F5.fiberMode = mp.F5.fiberCore.mode + mp.F5.fiberCladding.mode;
        fiberModeNorm = sqrt(sum(sum(abs(mp.F5.fiberMode).^2)));
        mp.F5.fiberMode = mp.F5.fiberMode./fiberModeNorm;
    else
                
        for nfib = 1:mp.Fend.Nfiber    
            maskFiberCore.pixresFP = mp.Fend.res;
            maskFiberCore.FOV = mp.Fend.FOV;
            maskFiberCore.xi_cen = mp.Fend.x_fiber(nfib);
            maskFiberCore.eta_cen = mp.Fend.y_fiber(nfib);
            [mp.Fend.fiberCore.mask, ~, ~] = falco_gen_SW_mask(maskFiberCore);

            maskFiberCladding.pixresFP = mp.Fend.res;
            maskFiberCladding.FOV = mp.Fend.FOV;
            maskFiberCladding.xi_cen = mp.Fend.x_fiber(nfib);
            maskFiberCladding.eta_cen = mp.Fend.y_fiber(nfib);
            [mp.Fend.fiberCladding.mask, ~, ~] = falco_gen_SW_mask(maskFiberCladding);

            [FENDXIS, FENDETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
            
            FENDRHOS = sqrt((FENDXIS - mp.Fend.x_fiber(nfib)).^2 + (FENDETAS - mp.Fend.y_fiber(nfib)).^2);
            mp.Fend.fiberCore.mode = mp.Fend.fiberCore.mask.*besselj(0, U.*FENDRHOS/mp.fiber.a)./besselj(0,U);
            mp.Fend.fiberCladding.mode = mp.Fend.fiberCladding.mask.*besselk(0, W.*FENDRHOS/mp.fiber.a)./besselk(0,W);
            mp.Fend.fiberCladding.mode(isnan(mp.Fend.fiberCladding.mode)) = 0;
            fiberModeTemp = mp.Fend.fiberCore.mode + mp.Fend.fiberCladding.mode;
            fiberModeNorm = sqrt(sum(sum(abs(fiberModeTemp).^2)));
            mp.Fend.fiberMode(:,:,:,nfib) = fiberModeTemp./fiberModeNorm;
        end
        
    end
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
if(isfield(mp.Fend,'clockAngDeg'));  maskScore.clockAngDeg = mp.Fend.clockAngDeg;  end

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
    if(mp.flagLenslet)
        mp.Fend.corr.inds = find(sum(mp.Fend.lenslet.mask,3)~=0);
        mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);
    else
        mp.Fend.corr.inds = find(mp.Fend.corr.mask~=0);
        mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);
    end
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
    if(mp.flagLenslet)
        mp.WspatialVec = ones(mp.Fend.Nlens,1);
    else
        mp = falco_config_spatial_weights(mp);
        %--Extract the vector of weights at the pixel locations of the dark hole pixels.
        mp.WspatialVec = mp.Wspatial(mp.Fend.corr.inds);
    end
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