% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

function mp = falco_set_optional_variables(mp)

%% Intializations of structures (if they don't exist yet)
mp.jac.dummy = 1;
mp.est.dummy = 1;
mp.compact.dummy = 1;
mp.full.dummy = 1;
mp.dm1.dummy = 1;
mp.dm2.dummy = 1;
mp.Fend.eval.dummy = 1;
mp.path.dummy = 1;

%% Default File Paths for Data Storage (all excluded from git)

%--Get the falco path for making the other default paths
[filepath, name, ext] = fileparts(mfilename('fullpath'));
mp.path.falco = filepath(1:end-5); % remove "setup" from the end of the path

%--Store minimal data to re-construct the data from the run: the config files and "out" structure after a trial go here
if(isfield(mp.path,'config')==false);  mp.path.config = [mp.path.falco filesep 'data' filesep 'brief' filesep];  end

%--Entire final workspace from FALCO gets saved here.
if(isfield(mp.path,'ws')==false);  mp.path.ws = [mp.path.falco filesep 'data' filesep 'ws' filesep];  end

if(isfield(mp.path,'ws')==false); mp.path.ws = [mp.path.falco 'data' filesep 'ws' filesep]; end % Store final workspace data here
if(isfield(mp.path,'maps')==false); mp.path.falcoaps = [mp.path.falco 'maps' filesep]; end % Maps go here
if(isfield(mp.path,'jac')==false); mp.path.jac = [mp.path.falco 'data' filesep 'jac' filesep]; end % Store the control Jacobians here
if(isfield(mp.path,'images')==false); mp.path.images = [mp.path.falco 'data' filesep 'images' filesep]; end % Store all full, reduced images here
if(isfield(mp.path,'dm')==false); mp.path.dm = [mp.path.falco 'data' filesep 'DM' filesep]; end % Store DM command maps here
if(isfield(mp.path,'wsInProgress')==false); mp.path.wsInProgress = [mp.path.falco 'data' filesep 'wsInProgress' filesep]; end % Store in progress workspace data here

%% Optional/hidden boolean flags
%--Saving data
if(isfield(mp,'flagSaveWS')==false);  mp.flagSaveWS = false;  end  %--Whether to save out the entire workspace at the end of the trial. Can take up lots of space.
if(isfield(mp,'flagSaveEachItr')==false);  mp.flagSaveEachItr = false;  end  %--Whether to save out the performance at each iteration. Useful for long trials in case it crashes or is stopped early.
if(isfield(mp,'flagSVD')==false);  mp.flagSVD = false;  end    %--Whether to compute and save the singular mode spectrum of the control Jacobian (each iteration)
%--Jacobian or controller related
if(isfield(mp,'flagTrainModel')==false);  mp.flagTrainModel = false;  end  %--Whether to call the Expectation-Maximization (E-M) algorithm to improve the linearized model. 
if(isfield(mp,'flagUseLearnedJac')==false);  mp.flagUseLearnedJac = false;  end  %--Whether to load and use an improved Jacobian from the Expectation-Maximization (E-M) algorithm 
if(isfield(mp.est,'flagUseJac')==false);  mp.est.flagUseJac = false;  end   %--Whether to use the Jacobian or not for estimation. (If not using Jacobian, model is called and differenced.)
if(isfield(mp.ctrl,'flagUseModel')==false);  mp.ctrl.flagUseModel = false;  end %--Whether to perform a model-based (instead of empirical) grid search for the controller

%--Model options (Very specialized cases--not for the average user)
if(isfield(mp,'flagFiber')==false);  mp.flagFiber = false;  end  %--Whether to couple the final image through lenslets and a single mode fiber.
if(isfield(mp,'flagLenslet')==false); mp.flagLenslet = false; end %--Flag to propagate through a lenslet array placed in Fend before coupling light into fibers
if(isfield(mp,'flagDMwfe')==false);  mp.flagDMwfe = false;  end  %--Temporary for BMC quilting study. Adds print-through to the DM surface.
if(isfield(mp,'flagWFS')==false);  mp.flagWFS = false;  end  %--Whether to activate the WFS mode 

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

%% Optional/hidden variables
if(isfield(mp.full,'pol_conds')==false);  mp.full.pol_conds = 0;  end %--Vector of which polarization state(s) to use when creating images from the full model. Currently only used with PROPER full models from John Krist.
if(isfield(mp,'propMethodPTP')==false);  mp.propMethodPTP = 'fft';  end %--Propagation method for postage stamps around the influence functions. 'mft' or 'fft'
if(isfield(mp,'apodType')==false);  mp.apodType = 'none';  end %--Type of apodizer. Only use this variable when generating the apodizer. Currently only binary-ring or grayscale apodizers can be generated.

%--Sensitivities to Zernike-Mode Perturbations
if(isfield(mp.full,'ZrmsVal')==false);  mp.full.ZrmsVal = 1e-9;  end %--Amount of RMS Zernike mode used to calculate aberration sensitivities [meters]. WFIRST CGI uses 1e-9, and LUVOIR and HabEx use 1e-10. 
if(isfield(mp.eval,'Rsens')==false);  mp.eval.Rsens = [];   end
if(isfield(mp.eval,'indsZnoll')==false);  mp.eval.indsZnoll = [2,3];   end

%--Deformable mirror settings
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

%--Off-axis, incoherent point source (exoplanet). Used if modvar.whichSource = 'exoplanet'
if(isfield(mp,'c_planet')==false);  mp.c_planet = 1e-10;  end % flux ratio of of exoplanet to star
if(isfield(mp,'x_planet')==false);  mp.x_planet = 5;  end % xi position of exoplanet in lambda0/D
if(isfield(mp,'y_planet')==false);  mp.y_planet = 1;  end % eta position of exoplanet in lambda0/D

%--Control
if(isfield(mp.jac,'zerns')==false); mp.jac.zerns = 1; else; mp.jac.Zcoef = 1; end %--Zernike modes in Jacobian
if(isfield(mp,'WspatialDef')==false);  mp.WspatialDef = [];  end %--spatial weights for the Jacobian

%--Performance Evaluation
if(isfield(mp.Fend.eval,'res')==false);  mp.Fend.eval.res = 10;  end % pixels per lambda0/D in compact evaluation model's final focus
mp.mas2lam0D = 1/(mp.lambda0/mp.P1.D*180/pi*3600*1000); %% Conversion factor: milliarcseconds (mas) to lambda0/D

%--Training Data: mp.NitrTrain = 5;  %--The number of correction iterations to use per round of training data for the adaptive Jacobian (E-M) algorithm.
%--Zernike sensitivities to 1nm RMS: which noll indices in which annuli, given by mp.eval.indsZnoll and mp.eval.Rsens 
%--Tied actuator pair definitions: See Section with variables mp.dmX.tied for X=1:9
%--Quantization of DM actuation steps based on least significant bit of the
% DAC (digital-analog converter). In height, so called HminStep. If HminStep (minimum step in H) is defined, then quantize the DM voltages
% Variables to define if wanted: mp.dm1.HminStep, mp.dm2.HminStep

%% Initialize some basic attributes for all DMs (which include hybrid FPMs).
mp.dm1.NactTotal=0; mp.dm2.NactTotal=0; mp.dm3.NactTotal=0; mp.dm4.NactTotal=0; mp.dm5.NactTotal=0; mp.dm6.NactTotal=0; mp.dm7.NactTotal=0; mp.dm8.NactTotal=0; mp.dm9.NactTotal=0; 
mp.dm1.Nele=0; mp.dm2.Nele=0; mp.dm3.Nele=0; mp.dm4.Nele=0; mp.dm5.Nele=0; mp.dm6.Nele=0; mp.dm7.Nele=0; mp.dm8.Nele=0; mp.dm9.Nele=0; %--Initialize for Jacobian calculations later. 

%--Intialize delta DM voltages. Needed for Kalman filters.
if(any(mp.dm_ind==1));  mp.dm1.dV = 0;  end
if(any(mp.dm_ind==2));  mp.dm2.dV = 0;  end
if(any(mp.dm_ind==3));  mp.dm3.dV = 0;  end
if(any(mp.dm_ind==4));  mp.dm4.dV = 0;  end
if(any(mp.dm_ind==5));  mp.dm5.dV = 0;  end
if(any(mp.dm_ind==6));  mp.dm6.dV = 0;  end
if(any(mp.dm_ind==7));  mp.dm7.dV = 0;  end
if(any(mp.dm_ind==8));  mp.dm8.dV = 0;  end
if(any(mp.dm_ind==9));  mp.dm9.dV = 0;  end

%%--Intialize tied actuator pairs if not already defined. 
% Dimensions of the pair list is [Npairs x 2]
if(any(mp.dm_ind==3)); if(isfield(mp.dm3,'tied')==false); mp.dm3.tied = zeros(0,2); end; end
if(any(mp.dm_ind==4)); if(isfield(mp.dm4,'tied')==false); mp.dm4.tied = zeros(0,2); end; end
if(any(mp.dm_ind==5)); if(isfield(mp.dm5,'tied')==false); mp.dm5.tied = zeros(0,2); end; end
if(any(mp.dm_ind==6)); if(isfield(mp.dm6,'tied')==false); mp.dm6.tied = zeros(0,2); end; end
if(any(mp.dm_ind==7)); if(isfield(mp.dm7,'tied')==false); mp.dm7.tied = zeros(0,2); end; end
if(any(mp.dm_ind==8)); if(isfield(mp.dm8,'tied')==false); mp.dm8.tied = zeros(0,2); end; end
if(any(mp.dm_ind==9)); if(isfield(mp.dm9,'tied')==false); mp.dm9.tied = zeros(0,2); end; end

end