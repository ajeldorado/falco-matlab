% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

function mp = falco_set_optional_variables(mp)

%% Intializations of structures (if they don't exist yet)
mp.jac.dummy = 1;
mp.est.dummy = 1;
mp.compact.dummy = 1;
mp.full.dummy = 1;
mp.dm1.dummy = 1;
mp.dm2.dummy = 1;

%% Optional/hidden boolean flags
%--Saving data
if(isfield(mp,'flagSaveWS')==false);  mp.flagSaveWS = false;  end  %--Whehter to save otu the entire workspace at the end of the trial. Can take up lots of space.
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

%--Sensitivities Zernike-Mode Perturbations
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

%--Training Data: mp.NitrTrain = 5;  %--The number of correction iterations to use per round of training data for the adaptive Jacobian (E-M) algorithm.
%--Zernike sensitivities to 1nm RMS: which noll indices in which annuli, given by mp.eval.indsZnoll and mp.eval.Rsens 
%--Tied actuator pair definitions: See Section with variables mp.dmX.tied for X=1:9
%--Quantization of DM actuation steps based on least significant bit of the
% DAC (digital-analog converter). In height, so called HminStep. If HminStep (minimum step in H) is defined, then quantize the DM voltages
% Variables to define if wanted: mp.dm1.HminStep, mp.dm2.HminStep

end