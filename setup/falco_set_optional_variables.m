% Copyright 2018-2021 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Define necessary, lesser-used variables if not already defined.
% This function exists to enable backwards compatibility when adding new
% features.
%
% INPUTS
% mp : structure of model parameters
%
% OUTPUTS
% mp : structure of model parameters

function mp = falco_set_optional_variables(mp)

%% Intializations of structures (if they don't exist yet)
mp.jac.dummy = 1;
mp.est.dummy = 1;
mp.star.dummy = 1;
mp.compact.star.dummy = 1;
mp.jac.star.dummy = 1;
mp.full.dummy = 1;
mp.dm1.dummy = 1;
mp.dm2.dummy = 1;
mp.Fend.eval.dummy = 1;
mp.path.dummy = 1;
mp.detector.dummy = 1;
mp.scc.dummy = 1;
mp.iefc.dummy = 1;

%% Default File Paths for Data Storage (all excluded from git)

%--Get the falco path for making the other default paths
mp.path.falco = fileparts(fileparts(mfilename('fullpath')));

%--Store minimal data to re-construct the data from the run: the config files and "out" structure after a trial go here
if(isfield(mp.path,'config')==false);  mp.path.config = [mp.path.falco filesep 'data' filesep 'brief' filesep];  end

%--Entire final workspace from FALCO gets saved here.
if(isfield(mp.path,'ws')==false);  mp.path.ws = [mp.path.falco filesep 'data' filesep 'ws' filesep];  end % Store final workspace data here
if(isfield(mp.path,'maps')==false); mp.path.maps = [mp.path.falco filesep 'maps' filesep]; end % Maps go here
if(isfield(mp.path,'jac')==false); mp.path.jac = [mp.path.falco filesep 'data' filesep 'jac' filesep]; end % Store the control Jacobians here
if(isfield(mp.path,'images')==false); mp.path.images = [mp.path.falco filesep 'data' filesep 'images' filesep]; end % Store all full, reduced images here
if(isfield(mp.path,'dm')==false); mp.path.dm = [mp.path.falco filesep 'data' filesep 'DM' filesep]; end % Store DM command maps here

%% Optional/hidden boolean flags

%--Debugging mode
if ~isfield(mp, 'debug'); mp.debug = false; end

%--Saving data
if(isfield(mp,'flagSaveWS')==false);  mp.flagSaveWS = false;  end  %--Whether to save out the entire workspace at the end of the trial. Can take up lots of space.
if(isfield(mp,'flagSVD')==false);  mp.flagSVD = false;  end    %--Whether to compute and save the singular mode spectrum of the control Jacobian (each iteration)
if(isfield(mp.jac,'fn')==false);  mp.jac.fn = 'default_jac_name.mat';  end  % Name of the file to save the Jacobian into. Do not include the path.

%--Optical model/layout
if(isfield(mp.full,'flagPROPER')==false);  mp.full.flagPROPER = false;  end %--Whether to use a full model written in PROPER. If true, then load (don't generate) all masks for the full model
if(isfield(mp,'flagRotation')==false);  mp.flagRotation = true;  end %--Whether to have the E-field rotate 180 degrees from one pupil to the next. Does not apply to PROPER full models.

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

%--Whether to use an apodizer at all
if(isfield(mp,'flagApod')==false);  mp.flagApod = false;  end

%--Lyot stop symmetry (for WFIRST/Roman only)
if(isfield(mp.P4,'flagSymm')==false);  mp.P4.flagSymm = false;  end

%% Optional/hidden variables

if ~isfield(mp, 'runLabel');  mp.runLabel = 'default_label_';  end

% How many stars to use and their positions
% mp.star is for the full model, and mp.compact.star is for the compact and
% Jacobian models.
if ~isfield(mp.star, 'count');  mp.star.count = 1;  end
if ~isfield(mp.star, 'xiOffsetVec');  mp.star.xiOffsetVec = 0;  end
if ~isfield(mp.star, 'etaOffsetVec');  mp.star.etaOffsetVec = 0;  end
if ~isfield(mp.star, 'weights');  mp.star.weights = 1;  end
if ~isfield(mp.compact.star, 'count');  mp.compact.star.count = 1;  end
if ~isfield(mp.compact.star, 'xiOffsetVec');  mp.compact.star.xiOffsetVec = 0;  end
if ~isfield(mp.compact.star, 'etaOffsetVec');  mp.compact.star.etaOffsetVec = 0;  end
if ~isfield(mp.compact.star, 'weights');  mp.compact.star.weights = 1;  end
if ~isfield(mp.jac.star, 'weights');  mp.jac.star.weights = ones(1, mp.compact.star.count);  end % Spatial weighting in the Jacobian by star

if(isfield(mp.full,'pol_conds')==false);  mp.full.pol_conds = 0;  end %--Vector of which polarization state(s) to use when creating images from the full model. Currently only used with PROPER full models from John Krist.

%--Propagation method
if(isfield(mp,'propMethodPTP')==false);  mp.propMethodPTP = 'fft';  end %--Propagation method for postage stamps around the influence functions. 'mft' or 'fft'

%--FPM Errors
if(isfield(mp.F3, 'full')==false);  mp.F3.full.dummy = 1;  end
if(isfield(mp.F3.full, 'flagErrors')==false);  mp.F3.full.flagErrors = false;  end

%--Vortex or other azithumal, phase-only FPMs
if(isfield(mp.jac, 'mftToVortex')==false);  mp.jac.mftToVortex = false;  end  %--Whether to use MFTs to propagate to/from the vortex FPM
if(isfield(mp.F3, 'VortexSpotDiam')==false);  mp.F3.VortexSpotDiam = 0;  end  %--Diameter of the opaque spot at the center of the vortex. [lambda0/D]
if(isfield(mp.F3, 'VortexSpotOffsets')==false);  mp.F3.VortexSpotOffsets = [0 0];  end  %--Offsets for the opaque spot at the center of the vortex. [lambda0/D]
if(isfield(mp.F3, 'phaseMaskType')==false);  mp.F3.phaseMaskType = 'vortex';  end  % Type of phase FPMs allowed: 'vortex', 'cos', 'sectors', and 'staircase'.
if(isfield(mp.F3, 'NstepStaircase')==false);  mp.F3.NstepStaircase = 6;  end  % Number of discrete steps per 2*pi radians of phase for a staircase phase mask at F3.
if(isfield(mp.F3, 'clocking')==false);  mp.F3.clocking = 0;  end  % Counterclockwise clocking of the phase FPM [degrees].
if(isfield(mp.F3, 'phaseScaleFac')==false);  mp.F3.phaseScaleFac = 1;  end  % Factor to apply to the phase in the phase FPM. Use a vector to add chromaticity to the model. 
if(isfield(mp.F3, 'inVal')==false);  mp.F3.inVal = 0.3;  end  % Inner radius to start the Tukey window for azimuthal phase FPMs. Units of lambda0/D.
if(isfield(mp.F3, 'outVal')==false);  mp.F3.outVal = 5.0;  end  % Out radius to end the Tukey window for azimuthal phase FPMs. Units of lambda0/D.
if(isfield(mp.F3, 'flagDimple')==false);  mp.F3.flagDimple = false;  end
if(isfield(mp.F3, 'roddierphase')==false);  mp.F3.roddierphase = 0.50;  end  % [waves]
if(isfield(mp.F3, 'roddierradius')==false);  mp.F3.roddierradius = 0.53;  end % [lambda/D]

%--HLC FPM materials
if(isfield(mp.F3, 'substrate')==false);  mp.F3.substrate = 'FS';  end % name of substrate material  [FS or N-BK7]
if(isfield(mp.F3, 'metal')==false);  mp.F3.metal = 'nickel';  end % name of the metal used in the FPM  ['nickel']
if(isfield(mp.F3, 'dielectric')==false);  mp.F3.dielectric = 'PMGI';  end % name of the dielectric used in the FPM  ['PMGI' or 'MgF2']

%--Sensitivities to Zernike-Mode Perturbations
if(isfield(mp.full,'ZrmsVal')==false);  mp.full.ZrmsVal = 1e-9;  end %--Amount of RMS Zernike mode used to calculate aberration sensitivities [meters]. WFIRST CGI uses 1e-9, and LUVOIR and HabEx use 1e-10. 
if(isfield(mp.eval,'Rsens')==false);  mp.eval.Rsens = [];   end
if(isfield(mp.eval,'indsZnoll')==false);  mp.eval.indsZnoll = [2,3];   end

%--Deformable mirror settings
if(isfield(mp.dm1,'Nactbeam')==false);  mp.dm1.Nactbeam = mp.dm1.Nact;  end % Number of actuators across the beam (approximate)
if(isfield(mp.dm1,'basisType')==false);  mp.dm1.basisType = 'actuator';  end %--Basis set used for control Options: 'actuator', 'fourier'
if(isfield(mp.dm1,'orientation')==false);  mp.dm1.orientation = 'rot0';  end %--Change to mp.dm1.V orientation before generating DM surface. Options: rot0, rot90, rot180, rot270, flipxrot0, flipxrot90, flipxrot180, flipxrot270
if(isfield(mp.dm1,'fitType')==false);  mp.dm1.fitType = 'linear';  end %--Type of response for displacement vs voltage. Options are 'linear', 'quadratic', and 'fourier2'.
if(isfield(mp.dm1,'dead')==false); mp.dm1.dead = []; end % Vector of linear indices of all dead actuators (those stuck at 0V absolute). This should stay fixed.
if(isfield(mp.dm1,'pinned')==false);  mp.dm1.pinned = [];  end %--Indices of pinned/railed actuators
if(isfield(mp.dm1,'Vpinned')==false);  mp.dm1.Vpinned = zeros(size(mp.dm1.pinned));  end %--(Fixed) relative voltage commands of pinned/railed actuators
if(isfield(mp.dm1,'tied')==false);  mp.dm1.tied = zeros(0, 2);  end %--Indices of tied actuator pairs sharing the same voltage. Two indices per row.
if mp.flagSim
    if(isfield(mp.dm1,'Vmin')==false);  mp.dm1.Vmin = 0;  end %--Min allowed absolute voltage command
    if(isfield(mp.dm1,'Vmax')==false);  mp.dm1.Vmax = 1000;  end %--Max allowed absolute voltage command
else
    if(isfield(mp.dm1,'Vmin')==false);  mp.dm1.Vmin = 0;  end %--Min allowed absolute voltage command
    if(isfield(mp.dm1,'Vmax')==false);  mp.dm1.Vmax = 100;  end %--Max allowed absolute voltage command
end
if(isfield(mp.dm1,'dVnbrLat')==false); mp.dm1.dVnbrLat = mp.dm1.Vmax; end % max voltage difference allowed between laterally-adjacent DM actuators
if(isfield(mp.dm1,'dVnbrDiag')==false); mp.dm1.dVnbrDiag = mp.dm1.Vmax; end % max voltage difference allowed between diagonally-adjacent DM actuators
if(isfield(mp.dm1,'biasMap')==false);  mp.dm1.biasMap = mp.dm1.Vmax/2*ones(mp.dm1.Nact, mp.dm1.Nact);  end  %--Bias voltage. Needed prior to WFSC to allow + and - voltages. Total voltage is mp.dm1.biasMap + mp.dm1.V
if(isfield(mp.dm1,'facesheetFlatmap')==false);  mp.dm1.facesheetFlatmap = mp.dm1.biasMap;  end %--Voltage map that produces a flat DM1 surface. Used when enforcing the neighbor rule.
if(isfield(mp.dm1,'comovingGroups')==false);  mp.dm1.comovingGroups = {};  end % Cell array with each index containing a vector of linear indices for actuators that move together. The vectors can be any length.
if(isfield(mp.dm1,'tolNbrRule')==false); mp.dm1.marginNbrRule = 0.001; end % voltage tolerance used when checking neighbor rule and bound limits. Units of volts.

if(isfield(mp.dm2,'Nactbeam')==false);  mp.dm2.Nactbeam = mp.dm2.Nact;  end % Number of actuators across the beam (approximate)
if(isfield(mp.dm2,'basisType')==false);  mp.dm2.basisType = 'actuator';  end %--Basis set used for control Options: 'actuator', 'fourier'
if(isfield(mp.dm2,'orientation')==false);  mp.dm2.orientation = 'rot0';  end %--Change to mp.dm2.V orientation before generating DM surface. Options: rot0, rot90, rot180, rot270, flipxrot0, flipxrot90, flipxrot180, flipxrot270
if(isfield(mp.dm2,'fitType')==false);  mp.dm2.fitType = 'linear';  end %--Type of response for displacement vs voltage. Options are 'linear', 'quadratic', and 'fourier2'.
if(isfield(mp.dm2,'dead')==false); mp.dm2.dead = []; end % Vector of linear indices of all dead actuators (those stuck at 0V absolute). This should stay fixed.
if(isfield(mp.dm2,'pinned')==false);  mp.dm2.pinned = [];  end %--Indices of pinned/railed actuators
if(isfield(mp.dm2,'Vpinned')==false);  mp.dm2.Vpinned = zeros(size(mp.dm2.pinned));  end %--(Fixed) relative voltage commands of pinned/railed actuators
if(isfield(mp.dm2,'tied')==false);  mp.dm2.tied = zeros(0, 2);  end %--Indices of tied actuator pairs sharing the same voltage. Two indices per row.
if mp.flagSim
    if(isfield(mp.dm2,'Vmin')==false);  mp.dm2.Vmin = 0;  end %--Min allowed absolute voltage command
    if(isfield(mp.dm2,'Vmax')==false);  mp.dm2.Vmax = 1000;  end %--Max allowed absolute voltage command
else
    if(isfield(mp.dm2,'Vmin')==false);  mp.dm2.Vmin = 0;  end %--Min allowed absolute voltage command
    if(isfield(mp.dm2,'Vmax')==false);  mp.dm2.Vmax = 100;  end %--Max allowed absolute voltage command
end
if(isfield(mp.dm2,'dVnbrLat')==false); mp.dm2.dVnbrLat = mp.dm2.Vmax; end % max voltage difference allowed between laterally-adjacent DM actuators
if(isfield(mp.dm2,'dVnbrDiag')==false); mp.dm2.dVnbrDiag = mp.dm2.Vmax; end % max voltage difference allowed between diagonally-adjacent DM actuators
if(isfield(mp.dm2,'biasMap')==false);  mp.dm2.biasMap = mp.dm2.Vmax/2*ones(mp.dm2.Nact, mp.dm2.Nact);  end %--Bias voltage. Needed prior to WFSC to allow + and - voltages. Total voltage is mp.dm2.biasMap + mp.dm2.V
if(isfield(mp.dm2,'facesheetFlatmap')==false);  mp.dm2.facesheetFlatmap = mp.dm2.biasMap;  end %--Voltage map that produces a flat DM2 surface. Used when enforcing the neighbor rule.
if(isfield(mp.dm2,'comovingGroups')==false);  mp.dm2.comovingGroups = {};  end % Cell array with each index containing a vector of linear indices for actuators that move together. The vectors can be any length.
if(isfield(mp.dm2,'tolNbrRule')==false); mp.dm2.marginNbrRule = 0.001; end % voltage tolerance used when checking neighbor rule and bound limits. Units of volts.


%--Control
if(isfield(mp.jac,'zerns')==false); mp.jac.zerns = 1; end %--Zernike modes in Jacobian
if(isfield(mp.jac,'Zcoef')==false); mp.jac.Zcoef = 1; end %--coefficients (i.e., weights) of Zernike modes in Jacobian. Weight for piston is always 1.
if(isfield(mp,'WspatialDef')==false);  mp.WspatialDef = [];  end %--spatial weights for the Jacobian
if(isfield(mp.jac,'minimizeNI')==false); mp.jac.minimizeNI = false; end %--Have EFC minimize normalized intensity instead of intensity
    
%--Estimation
if ~isfield(mp.est, 'probeSchedule'); mp.est.probeSchedule = ProbeSchedule; end %--Schedule of per-iteration values for pairwise probing. Default is empty vectors, meaning they aren't used.
if(isfield(mp.est,'InormProbeMax')==false); mp.est.InormProbeMax = 1e-4; end %--Max allowed probe intensity
if(isfield(mp.est,'Ithreshold')==false); mp.est.Ithreshold = 1e-2; end %--Reduce estimated intensities to this value if they exceed this (probably due to a bad inversion)
if(isfield(mp.scc,'modeCoef')==false); mp.scc.modeCoef = 1; end %--Gain coefficient to apply to the normalized DM basis sets for the empirical SCC calibration.

%--IEFC
if(isfield(mp.iefc,'modeCoef')==false); mp.iefc.modeCoef = 1e-3; end %--Gain coefficient to apply to the normalized DM basis sets for the empirical SCC calibration.
if(isfield(mp.iefc,'probeCoef')==false); mp.iefc.probeCoef = 1e-2; end %--Gain coefficient to apply to the stored probe commands used for IEFC state estimation.
if(isfield(mp.iefc,'probeDM')==false); mp.iefc.probeDM = 1; end %--Which DM to use when probing for IEFC.

%--Performance Evaluation
if(isfield(mp.Fend.eval,'res')==false);  mp.Fend.eval.res = 10;  end % pixels per lambda0/D in compact evaluation model's final focus
mp.mas2lam0D = 1/(mp.lambda0/mp.P1.D*180/pi*3600*1000); %% Conversion factor: milliarcseconds (mas) to lambda0/D
if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0; end % Needed for computing RMS DM surface actuation

%--Training Data: mp.NitrTrain = 5;  %--The number of correction iterations to use per round of training data for the adaptive Jacobian (E-M) algorithm.
%--Zernike sensitivities to 1nm RMS: which noll indices in which annuli, given by mp.eval.indsZnoll and mp.eval.Rsens 
%--Tied actuator pair definitions: See Section with variables mp.dmX.tied for X=1:9
%--Quantization of DM actuation steps based on least significant bit of the
% DAC (digital-analog converter). In height, so called HminStep. If HminStep (minimum step in H) is defined, then quantize the DM voltages
% Variables to define if wanted: mp.dm1.HminStep, mp.dm2.HminStep

%% Detector properties for adding noise to images

% Default values are for the Andor Neo sCMOS detector and testbed flux
if ~isfield(mp, 'flagImageNoise'); mp.flagImageNoise = false; end % whether to include noise in the images
if ~isfield(mp.detector, 'gain'); mp.detector.gain = 1.0; end % detector gain [e-/count]
if ~isfield(mp.detector, 'darkCurrentRate'); mp.detector.darkCurrentRate = 0.015; end % [e-/pixel/second]
if ~isfield(mp.detector, 'readNoiseStd'); mp.detector.readNoiseStd = 1.7; end % standard deviation of Gaussian read noise [e-]
if ~isfield(mp.detector, 'wellDepth'); mp.detector.wellDepth = 3e4; end % well depth of the detector [e-]
if ~isfield(mp.detector, 'peakFluxVec'); mp.detector.peakFluxVec = 1e8 * ones(mp.Nsbp, 1); end % [counts/pixel/second]
if ~isfield(mp.detector, 'tExpVec'); mp.detector.tExpVec = 1.0 * ones(mp.Nsbp, 1); end % exposure times for images in each subband. Set to the unprobed or probed values depending on the situation. [seconds]
if ~isfield(mp.detector, 'tExpUnprobedVec'); mp.detector.tExpUnprobedVec = 1.0 * ones(mp.Nsbp, 1); end % exposure times for unprobed images in each subband [seconds]
if ~isfield(mp.detector, 'tExpProbedVec'); mp.detector.tExpProbedVec = 1.0 * ones(mp.Nsbp, 1); end % exposure times for probed images in each subband [seconds]
if ~isfield(mp.detector, 'Nexp'); mp.detector.Nexp = 1; end % number of exposures to stack for each combined frame

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
