% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Initialize arrays to store performance metrics from each WFSC iteration.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% out : structure of output variables

function out = falco_init_storage_arrays(mp)

    out.Itr = 1;
    out.Nitr = mp.Nitr;

    %--Normalization history
    out.I00compactZeros = zeros(mp.Nsbp, 1); % Star normalization for mp.dm1.V = mp.dm2.V = zeros
    out.I00fullZeros = zeros(mp.Nsbp, mp.Nwpsbp); % Star normalization for mp.dm1.V = mp.dm2.V = zeros
    out.I00compactHist = zeros(mp.Nsbp, mp.Nitr+1);
    out.I00fullHist = zeros(mp.Nsbp, mp.Nwpsbp, mp.Nitr+1);
    out.I00compactRatioHist = zeros(mp.Nsbp, mp.Nitr+1); % Stellar peak ratio
    out.I00fullRatioHist = zeros(mp.Nsbp, mp.Nwpsbp, mp.Nitr+1); % Stellar peak ratio

    %--EFC regularization history
    out.log10regHist = zeros(mp.Nitr, 1);  % Keeping for backwards compatibility
    out.ctrl.log10regHist = zeros(mp.Nitr, 1);
    out.ctrl.dmfacHist = zeros(mp.Nitr, 1);    
    out.ctrl.relinHist = false(mp.Nitr, 1);

    %--Peak-to-Valley DM voltages
    out.dm1.Vpv = zeros(mp.Nitr, 1);
    out.dm2.Vpv = zeros(mp.Nitr, 1);
    out.dm8.Vpv = zeros(mp.Nitr, 1);
    out.dm9.Vpv = zeros(mp.Nitr, 1);

    %--Peak-to-Valley DM surfaces
    out.dm1.Spv = zeros(mp.Nitr, 1);
    out.dm2.Spv = zeros(mp.Nitr, 1);
    out.dm8.Spv = zeros(mp.Nitr, 1);
    out.dm9.Spv = zeros(mp.Nitr, 1);

    %--RMS DM surfaces
    out.dm1.Srms = zeros(mp.Nitr, 1);
    out.dm2.Srms = zeros(mp.Nitr, 1);
    out.dm8.Srms = zeros(mp.Nitr, 1);
    out.dm9.Srms = zeros(mp.Nitr, 1);
    
    %--Delta DM Surface Metrics
    out.dm1.DeltaSpv = zeros(mp.Nitr, 1);
    out.dm2.DeltaSpv = zeros(mp.Nitr, 1);
    out.dm1.DeltaSrms = zeros(mp.Nitr, 1);
    out.dm2.DeltaSrms = zeros(mp.Nitr, 1);
    
    %--DM constraints
    % pinned is the vector of linear indices for dead, railed, or otherwise stuck actuators.
    % Vpinned is the vector of relative voltages for all pinned actuators.
    % comovingGroups is the cell array of vectors of (linear indices of) comoving actuator
    % groups, either because of electrical ties or neighbor rule ties. 
    out.dm1.pinned = {};
    out.dm1.Vpinned = {};
    out.dm1.comovingGroups = {};
    out.dm1.Npinned = zeros(mp.Nitr, 1);
    out.dm1.Ncomoving = zeros(mp.Nitr, 1);
    
    out.dm2.Vpinned = {};    
    out.dm2.pinned = {};
    out.dm2.comovingGroups = {};
    out.dm2.Npinned = zeros(mp.Nitr, 1);
    out.dm2.Ncomoving = zeros(mp.Nitr, 1);
    
    %--Sensitivities Zernike-Mode Perturbations
    Nannuli = size(mp.eval.Rsens, 1);
    Nzern = length(mp.eval.indsZnoll);
    out.Zsens = zeros(Nzern, Nannuli, mp.Nitr);

    %--Store the DM commands at each iteration
    if(isfield(mp, 'dm1')); if(isfield(mp.dm1, 'V'));  out.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, mp.Nitr+1);  end; end
    if(isfield(mp, 'dm2')); if(isfield(mp.dm2, 'V'));  out.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, mp.Nitr+1); end; end
    if(isfield(mp, 'dm8')); if(isfield(mp.dm8, 'V'));  out.dm8.Vall = zeros(mp.dm8.NactTotal, mp.Nitr+1); end; end
    if(isfield(mp, 'dm9')); if(isfield(mp.dm9, 'V'));  out.dm9.Vall = zeros(mp.dm9.NactTotal, mp.Nitr+1); end; end

    %--Delta electric field performance metrics
    out.complexProjection = zeros(mp.Nitr-1, mp.Nsbp); % Metric to compare magnitude of the correction step taken to the expected one
    out.complexCorrelation = zeros(mp.Nitr-1, mp.Nsbp); % Metric to compare the morphology of the delta E-field estimated vs expected in the model
    
    %--Intensity history at each iteration
    if mp.flagFiber
        out.InormFiberHist = zeros(mp.Nitr+1, mp.Fend.Nfiber);
    end
    if mp.flagFiber
        out.InormHist = zeros(mp.Nitr+1, mp.Fend.Nfiber); % Measured, mean raw NI in correction region of dark hole.
        out.IrawCorrHist = zeros(mp.Nitr+1, mp.Fend.Nfiber); % Measured, mean raw NI in correction region of dark hole.
        out.IrawScoreHist = zeros(mp.Nitr+1, mp.Fend.Nfiber); % Measured, mean raw NI in scoring region of dark hole.
        out.IestCorrHist = zeros(mp.Nitr, mp.Fend.Nfiber); % Mean estimated coherent NI in correction region of dark hole.
        out.IestScoreHist = zeros(mp.Nitr, mp.Fend.Nfiber); % Mean estimated coherent NI in scoring region of dark hole.
        out.IincoCorrHist = zeros(mp.Nitr, mp.Fend.Nfiber); % Mean estimated incoherent NI in correction region of dark hole.
        out.IincoScoreHist = zeros(mp.Nitr, mp.Fend.Nfiber); % Mean estimated incoherent NI in scoring region of dark hole.
    else
        out.InormHist = zeros(mp.Nitr+1, 1); % Measured, mean raw NI in correction region of dark hole.
        out.IrawCorrHist = zeros(mp.Nitr+1, 1); % Measured, mean raw NI in correction region of dark hole.
        out.IrawScoreHist = zeros(mp.Nitr+1, 1); % Measured, mean raw NI in scoring region of dark hole.
        out.IestCorrHist = zeros(mp.Nitr, 1); % Mean estimated coherent NI in correction region of dark hole.
        out.IestScoreHist = zeros(mp.Nitr, 1); % Mean estimated coherent NI in scoring region of dark hole.
        out.IincoCorrHist = zeros(mp.Nitr, 1); % Mean estimated incoherent NI in correction region of dark hole.
        out.IincoScoreHist = zeros(mp.Nitr, 1); % Mean estimated incoherent NI in scoring region of dark hole.
    end
    
    out.normIntMeasCorr = zeros(mp.Nitr, mp.Nsbp); % Measured raw NI in correction region of dark hole.
    out.normIntMeasScore = zeros(mp.Nitr, mp.Nsbp); % Measured raw NI in scoring region of dark hole.
    out.normIntModCorr = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated modulated NI in correction region of dark hole.
    out.normIntModScore = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated modulated NI in scoring region of dark hole.
    out.normIntUnmodCorr = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated unmodulated NI in correction region of dark hole.
    out.normIntUnmodScore = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated unmodulated NI in correction region of dark hole.
    
    %--Storage array for throughput at each iteration
    if mp.flagFiber
        out.thput = zeros(mp.Nitr+1, mp.Fend.Nfiber);
    else
        out.thput = zeros(mp.Nitr+1, 1);
    end
    
    %--Variables related to final image
    out.Fend.res = mp.Fend.res;
    out.Fend.xisDL = mp.Fend.xisDL;
    out.Fend.etasDL = mp.Fend.etasDL;
    out.Fend.scoreInCorr = mp.Fend.scoreInCorr;
    out.Fend.corr.maskBool = mp.Fend.corr.maskBool;
    out.Fend.score.maskBool = mp.Fend.score.maskBool;
    
    out.serialDateVec = zeros(mp.Nitr, 1); % start time of each iteration as a serial date number

end
