% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Initialize arrays to store data from each WFSC iteration.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% out : structure of output variables

function out = falco_init_storage_arrays(mp)

    %--EFC regularization history
    out.Nitr = mp.Nitr;
    out.log10regHist = zeros(mp.Nitr, 1);

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
    out.InormHist = zeros(mp.Nitr+1, 1); % Measured, mean raw NI in correction region of dark hole.
    out.IrawCorrHist = zeros(mp.Nitr+1, 1); % Measured, mean raw NI in correction region of dark hole.
    out.IrawScoreHist = zeros(mp.Nitr+1, 1); % Measured, mean raw NI in scoring region of dark hole.
    out.IestCorrHist = zeros(mp.Nitr, 1); % Mean estimated coherent NI in correction region of dark hole.
    out.IestScoreHist = zeros(mp.Nitr, 1); % Mean estimated coherent NI in scoring region of dark hole.
    out.IincoCorrHist = zeros(mp.Nitr, 1); % Mean estimated incoherent NI in correction region of dark hole.
    out.IincoScoreHist = zeros(mp.Nitr, 1); % Mean estimated incoherent NI in scoring region of dark hole.
    
    out.normIntMeasCorr = zeros(mp.Nitr, mp.Nsbp); % Measured raw NI in correction region of dark hole.
    out.normIntMeasScore = zeros(mp.Nitr, mp.Nsbp); % Measured raw NI in scoring region of dark hole.
    out.normIntModCorr = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated modulated NI in correction region of dark hole.
    out.normIntModScore = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated modulated NI in scoring region of dark hole.
    out.normIntUnmodCorr = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated unmodulated NI in correction region of dark hole.
    out.normIntUnmodScore = zeros(mp.Nitr, mp.Nsbp*mp.compact.star.count); % Estimated unmodulated NI in correction region of dark hole.
    
    %--Storage array for throughput at each iteration
    out.thput = zeros(mp.Nitr+1, 1);
    
    %--Variables related to final image
    out.Fend.res = mp.Fend.res;
    out.Fend.xisDL = mp.Fend.xisDL;
    out.Fend.etasDL = mp.Fend.etasDL;
    out.Fend.scoreInCorr = mp.Fend.scoreInCorr;
    out.Fend.corr.maskBool = mp.Fend.corr.maskBool;
    out.Fend.score.maskBool = mp.Fend.score.maskBool;
    
    out.serialDateVec = zeros(mp.Nitr, 1); % start time of each iteration as a serial date number

end
