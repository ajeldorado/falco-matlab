% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Evaluate the coronagraphic system with just estimation and no control.

function outSingle = falco_eval_without_control(mp)

    Nitr = mp.Nitr;

    mp.Nitr = 1; % change only for this function
    Itr = 1;
    
    outSingle = falco_init_storage_arrays(mp);

    mp = falco_compute_psf_norm_factor(mp);
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
    outSingle.thput(Itr) = thput;
    
    ev.dummy = 1;
    jacStruct = [];
    ev = falco_est(mp, ev, jacStruct);
    
    outSingle = store_intensities(mp, outSingle, ev, Itr);
    
    outSingle = falco_compute_dm_stats(mp, outSingle, Itr);
    
    fprintf('Mean NI:\t\t\t %.2e \n', outSingle.InormHist(Itr))
    
    %% Wrap up
    
    mp.Nitr = Nitr; % reset
    
end


function out = store_intensities(mp, out, ev, Itr)

    %% Calculate the average measured, coherent, and incoherent intensities
    
    % Apply subband weights and then sum over subbands
    Iest = abs(ev.Eest).^2;
    Iinco = ev.IincoEst;
    for iMode = 1:mp.jac.Nmode
        iSubband = mp.jac.sbp_inds(iMode);
        Iest(:, iMode) = mp.sbp_weights(iSubband) * Iest(:, iMode);
        Iinco(:, iMode) = mp.sbp_weights(iSubband) * Iinco(:, iMode);
    end
    IestAllBands = sum(Iest, 2);
    IincoAllBands = sum(Iinco, 2);
    
    
    % Put intensities back into 2-D arrays to use correct indexing of scoring region.
    % Modulated
    Iest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Iest2D(mp.Fend.corr.maskBool) = IestAllBands(:);
    out.IestScoreHist(Itr) = mean(Iest2D(mp.Fend.score.maskBool));
    out.IestCorrHist(Itr) = mean(Iest2D(mp.Fend.corr.maskBool));
    % Unmodulated
    Iinco2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Iinco2D(mp.Fend.corr.maskBool) = IincoAllBands(:);
    out.IincoScoreHist(Itr) = mean(Iinco2D(mp.Fend.score.maskBool));
    out.IincoCorrHist(Itr) = mean(Iinco2D(mp.Fend.corr.maskBool));

    % Measured
    out.IrawScoreHist(Itr) = mean(ev.Im(mp.Fend.score.maskBool));
    out.IrawCorrHist(Itr) = mean(ev.Im(mp.Fend.corr.maskBool));
    out.InormHist(Itr) = out.IrawCorrHist(Itr); % InormHist is a vestigial variable
    
    %% Calculate the measured, coherent, and incoherent intensities by subband
    
    % measured intensities
    for iSubband = 1:mp.Nsbp
        imageMeas = squeeze(ev.imageArray(:, :, 1, iSubband));
        out.normIntMeasCorr(Itr, iSubband) = mean(imageMeas(mp.Fend.corr.maskBool));
        out.normIntMeasScore(Itr, iSubband) = mean(imageMeas(mp.Fend.score.maskBool));
        clear im        
    end
    
    % estimated
    for iMode = 1:mp.jac.Nmode
        
        imageModVec = abs(ev.Eest(:, iMode)).^2;
        imageUnmodVec = ev.IincoEst(:, iMode);
        
        out.normIntModCorr(Itr, iMode) = mean(imageModVec);
        out.normIntModScore(Itr, iMode) = mean(imageModVec(mp.Fend.scoreInCorr));
        
        out.normIntUnmodCorr(Itr, iMode) = mean(imageUnmodVec);
        out.normIntUnmodScore(Itr, iMode) = mean(imageUnmodVec(mp.Fend.scoreInCorr));
    end
    
end
