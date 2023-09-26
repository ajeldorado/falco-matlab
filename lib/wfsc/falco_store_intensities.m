% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Calculate the average measured, coherent, and incoherent intensities after
% running the estimator.
%
% Requires falco_flesh_out_workspace() and the estimator to have been run first.

function out = falco_store_intensities(mp, out, ev, Itr)
    
    if strcmpi(mp.estimator, 'iefc')
        if ~mp.flagFiber
            ev.Eest = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
        else
            ev.Eest = zeros(mp.Fend.Nfiber, mp.jac.Nmode);
        end
    end

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
    if ~mp.flagFiber
        Iest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
        Iest2D(mp.Fend.corr.maskBool) = IestAllBands(:);
        out.IestScoreHist(Itr) = mean(Iest2D(mp.Fend.score.maskBool));
        out.IestCorrHist(Itr) = mean(Iest2D(mp.Fend.corr.maskBool));
        % Unmodulated
        Iinco2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
        Iinco2D(mp.Fend.corr.maskBool) = IincoAllBands(:);
        out.IincoScoreHist(Itr) = mean(Iinco2D(mp.Fend.score.maskBool));
        out.IincoCorrHist(Itr) = mean(Iinco2D(mp.Fend.corr.maskBool));
    end

    % Measured
    out.IrawScoreHist(Itr) = mean(ev.Im(mp.Fend.score.maskBool));
    out.IrawCorrHist(Itr) = mean(ev.Im(mp.Fend.corr.maskBool));
    out.InormHist(Itr) = out.IrawCorrHist(Itr); % InormHist is a vestigial variable
    if mp.flagFiber
        out.InormFiberHist(Itr,:) = ev.Ifiber;
    end

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
        if mp.flagFiber
            out.normIntModScore(Itr, iMode) = mean(imageModVec);
        else
            out.normIntModScore(Itr, iMode) = mean(imageModVec(mp.Fend.scoreInCorr));
        end

        out.normIntUnmodCorr(Itr, iMode) = mean(imageUnmodVec);
        if mp.flagFiber
            out.normIntUnmodScore(Itr, iMode) = mean(imageUnmodVec);
        else
            out.normIntUnmodScore(Itr, iMode) = mean(imageUnmodVec(mp.Fend.scoreInCorr));
        end
    end
    
end
