% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Wrapper function for the different estimator functions.
%
% INPUTS
% ------
% mp : structure of model parameters
% out : structure of output variables
% ev : structure of estimator variables
% jacStruct : structure containing the Jacobians
%
% OUTPUTS
% -------
% ev : structure of estimator variables

function [ev, out] = falco_est(mp, ev, out, jacStruct)

    if ~mp.est.flagUseJac
        clear jacStruct
    end

    switch lower(mp.estimator)
        case{'perfect'}
            ev.Eest  = falco_est_perfect_Efield_with_Zernikes(mp);
            ev.IincoEst = zeros(size(ev.Eest));
            ev.Im = falco_get_summed_image(mp);
            
        case{'pwp-bp-square', 'pwp-bp', 'pwp-kf'}
            if(mp.flagFiber && mp.flagLenslet)
				if mp.est.flagUseJac
					ev = falco_est_pairwise_probing_fiber(mp, jacStruct);
                else
					ev = falco_est_pairwise_probing_fiber(mp);
				end
			else
				if mp.est.flagUseJac
					ev = falco_est_pairwise_probing(mp, ev, jacStruct);
                else
					ev = falco_est_pairwise_probing(mp, ev);
				end
            end
    end
    
    %% Calculate the averaged coherent and incoherent intensities
    
    % Apply subband weights and then sum over subbands
    Iest = abs(ev.Eest).^2;
    Iinco = ev.IincoEst;
    for iMode = 1:mp.jac.Nmode
        sbpIndex = mp.jac.sbp_inds(iMode);
        Iest(:, iMode) = mp.sbp_weights(sbpIndex) * Iest(:, iMode);
        Iinco(:, iMode) = mp.sbp_weights(sbpIndex) * Iinco(:, iMode);
    end
    IestAllBands = sum(Iest, 2);
    IincoAllBands = sum(Iinco, 2);
    
    ev.IestCorrMean = mean(IestAllBands);
    ev.IincoCorrMean = mean(IincoAllBands);
    
    % Put intensities back into 2-D arrays to use correct indexing of scoring region.
    Iest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Iest2D(mp.Fend.corr.maskBool) = IestAllBands(:);
    ev.IestScoreMean = mean(Iest2D(mp.Fend.score.maskBool));
    
    Iinco2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Iinco2D(mp.Fend.corr.maskBool) = IincoAllBands(:);
    ev.IincoScoreMean = mean(Iinco2D(mp.Fend.score.maskBool));
    
    out.IestCorrHist(Itr) = ev.IestCorrMean;
    out.IestScoreHist(Itr) = ev.IestScoreMean;
    out.IincoCorrHist(Itr) = ev.IincoCorrMean;
    out.IincoScoreHist(Itr) = ev.IincoScoreMean;
    out.InormHist(Itr) = mean(ev.Im(mp.Fend.corr.maskBool));
    out.IrawCorrHist(Itr) = mean(ev.Im(mp.Fend.corr.maskBool));
    out.IrawScoreHist(Itr) = mean(ev.Im(mp.Fend.score.maskBool));

end
