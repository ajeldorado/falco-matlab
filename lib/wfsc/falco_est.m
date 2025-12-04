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
% ev : structure of estimator variables
% jacStruct : structure containing the Jacobians
%
% OUTPUTS
% -------
% ev : structure of estimator variables

function ev = falco_est(mp, ev, jacStruct)
    
    if ~mp.est.flagUseJac; clear jacStruct; end % save RAM

    switch lower(mp.estimator)
        case{'perfect'}
            ev  = falco_est_perfect_Efield_with_Zernikes(mp);
            if mp.flagFiber
                [ev.Im,ev.Ifiber] = falco_get_summed_image(mp);
            else
                ev.Im = falco_get_summed_image(mp);
            end
                
            
        case{'pwp-bp-square', 'pwp-bp', 'pairwise', 'pairwise-square', ...
                'pairwise-rect', 'pwp-kf', 'pairwise-kf', 'pairwise-rect-kf',...
                'pwp-fiber'}
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
            
        case{'borde-traub', 'bt', 'bt-kf', 'bt-rect', 'bt-rect-kf'}
            ev = falco_est_borde_traub_probing(mp, ev);
            
        case{'scc'}
            ev = falco_est_scc(mp);
            if mp.flagFiber
                [ev.Im,ev.Ifiber] = falco_get_summed_image(mp);
            else
                ev.Im = falco_get_summed_image(mp);
            end

        case{'iefc'}
            ev = falco_est_iefc(mp);
            if mp.flagFiber
                [ev.Im,ev.Ifiber] = falco_get_summed_image(mp);
            else
                ev.Im = falco_get_summed_image(mp);
            end

        case{'ekf_maintenance'}
            
            if ev.Itr == 1
                disp('starting ekf initialization')
                ev = initialize_ekf_maintenance(mp, ev, jacStruct);
                disp('done ekf initialization')
            end
            
            ev = falco_est_ekf_maintenance(mp,ev,jacStruct);
    end

end
