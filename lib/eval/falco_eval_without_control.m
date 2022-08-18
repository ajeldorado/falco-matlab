% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Evaluate the coronagraphic system with just estimation and no control.
%
% Requires falco_flesh_out_workspace() to have been run first.


function outSingle = falco_eval_without_control(mp)

    Nitr = mp.Nitr;

    mp.Nitr = 1; % change only for this function
    Itr = 1;
    ev.Itr = Itr;
    
    outSingle = falco_init_storage_arrays(mp);

    mp = falco_compute_psf_norm_factor(mp);
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
    outSingle.thput(Itr) = thput;
    
    ev.dummy = 1;
    jacStruct = [];
    ev = falco_est(mp, ev, jacStruct);
    
    outSingle = falco_store_intensities(mp, outSingle, ev, Itr);
    
    outSingle = falco_compute_dm_stats(mp, outSingle, Itr);
    
    fprintf('Mean NI:\t\t\t %.2e \n', outSingle.InormHist(Itr))
        
    mp.Nitr = Nitr; % reset
    
end
