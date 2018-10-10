% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-01 by extracting material from falco_wfsc_loop.m.
% ---------------

function [mp,thput] = falco_compute_thput(mp)

    ImTemp = falco_sim_image_compact_offaxis(mp, mp.thput_eval_x, mp.thput_eval_y,'eval');
    if(mp.flagPlot); figure(324); imagesc(mp.F4.eval.xisDL,mp.F4.eval.etasDL,ImTemp); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end

    maskHM = 0*mp.FP4.eval.RHOS;
    maskHM(ImTemp>=1/2*max(max(ImTemp))) = 1;
    mp.maskHMcore = maskHM.*mp.maskCore;
    % figure(325); imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,mp.maskCore); axis xy equal tight; drawnow;
    thput = sum(ImTemp(mp.maskHMcore==1))/mp.sumPupil*mean(mp.F4.eval.I00);
    fprintf('>=Half-max core throughput within a %.2f lambda/D radius = %.2f%% \tat separation = (%.1f, %.1f) lambda/D.\n',mp.thput_radius,100*thput,mp.thput_eval_x,mp.thput_eval_y);


end %--END OF FUNCTION