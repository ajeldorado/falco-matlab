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

function [mp,thput,varargout] = falco_compute_thput(mp)

if mp.flagFiber
    fprintf('Computing throughput over %i fibers (one propagation per fiber) \n',mp.Fend.Nfiber);
    thput = zeros(mp.Fend.Nfiber,1);
    for ii=1:mp.Fend.Nfiber
        [~, ImSimCompact] = falco_sim_image_compact_offaxis(mp, mp.Fend.x_fiber(ii), mp.Fend.y_fiber(ii));
        thput(ii) = mean(ImSimCompact(ii))/mp.sumPupil*mean(mp.Fend.eval.I00_fiber(ii));
        fprintf('Fiber (%i/%i) throughput = %.2f%% \tat separation = (%.1f, %.1f) lambda/D.\n', ii, mp.Fend.Nfiber, 100*thput(ii), mp.Fend.x_fiber(ii), mp.Fend.y_fiber(ii));
    end
else

    ImSimCompact = falco_sim_image_compact_offaxis(mp, mp.thput_eval_x, mp.thput_eval_y, 'eval');
%     if(mp.flagPlot); figure(324); imagesc(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL,ImSimCompact); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end

    switch lower(mp.thput_metric)
        
        case{'hmi'} %--Absolute energy within half-max isophote(s)
            
            maskHM = 0*mp.Fend.eval.RHOS;
            maskHM(ImSimCompact>=1/2*max(max(ImSimCompact))) = 1;
            % figure(325); imagesc(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL,maskHM); axis xy equal tight; drawnow;
            thput = sum(ImSimCompact(maskHM==1))/mp.sumPupil*mean(mp.Fend.eval.I00);
            fprintf('Core throughput within the half-max isophote(s) = %.2f%% \tat separation = (%.1f, %.1f) lambda0/D.\n',100*thput,mp.thput_eval_x,mp.thput_eval_y);
            
        case{'ee','e.e.'} %--Absolute energy encircled within a given radius

            % (x,y) location [lambda_c/D] in dark hole at which to evaluate throughput
            maskEE  = 0*mp.Fend.eval.RHOS;
            maskEE(mp.Fend.eval.RHOS<=mp.thput_radius) = 1;
            % figure(325); imagesc(mp.Fend.eval.xisDL,mp.Fend.eval.etasDL,maskEE); axis xy equal tight; drawnow;
            thput = sum(ImSimCompact(maskEE==1))/mp.sumPupil*mean(mp.Fend.eval.I00);
            fprintf('E.E. throughput within a %.2f lambda/D radius = %.2f%% \tat separation = (%.1f, %.1f) lambda/D.\n',mp.thput_radius,100*thput,mp.thput_eval_x,mp.thput_eval_y);
    end
    
end
varargout{1} = ImSimCompact;

end %--END OF FUNCTION