% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated broadband image for all wavelengths and
% tip/tilt settings provided.
%  -> Purpose is to be easier to use and faster than re-computing DM
%  surfaces for each call of model_full.m
%
% Modified on 2018-05-17 by A.J. Riggs to be a wrapper function so that
%   parfor can be used.
% Modified on 23 May 2017 by A.J. Riggs to not include any noise sources.
%   Now intended only for coronagraph design. New function name is
%   falco_get_image.m.
% Changed on 19 September 2016 by A.J. Riggs from hcil_* to
%  fpwc_getSimImage.m
% 7 May 2015: Fixed error in which I was zeroing out ImTot instead of
%  stacking it. Moved ImTot=0 outside qq loop and inserted line to divide
%  ImTot by number of images at the end.
%  getSimImage_HCIL_BB_v1.m - 26 Nov. 2014: 
%   -Now uses a flux and exposure time instead of desired avg counts to set
%   number of counts.
% getSimImage_HCIT_BB_v2.m - 23 Sept. 2014: Modified to create an image from several discrete
%  wavelengths instead of just 1. Also added a separate flag for the
%  planet.


function Isum = falco_gen_summed_image(mp, DM)

%--Compute the DM surfaces outside the full model to save lots of time
if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end


Isum = 0; % Initialize summed image

%--Loop over pairs of sub-bandpass and tip/tilt offsets
if(mp.flagParfor)
    parfor si=1:mp.Nsbp
        Itemp = falco_gen_summed_image_func(si, mp, DM);
        Isum = Isum + Itemp*mp.sbp_weights(si);
    end
else
    fprintf('Generating summed image for mode:\n');
    for si=1:mp.Nsbp
        fprintf('\t%4d/%4d\n',si,mp.Nsbp);
        Itemp = falco_gen_summed_image_func(si, mp, DM);
        Isum = Isum + Itemp*mp.sbp_weights(si);
    end
end


% %--Loop over pairs of sub-bandpass and tip/tilt offsets
% if(mp.flagParfor)
%     parfor tsi=1:mp.jac.Nmode
%         Ittlam = falco_gen_summed_image_func(tsi, mp, DM);
%         Isum = Isum + Ittlam;
%     end
% else
%     fprintf('Generating summed image for mode:\n');
%     for tsi=1:mp.jac.Nmode
%         fprintf('\t%4d/%4d\n',tsi,mp.jac.Nmode);
%         Ittlam = falco_gen_summed_image_func(tsi, mp, DM);
%         Isum = Isum + Ittlam;
%     end
% end
    
end %--END OF FUNCTION

