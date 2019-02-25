% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated broadband image for all wavelengths.
% Tip/tilt no longer included.
%  -> Purpose is to be easier to use and faster than re-computing DM
%  surfaces for each call of model_full.m
%
% Modified on 2018-08-31 by A.J. Riggs to use the structure mp.entireBWimage and not to include tip/tilt
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


function [Isum, IfiberSum] = falco_get_summed_image(mp)

%--Compute the DM surfaces outside the full model to save lots of time
if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end
if(any(mp.dm_ind==2)); mp.dm2.surfM = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,mp.dm2.NdmPad); end
if(any(mp.dm_ind==9)); mp.dm9.phaseM = falco_dm_surf_from_cube(mp.dm9,mp.dm9); end

Isum = 0; % Initialize image
IfiberSum = 0;

for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    
    for wi=1:mp.Nwpsbp
        modvar.whichSource = 'star';
        modvar.wpsbpIndex = wi;
        [Etemp, EfiberTemp] = model_full(mp, modvar);
        Isum = Isum + (abs(Etemp).^2)*mp.sbp_weights(si)*mp.full.lambda_weights(wi);
        if(mp.flagFiber)
            IfiberSum = IfiberSum + (abs(EfiberTemp).^2)*mp.sbp_weights(si)*mp.full.lambda_weights(wi);
        end
    end 
end

end %--END OF FUNCTION

