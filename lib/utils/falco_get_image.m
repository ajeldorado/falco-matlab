% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Wrapper function to obtain a real image from a testbed camera or
%  a simulated image with noise using the full model.
%  -For a single sub-bandpass only.
%
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


function Iout = falco_get_image(mp, modvar, DM)


    if(isfield(mp,'flagTestbed'))
        if(mp.flagTestbed) 
            %--This function will take care of image processing (dark subtraction, flat
            %fielding, etc.); crop, rotate, and flip the image as necessary; and normalize the intensity.
            Iout = falco_lab_image(mp, modvar, DM); %--Needs more structures as inputs
            
        else %--Get image in simulation
            Iout = falco_sim_image(mp, modvar, DM); 
        end
    else %--Get image in simulation
        Iout = falco_sim_image(mp, modvar, DM);
    end


end
