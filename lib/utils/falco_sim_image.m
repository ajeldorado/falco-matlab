% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated, single-sub-bandpass image with optional planet light added.
% -Measurement noise terms not added in yet
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

function Iout = falco_sim_image(mp, modvar, DM)

modvar.flagGetNormVal = false; % False for all cases

if( strcmpi(modvar.whichSource,'offaxis') ) %--For throughput calculations by moving star off axis.
    modvar.wpsbpIndex = mp.wi_ref; %--At center wavelength only
    Eout = model_full(mp, DM, modvar);
    Iout = abs(Eout).^2;
    
else %--Regular
    Isum = 0; % Initialize unstacked image
    
    for wi=1:mp.Nwpsbp  % Add intensities from all incoherent sources separately
        % Starlight
        modvar.wpsbpIndex = wi;
        modvar.whichSource = 'star';
        Eout = model_full(mp, DM, modvar);
        Isum = Isum + (abs(Eout).^2);
    end

    Iout = Isum/mp.Nwpsbp;
end