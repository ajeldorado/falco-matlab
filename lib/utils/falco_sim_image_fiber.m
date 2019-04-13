% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a simulated, single-sub-bandpass image through the 
% fiber.

function [Iout, Ifiber] = falco_sim_image_fiber(mp, modvar, DM)

modvar.flagGetNormVal = false; % False for all cases

if( strcmpi(modvar.whichSource,'offaxis') ) %--For throughput calculations by moving star off axis.
    modvar.wpsbpIndex = mp.wi_ref; %--At center wavelength only
    [Eout, Efiber] = model_full(mp, DM, modvar);
    Iout = abs(Eout).^2;
    Ifiber = abs(Efiber).^2;
    
else %--Regular
    Isum = 0; % Initialize unstacked image
    
    for wi=1:mp.Nwpsbp  % Add intensities from all incoherent sources separately
        % Starlight
        modvar.wpsbpIndex = wi;
        modvar.whichSource = 'star';
        [Eout, Efiber] = model_full(mp, DM, modvar);
        Isum = Isum + (abs(Eout).^2);
        Isumfiber = Isumfiber + (abs(Efiber).^2);
    end

    Iout = Isum/mp.Nwpsbp;
    Ifiber = Isumfiber/mp.Nwpsbp;
end