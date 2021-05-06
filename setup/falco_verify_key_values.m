% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Verify that some key variables have allowed values.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_verify_key_values(mp)

    mp.allowedCenterings = {'pixel', 'interpixel'};
    mp.allowedCoronagraphTypes = {'VC', 'VORTEX', 'LC', 'APLC', 'FLC', 'SPLC', 'HLC'};
    mp.allowedLayouts = {'fourier', 'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'};

    %--Check centering
    mp.centering = lower(mp.centering);
    if ~any(strcmp(mp.allowedCenterings, mp.centering))
        error('Error: %s is not an allowed value of mp.centering.', mp.centering)
    end
    
    %--Check coronagraph type
    mp.coro = upper(mp.coro);
    if ~any(strcmp(mp.allowedCoronagraphTypes, mp.coro))
        error('Error: %s is not an allowed value of mp.coro.', mp.coro)
    end
    
    %--Check optical layout
    mp.layout = lower(mp.layout);
    if ~any(strcmp(mp.allowedLayouts, mp.layout))
        error('Error: %s is not an allowed value of mp.layout.', mp.layout)
    end

end
