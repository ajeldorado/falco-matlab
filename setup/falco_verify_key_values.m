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

    allowedCenterings = {'pixel', 'interpixel'};
    allowedCoronagraphTypes = {'VC', 'VORTEX', 'LC', 'APLC', 'FLC', 'SPLC', 'HLC'};
    allowedLayouts = {'fourier', 'fpm_scale', 'proper', 'roman_phasec_proper', 'wfirst_phaseb_proper'};

    %--Check centering
    mp.centering = lower(mp.centering);
    if ~any(strcmp(allowedCenterings, mp.centering))
        error('Error: %s is not an allowed value of mp.centering.', mp.centering)
    end
    
    %--Check coronagraph type
    mp.coro = upper(mp.coro);
    if ~any(strcmp(allowedCoronagraphTypes, mp.coro))
        error('Error: %s is not an allowed value of mp.coro.', mp.coro)
    end
    
    %--Check optical layout
    mp.layout = lower(mp.layout);
    if ~any(strcmp(allowedLayouts, mp.layout))
        error('Error: %s is not an allowed value of mp.layout.', mp.layout)
    end

end
