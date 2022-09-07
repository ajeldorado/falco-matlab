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
    mp.allowedEstimators = {'perfect', 'pairwise', 'pairwise-square', 'pwp-bp-square', 'pairwise-rect', 'pwp-bp', 'pwp-kf', 'scc'};
    mp.allowedControllers = {'gridsearchefc', 'plannedefc'};
    
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
    
    %--Check valid combinations of coronagraph type and layout
    if strcmpi(mp.layout, 'fpm_scale') && ~strcmpi(mp.coro, 'HLC')
        error("Error: mp.layout is only allowed as 'fpm_scale' when mp.coro = 'HLC'.")
    end
    
    %--Check estimator
    mp.estimator = lower(mp.estimator);
    if ~any(strcmp(mp.allowedEstimators, mp.estimator))
        error('Error: %s is not an allowed value of mp.estimator.', mp.estimator)
    end
    
    %--Check controller
    mp.controller = lower(mp.controller);
    if ~any(strcmp(mp.allowedControllers, mp.controller))
        error('Error: %s is not an allowed value of mp.controller.', mp.controller)
    end
    
end
