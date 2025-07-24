% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Wrapper for the simplified optical models used for the fast Jacobian calculation.
%  The first-order derivative of the DM pokes are propagated through the system.
%  Does not include unknown aberrations/errors that are in the full model.
%  For the SCC, the Jacobians are computed empirically from images.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% jacStruct : structure containing a separate control Jacobian for each DM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select which optical layout's Jacobian model to use and get the output E-field
% INPUTS
% ------
% mp : structure of model parameters
% vals_list : structure containing combinations of:
%     -imode : index of the DM, subband, and Zernike mode
%     -whichDM : which DM number
%     -iactSubset : jacobian index of actuator number
%
% OUTPUTS
% -------
% jacModeAct = Jacobian for the specified combo of DM, subband, and Zernike mode, and actuator.
%
function jacModeAct = model_Jacobian_middle_layer(mp, vals_list, icombo)

    imode = vals_list(1, icombo); %--index for Zernike-&-subbandpass pair
    whichDM = vals_list(2, icombo); %--number of the specified DM
    iactSubset = vals_list(3, icombo); %--jacobian index of actuator number
    
    % Convert to index of specified actuator on that DM
    if whichDM == 1; iact = mp.dm1.act_ele(iactSubset); end
    if whichDM == 2; iact = mp.dm2.act_ele(iactSubset); end
    if whichDM == 8; iact = mp.dm8.act_ele(iactSubset); end
    if whichDM == 9; iact = mp.dm9.act_ele(iactSubset); end

    switch upper(mp.coro)
        
        case{'LC', 'APLC', 'SPLC', 'FLC', 'HLC'}

            jacModeAct = model_Jacobian_lyot(mp, imode, whichDM, iact); 
            
        case{'VORTEX', 'VC'}
            jacModeAct = model_Jacobian_VC(mp, imode, whichDM, iact);
            
        case{'EHLC'} %--Extended HLC: DMs, extended FPM with nickel and dielectric modulation, and LS.
            jacModeAct = model_Jacobian_EHLC(mp, imode, whichDM, iact);
            
        otherwise
            error(sprintf('No Jacobian function is defined for mp.coro = %s', mp.coro));        
    end

end %--END OF FUNCTION
