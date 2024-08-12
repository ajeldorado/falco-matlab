% Copyright 2018, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Measure the cross term of the E-field and the probe for all subbands and stars. Used for IEFC.
%
% Intended for one-off estimates of Delta I. For computing Delta I for each
% mode in a testbed in broadband, do not loop over this function because it
% would require changing the subbands many times and cause unnecessary
% wear on the laser wavelength selector (e.g., a VARIA).
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% ev : structure of estimation variables

function ev = falco_est_iefc(mp)

    
    Nprobes = size(mp.iefc.probeCube, 3);
    
    ev.Eest = zeros(Nprobes*mp.Fend.corr.Npix, mp.jac.Nmode);
    ev.IincoEst = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
    ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    ev.imageArray = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1+2*Nprobes, mp.jac.Nmode);

    for iJacMode = 1:mp.jac.Nmode

        for iProbe = 1:Nprobes
            
            % Get the unprobed image 
            sbpIndex = mp.jac.sbp_inds(iJacMode);
            ev.imageArray(:, :, 1, iJacMode) = falco_get_sbp_image(mp, sbpIndex);
            
            % Get the probed images and E-field estimate 
            mp.isProbing = true;
            dV = mp.iefc.probeCoef * mp.iefc.probeCube(:, :, iProbe);
            iefcStruct = falco_est_delta_intensity(mp, mp.iefc.probeDM, dV, iJacMode);
            ev.Eest((iProbe-1)*mp.Fend.corr.Npix+1:iProbe*mp.Fend.corr.Npix, iJacMode) = iefcStruct.DeltaI / mp.iefc.probeCoef;
            ev.imageArray(:, :, 1+2*(iProbe-1)+1, iJacMode) = iefcStruct.imagePlus;
            ev.imageArray(:, :, 1+2*(iProbe-1)+2, iJacMode) = iefcStruct.imageMinus;
            mp.isProbing = false;  
            
        end      
         
    end
   
end %--END OF FUNCTION
