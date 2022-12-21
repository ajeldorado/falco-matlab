% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Empirically compute the Jacobian for the SCC.
%
% INPUTS
% ------
% mp : structure of model parameters
% whichDM : which DM number
%
% OUTPUTS
% -------
% jac : Jacobian for the specified DM.

function jac = model_Jacobian_SCC(mp, whichDM)

    if whichDM == 1
        
        jac = zeros(mp.Fend.corr.Npix, mp.dm1.Nele, mp.jac.Nmode); %--Initialize the Jacobian
        V0 = mp.dm1.V;

        for ibm = 1:mp.dm1.NbasisModes

            fprintf('Empirically obtaining DM1 response matrix for basis mode: %d/%d\n', ibm, mp.dm1.NbasisModes);
            dV = mp.scc.modeCoef * mp.dm1.basisCube(:, :, ibm);
            mp.dm1.V = V0 + dV;
            evTempPos = falco_est_scc(mp);
            mp.dm1.V = V0 - dV;
            evTempNeg = falco_est_scc(mp);

            % falco_est_scc() obtains the estimates at all wavelengths, so
            % unpack here.
            for iMode = 1:mp.jac.Nmode
                Ediff = (evTempPos.Eest(:, iMode) - evTempNeg.Eest(:, iMode))/(2*mp.scc.modeCoef);
                jac(:, ibm, iMode) = Ediff;
            end

        end

         jac = mp.dm1.weight * jac;
   
    end
    
    if whichDM == 2
        
        jac = zeros(mp.Fend.corr.Npix, mp.dm2.Nele, mp.jac.Nmode); %--Initialize the Jacobian
        V0 = mp.dm2.V;

        for ibm = 1:mp.dm2.NbasisModes

            fprintf('Empirically obtaining DM2 response matrix for basis mode: %d/%d\n', ibm, mp.dm2.NbasisModes);
            dV = mp.scc.modeCoef * mp.dm2.basisCube(:, :, ibm);
            mp.dm2.V = V0 + dV;
            evTempPos = falco_est_scc(mp);
            mp.dm2.V = V0 - dV;
            evTempNeg = falco_est_scc(mp);

            % falco_est_scc() obtains the estimates at all wavelengths, so
            % unpack here.
            for iMode = 1:mp.jac.Nmode
                Ediff = (evTempPos.Eest(:, iMode) - evTempNeg.Eest(:, iMode))/(2*mp.scc.modeCoef);
                jac(:, ibm, iMode) = Ediff;
            end

        end

         jac = mp.dm2.weight * jac;
   
    end
    

end
