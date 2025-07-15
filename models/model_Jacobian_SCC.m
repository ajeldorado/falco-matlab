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
    
    mp.flagPlot = false;
    flagPlot = true; % plotting for debugging
    if whichDM == 1
        
        jac = zeros(mp.Fend.corr.Npix, mp.dm1.Nele, mp.jac.Nmode); %--Initialize the Jacobian
        V0 = mp.dm1.V;

        for ibm = 1:mp.dm1.NbasisModes

            fprintf('Empirically obtaining DM1 response matrix for basis mode: %d/%d\n', ibm, mp.dm1.NbasisModes);
            dV = mp.scc.modeCoef * mp.dm1.basisCube(:, :, ibm);
            mp.scc.order = 'open_close';
            mp.dm1.V = V0 + dV;
            evTempPos = falco_est_scc(mp);
            mp.scc.order = 'close_open';
            mp.dm1.V = V0 - dV;
            evTempNeg = falco_est_scc(mp);
            

            % falco_est_scc() obtains the estimates at all wavelengths, so
            % unpack here.
            for iMode = 1:mp.jac.Nmode
                Ediff = (evTempPos.Eest_full(:, :, 1, iMode) - evTempNeg.Eest_full(:, :, 1, iMode))/(2*mp.scc.modeCoef);
                jac(:, ibm, iMode) = Ediff(mp.Fend.corr.maskBool);
                
                if flagPlot
                    Diff_im = Ediff.*mp.Fend.corr.mask;
                    figure(1020); imagesc(real(Diff_im)); axis xy equal tight; colorbar;
                    title(sprintf('DM1 response matrix for basis mode: %d/%d\n', ibm, mp.dm1.NbasisModes))
                    drawnow;
                end
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
            mp.scc.order = 'open_close';
            mp.dm1.V = V0 + dV;
            evTempPos = falco_est_scc(mp);
            mp.scc.order = 'close_open';
            mp.dm1.V = V0 - dV;
            evTempNeg = falco_est_scc(mp);
            

            % falco_est_scc() obtains the estimates at all wavelengths, so
            % unpack here.
            for iMode = 1:mp.jac.Nmode
                Ediff = (evTempPos.Eest_full(:, :, 1, iMode) - evTempNeg.Eest_full(:, :, 1, iMode))/(2*mp.scc.modeCoef);
                jac(:, ibm, iMode) = Ediff(mp.Fend.corr.maskBool);
                
                if flagPlot
                    Diff_im = Ediff.*mp.Fend.corr.mask;
                    figure(1020); imagesc(real(Diff_im)); axis xy equal tight; colorbar;
                    title(sprintf('DM2 response matrix for basis mode: %d/%d\n', ibm, mp.dm2.NbasisModes))
                    drawnow;
                end
            end

        end

         jac = mp.dm2.weight * jac;
   
    end
    
    mp.scc.order = 'open_close';
    flagPlot = false; % plotting for debugging
    mp.flagPlot = true;
end
