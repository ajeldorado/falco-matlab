% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%  Empirically compute the Jacobian for implicit EFC (IEFC).
%
% INPUTS
% ------
% mp : structure of model parameters
% whichDM : which DM number
%
% OUTPUTS
% -------
% jac : Jacobian for the specified DM.

function jac = model_Jacobian_IEFC(mp, whichDM)

    if ~mp.flagSim % jllopsay
        mp.tint = mp.tint_est;
    end
    flagPlot = true; % plotting for debugging
    Nprobes = size(mp.iefc.probeCube, 3);
    
    if whichDM == 1
        
        V0 = mp.dm1.V;
        jac = zeros(Nprobes*mp.Fend.corr.Npix, mp.dm1.Nele, mp.jac.Nmode); %--Initialize the Jacobian

        for iJacMode = 1:mp.jac.Nmode

            for iBasisMode = 1:mp.dm1.NbasisModes

                fprintf('Empirically obtaining DM1 response matrix for basis mode: %d/%d\n', iBasisMode, mp.dm1.NbasisModes);
                dVbasis = mp.iefc.modeCoef * mp.dm1.basisCube(:, :, iBasisMode);

                for iProbe = 1:Nprobes

                    dVprobe = mp.iefc.probeCoef * mp.iefc.probeCube(:, :, iProbe);
 
                    mp.dm1.V = V0 + dVbasis;
                    plusStruct = falco_est_delta_intensity(mp, mp.iefc.probeDM, dVprobe, iJacMode);
                    mp.dm1.V = V0 - dVbasis;
                    minusStruct = falco_est_delta_intensity(mp, mp.iefc.probeDM, dVprobe, iJacMode);
                    mp.dm1.V = V0; % reset

                    jac((iProbe-1)*mp.Fend.corr.Npix+1:iProbe*mp.Fend.corr.Npix, iBasisMode, iJacMode) = (plusStruct.DeltaI - minusStruct.DeltaI)/(2*mp.iefc.modeCoef*mp.iefc.probeCoef);

                    if flagPlot
                        figure(30); imagesc(dVbasis); axis xy equal tight; colorbar; 
                        title(sprintf('Basis Mode %d/%d', iBasisMode, mp.dm1.NbasisModes));
                        drawnow;

                        E2D = zeros(mp.Fend.Nxi, mp.Fend.Neta);
                        E2D(mp.Fend.corr.maskBool) = jac(1:mp.Fend.corr.Npix, iBasisMode, iJacMode);
                        figure(31); imagesc(log10(abs(E2D))); axis xy equal tight; colorbar; %,[-1.5 0.5]
                        
                        title(sprintf('Basis Mode %d/%d', iBasisMode, mp.dm1.NbasisModes));
                        drawnow;
                    end

                end  

            end
        
        end
        
        jac = mp.dm1.weight * jac;

        
    elseif whichDM == 2
        
        V0 = mp.dm2.V;
        jac = zeros(Nprobes*mp.Fend.corr.Npix, mp.dm2.Nele, mp.jac.Nmode); %--Initialize the Jacobian

        for iJacMode = 1:mp.jac.Nmode

            for iBasisMode = 1:mp.dm2.NbasisModes

                fprintf('Empirically obtaining DM2 response matrix for basis mode: %d/%d\n', iBasisMode, mp.dm2.NbasisModes);
                dVbasis = mp.iefc.modeCoef * mp.dm2.basisCube(:, :, iBasisMode);

                for iProbe = 1:Nprobes

                    dVprobe = mp.iefc.probeCoef * mp.iefc.probeCube(:, :, iProbe);
 
                    mp.dm2.V = V0 + dVbasis;
                    plusStruct = falco_est_delta_intensity_subband(mp, mp.iefc.probeDM, dVprobe, iJacMode);
                    mp.dm2.V = V0 - dVbasis;
                    minusStruct = falco_est_delta_intensity_subband(mp, mp.iefc.probeDM, dVprobe, iJacMode);
                    mp.dm2.V = V0; % reset

                    jac((iProbe-1)*mp.Fend.corr.Npix+1:iProbe*mp.Fend.corr.Npix, iBasisMode, iJacMode) = (plusStruct.DeltaI - minusStruct.DeltaI)/(2*mp.iefc.modeCoef*mp.iefc.probeCoef);

                    if flagPlot
                        figure(30); imagesc(dVmode); axis xy equal tight; colorbar; 
                        title(sprintf('Basis Mode %d/%d', iBasisMode, Nbasis));
                        drawnow;

                        E2D = zeros(mp.Fend.Nxi, mp.Fend.Neta);
                        E2D(mp.Fend.corr.maskBool) = jac(1:mp.Fend.corr.Npix, iBasisMode, iJacMode);
                        figure(31); imagesc(log10(abs(E2D))); axis xy equal tight; colorbar; 
                        title(sprintf('Basis Mode %d/%d', iBasisMode, Nbasis));
                        drawnow;
                    end

                end  

            end
        
        end
        
        jac = mp.dm2.weight * jac;
    
    end

     

end
