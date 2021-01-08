% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%
% Plot the model-based (simulated) and estimated change in electric field
% at each subband.

function out = falco_plot_DeltaE(mp, out, EfieldMeas, EfieldMeasPrev, EfieldSim, EfieldSimPrev, Itr)

    if(Itr > 1 && ~any(mp.ctrl.dmfacVec == 0))
        for si = 1:mp.Nsbp
            dEmeas = squeeze(EfieldMeas(:, si) - EfieldMeasPrev(:, si));
            dEmeas2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
            dEmeas2D(mp.Fend.corr.maskBool) = dEmeas; % 2-D for plotting
            
            dEsim = squeeze(EfieldSim(:, si) - EfieldSimPrev(:, si));
            dEsim2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
            dEsim2D(mp.Fend.corr.maskBool) = dEsim;  % 2-D for plotting
            
            dEmax = max(abs(dEsim)); % max value in plots
            out.complexProjection(Itr-1, si) = abs(dEsim'*dEmeas) / abs(dEsim'*dEsim);
            fprintf('Complex projection of deltaE is %3.2f    for subband %d/%d\n', out.complexProjection(Itr-1, si), si, mp.Nsbp);
            out.complexCorrelation(Itr-1, si) = abs(dEsim'*dEmeas/(sqrt(abs(dEmeas'*dEmeas))*sqrt(abs(dEsim'*dEsim))));
            fprintf('Complex correlation of deltaE is %3.2f    for subband %d/%d\n', out.complexCorrelation(Itr-1, si), si, mp.Nsbp);

            if mp.flagPlot            
                figure(50+si); set(gcf, 'Color', 'w');
                fs = 18;

                hModelAmp = subplot(2,2,1); % Save the handle of the subplot
                imagesc(mp.Fend.xisDL, mp.Fend.etasDL, abs(dEsim2D)); axis xy equal tight; colorbar; colormap(hModelAmp, 'parula');
                title('abs(dE_{model})', 'Fontsize', fs); 
                set(gca,'FontSize', fs); %,'FontName','Times','FontWeight','Normal')

                hMeasAmp = subplot(2,2,2); % Save the handle of the subplot
                imagesc(mp.Fend.xisDL, mp.Fend.etasDL, abs(dEmeas2D), [0, dEmax]); axis xy equal tight; colorbar; colormap(hMeasAmp, 'parula');
                title('abs(dE_{meas})', 'Fontsize', fs); 
                set(gca,'FontSize', fs); %,'FontName','Times','FontWeight','Normal')

                hModelPh = subplot(2,2,3); % Save the handle of the subplot
                imagesc(mp.Fend.xisDL, mp.Fend.etasDL, angle(dEsim2D)); axis xy equal tight; colorbar; colormap(hModelPh, 'hsv');
                title('angle(dE_{model})', 'Fontsize', fs); 
                set(gca,'FontSize', fs); %,'FontName','Times','FontWeight','Normal')

                hMeasPh = subplot(2,2,4); % Save the handle of the subplot
                imagesc(mp.Fend.xisDL, mp.Fend.etasDL, angle(dEmeas2D)); axis xy equal tight; colorbar; colormap(hMeasPh, 'hsv');
                title('angle(dE_{meas})', 'Fontsize', fs); 
                set(gca,'FontSize', fs); %,'FontName','Times','FontWeight','Normal')
                drawnow;
            end
        end
    end

end %--END OF FUNCTION