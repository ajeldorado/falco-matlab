% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%
% Plot the singular mode spectrum of the E-field.

function out = falco_plot_singular_mode_spectrum_of_Efield(mp, out, jacStruct, Eest, Itr)
    
    % Initialize Gall and Eall
    iMode = 1;
    Gcomplex = [jacStruct.G1(:,:,iMode), jacStruct.G2(:,:,iMode), jacStruct.G3(:,:,iMode), jacStruct.G4(:,:,iMode), jacStruct.G5(:,:,iMode), jacStruct.G6(:,:,iMode), jacStruct.G7(:,:,iMode), jacStruct.G8(:,:,iMode), jacStruct.G9(:,:,iMode)];
    Gall = zeros(mp.jac.Nmode*size(Gcomplex, 1), size(Gcomplex, 2));
    Eall = zeros(mp.jac.Nmode*size(Eest, 1), 1);

    for iMode = 1:mp.jac.Nmode
        N = size(Gcomplex, 1);
        inds = (iMode-1)*N+1 : iMode*N;
        Gcomplex = [jacStruct.G1(:,:,iMode), jacStruct.G2(:,:,iMode), jacStruct.G3(:,:,iMode), jacStruct.G4(:,:,iMode), jacStruct.G5(:,:,iMode), jacStruct.G6(:,:,iMode), jacStruct.G7(:,:,iMode), jacStruct.G8(:,:,iMode), jacStruct.G9(:,:,iMode)];
        Gall(inds, :) = Gcomplex;
        Eall(inds) = Eest(:, iMode);
    end

    Eri = [real(Eall); imag(Eall)];
    alpha2 = max(diag(real(Gall' * Gall)));
    Gri = [real(Gall); imag(Gall)];
    [U, S, ~] = svd(Gri, 'econ');
    s = diag(S);

    EriPrime = U' * Eri;
    IriPrime = abs(EriPrime).^2;

    % Store data for later analysis
    out.EforSpectra{Itr} = EriPrime;
    out.smspectra{Itr} = IriPrime;
    out.sm{Itr} = s;
    out.alpha2{Itr} = alpha2;

    if(mp.flagPlot)
        figure(401);
        if(Itr == 1); hold off; end
        loglog(out.sm{Itr}.^2/out.alpha2{Itr},smooth(out.smspectra{Itr},31), 'Linewidth', 2,...
            'Color', [0.3, 1-(0.2+Itr/mp.Nitr)/(1.3), 1]);
        grid on; set(gca,'minorgridlines','none')
        set(gcf,'Color',[1 1 1]);
        title('Singular Mode Spectrum')%,'Fontsize',20)
        xlim([1e-10, 2*max(s.^2/alpha2)])
        ylim([1e-12, 1e-0]) 
        drawnow;
        hold on;
    end

    clear Gcomplex Gall Eall Eri alpha2 Uri U S EriPrime IriPrime
        
end %--END OF FUNCTION