% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to plot the coronagraphic PSF

function fpwc_imagescPSF(plotTitle,xis,etas,data,c_range)

    imagesc(xis,etas,data,c_range); 
    colorbar; colormap parula;
    title(plotTitle,'FontSize',24,'Interpreter','LaTeX');
    xlabel('$\lambda_0$/D','FontSize',18,'Interpreter','LaTeX'); 
    ylabel('$\lambda_0$/D','FontSize',18,'Interpreter','LaTeX');
    axis xy equal tight;
    set(gca,'FontSize',18,'FontName','Times','FontWeight','Normal')
