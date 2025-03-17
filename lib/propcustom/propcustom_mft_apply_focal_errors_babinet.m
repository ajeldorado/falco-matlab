% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Apply an error map to a focal plane. Use Babinet's principle to avoid
% windowing the focal plane.
%
% Propates from the pupil before the focal plane, applies the error map,
% then back propagates to the original pupil plane.
%
% INPUTS
% Epupil : 2-D array of the electric field at the pupil from which to start
% EerrorMap : 2-D arry of the complex-valued errors to multiply at the focal plane 
% resErrorMap : Resolution of the error map in pixels per lambda0/D.
% beamDiamPix : Number of pixels across the beam diameter in the pupil.
%
% OUTPUTS
% Epupil : 2-D E-field array back at the original pupil plane but
%   containing the errors as applied at the focal plane

function Epupil = propcustom_mft_apply_focal_errors_babinet(Epupil, EerrorMap, resErrorMap, beamDiamPix)
    
    NarrayPupil = size(Epupil, 1);
    dx = 1/beamDiamPix;

    dxi = 1/resErrorMap;
    deta = dxi;
    Neta = size(EerrorMap, 1);
    Nxi = size(EerrorMap, 2);
    
    % Normalized units since resErrorMap is in pixels per lambda/D.
    fl = 1;
    lambda = 1;

    Efocus = propcustom_mft_PtoF(Epupil, fl, lambda, dx, dxi, Nxi, deta, Neta);
%     figure(222); semilogy(abs(Efocus(Neta/2+1,Neta/2+1:end)).^2);  title('Roddier focal plane: PSF');
%     figure(223); plot(angle(EerrorMap(Neta/2+1,Neta/2+1:end)));  title('Roddier focal plane: mask');
    
    % Backwards MFT to the original pupil plane
    deltaE = propcustom_mft_FtoP(Efocus .* (1-EerrorMap), -fl, lambda, dxi, deta, dx, NarrayPupil);
    
    % figure(14);
    % imagesc(abs(deltaE));
    % axis xy equal tight; colorbar;
    % colormap gray;
    % 
    % figure(15);
    % imagesc(angle(deltaE));
    % axis xy equal tight; colorbar;
    % colormap hsv;



    Epupil = Epupil + deltaE; %This should be a plus sign! (Always check plots to verify)

end
