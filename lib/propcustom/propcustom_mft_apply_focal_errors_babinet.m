% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
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
% resErrorMap : pixels per lambda/D
% 

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
    
    % Backwards MFT to the original pupil plane
    deltaE = propcustom_mft_FtoP(Efocus .* (1-EerrorMap), -fl, lambda, dxi, deta, dx, NarrayPupil);
    
    Epupil = Epupil - deltaE;

end