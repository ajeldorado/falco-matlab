% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Formerly MakeSincSincSinOffset.m
%
function ProbeSurf = fpwc_makeProbeSurf(ProbeArea,D,lambda,psi,offsetX,offsetY,XS,YS,DesiredCont)

mx = (ProbeArea(2)-ProbeArea(1))/D;
my = (ProbeArea(4)-ProbeArea(3))/D;
wx = (ProbeArea(2)+ProbeArea(1))/2;
wy = (ProbeArea(4)+ProbeArea(3))/2;
SincAmp = lambda*sqrt(DesiredCont)*2.51;  % sqrt(2*pi) = 2.51 % Amplitude is in meters!!!
ProbeSurf = SincAmp*sinc(mx*(XS+offsetX)).*sinc(my*(YS+offsetY)).*cos(2*pi*wx/D*XS+ psi).*cos(2*pi*wy/D*YS);

% Offsets move the main lobe away from the center of the PSF. This is
% crucial for AFTA and other telescopes with large secondaries obscuring
% the center of the pupil.
