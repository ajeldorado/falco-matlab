% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute a pair-wise probe shape for batch process estimation
%  of the electric field in the final focal plane.
%
%--INPUTS
%  ProbeArea:
%  Nact: number of actuators across the (square-grid) actuator array
%  psi: phase shift of the sincs
%  offsetX: offset of probe center in x (actuators)
%  offsetY: offset of probe center in y (actuators)
%  lambda: center wavelength
%  InormDes: desired normalized intensity of the probes in the image
%
%--OUTPUTS
%  probeCmd: Nact x Nact array of DM actuator commands to make a probe
%
%--VERSION CHANGE HISTORY
% Modified on 2018-04-23 by A.J. Riggs at JPL.
% Written on 2014-02-18 by A.J. Riggs at Princeton University.

function probeCmd = falco_gen_probe_command(ProbeArea,Nact,psi,offsetX,offsetY,lambda,InormDes)
% function probeCmd = falco_gen_probe_command(ProbeArea,D,lambda,psi,offsetX,offsetY,XS,YS,DesiredCont)

% %--Hard-coded values for debugging
% xs = (-(Nact-1)/2:(Nact-1)/2)/Nact;
% D = 1;
% psi = pi/2;
% Nact = 48;
% offsetX = 10; %--offset in y (actuators)
% offsetY = 10; %--offset of probe center in y (actuators)
% ProbeArea = [0 Nact/2,-Nact/2, Nact/2]; % [x_inner, x_outer, y_lower, y_upper

xs = (-(Nact-1)/2:(Nact-1)/2)/Nact;
[XS,YS] = meshgrid(xs);

mx = (ProbeArea(2)-ProbeArea(1))/D;
my = (ProbeArea(4)-ProbeArea(3))/D;
wx = (ProbeArea(2)+ProbeArea(1))/2;
wy = (ProbeArea(4)+ProbeArea(3))/2;

magn = lambda*sqrt(2*pi)*sqrt(InormDes);   % surface height (meters) to get desired intensity
probeCmd = magn*sinc(mx*(XS-offsetX/Nact)).*sinc(my*(YS-offsetY/Nact)).*cos(2*pi*wx/D*XS+ psi).*cos(2*pi*wy/D*YS);

% figure(1); imagesc(probeCmd); axis xy equal tight; colorbar;






% function ProbeSurf = falco_gen_probe_command(ProbeArea,D,lambda,psi,offsetX,offsetY,XS,YS,DesiredCont)
% 
% mx = (ProbeArea(2)-ProbeArea(1))/D;
% my = (ProbeArea(4)-ProbeArea(3))/D;
% wx = (ProbeArea(2)+ProbeArea(1))/2;
% wy = (ProbeArea(4)+ProbeArea(3))/2;
% sincAmp = lambda*sqrt(2*pi)*sqrt(DesiredCont);   % surface height in meters
% ProbeSurf = sincAmp*sinc(mx*(XS+offsetX)).*sinc(my*(YS+offsetY)).*cos(2*pi*wx/D*XS+ psi).*cos(2*pi*wy/D*YS);

% Offsets move the main lobe away from the center of the PSF. This is
% crucial for AFTA and other telescopes with large secondaries obscuring
% the center of the pupil.


%%





