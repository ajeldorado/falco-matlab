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
%  psi: phase shift of the sine wave [radians]
%  InormDes: desired normalized intensity of the probes in the image
%  mp: structure of model parameters
%    -mp.est.probe.whichDM:      Which DM # to use for probing. 1 or 2. Default is 1
%    -mp.est.probe.radius:  Max x/y extent of probed region [actuators].
%    [Equivalent to lambda_central/D if DM is fully illuminated]
%    -mp.est.probe.offsetX: offset of probe center in x [actuators]. Use
%     to avoid central obscurations.
%    -mp.est.probe.offsetY: offset of probe center in y [actuators]. Use
%     to avoid central obscurations.
%    -mp.lambda0:       central wavelength of bandpass [meters]
%    -mp.est.probe.axis:    which axis to have the phase discontinuity
%    along [x or y]
% % % % %  Nact: number of actuators across the (square-grid) actuator array
%
%--OUTPUTS
%  probeCmd: Nact x Nact array of DM actuator commands to make a probe
%
%--VERSION CHANGE HISTORY
% Modified on 2019-01-23 by A.J. Riggs to use structures and a better fit.
% Modified on 2018-04-23 by A.J. Riggs at JPL.
% Written on 2014-02-18 by A.J. Riggs at Princeton University.
%
%
%--New variables
%    -mp.est.probe.whichDM:      Which DM # to use for probing. 1 or 2. Default is 1
%    -mp.est.probe.radius:  Max x/y extent of probed region [actuators].
%    -mp.est.probe.offsetX: offset of probe center in x [actuators]. Use to avoid central obscurations.
%    -mp.est.probe.offsetY: offset of probe center in y [actuators]. Use to avoid central obscurations.
%    -mp.lambda0:       central wavelength of bandpass [meters]
%    -mp.est.probe.axis:    which axis to have the phase discontinuity along [x or y]


function probeCmd = falco_gen_pairwise_probe(mp,InormDes,psi)
% function probeCmd = falco_gen_pairwise_probe(ProbeArea,Nact,psi,offsetX,offsetY,lambda,InormDes)
% function probeCmd = falco_gen_probe_command(ProbeArea,D,lambda,psi,offsetX,offsetY,XS,YS,DesiredCont)

%--Number of actuators across DM surface (independent of beam for time being)
if(mp.est.probe.whichDM==1)
    Nact = mp.dm1.Nact;
    dm = mp.dm1;
elseif(mp.est.probe.whichDM==2)
    Nact = mp.dm2.Nact;
    dm = mp.dm2;
end

xs = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.offsetX/Nact);
ys = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.offsetY/Nact);
[XS,YS] = meshgrid(xs,ys);

magn = mp.lambda0*sqrt(2*pi)*sqrt(InormDes);   % surface height to get desired intensity [meters]

switch lower(mp.est.probe.axis)
    case 'y'
        mX = mp.est.probe.radius;
        mY = 2*mp.est.probe.radius;
        omegaX = mp.est.probe.radius/2;        
        probeCmd = magn*sinc(mX*XS).*sinc(mY*YS).*cos(2*pi*omegaX*XS + psi);

    case 'x'
        mX = 2*mp.est.probe.radius;
        mY = mp.est.probe.radius;
        omegaY = mp.est.probe.radius/2;
        probeCmd = magn*sinc(mX*XS).*sinc(mY*YS).*cos(2*pi*omegaY*YS + psi);

end

%--Option to use just the sincs for a zero phase shift. This avoids the
% phase discontinuity along one axis.
if(psi==0)
    m = 2*mp.est.probe.radius;
    probeCmd = magn*sinc(m*XS).*sinc(m*YS);

end


probeCmd = falco_fit_dm_surf(dm,probeCmd);
% figure(1); imagesc(probeCmd); axis xy equal tight; colorbar;

end %--END OF FUNCTION




% OLD CODE (for reference only)
%
% function ProbeSurf = falco_gen_probe_command(ProbeArea,D,lambda,psi,offsetX,offsetY,XS,YS,DesiredCont)
%
% %--Hard-coded values for debugging
% xs = (-(Nact-1)/2:(Nact-1)/2)/Nact;
% D = 1;
% psi = pi/2;
% Nact = 48;
% offsetX = 10; %--offset in y (actuators)
% offsetY = 10; %--offset of probe center in y (actuators)
% ProbeArea = [0 Nact/2,-Nact/2, Nact/2]; % [x_inner, x_outer, y_lower, y_upper
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





