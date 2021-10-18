% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute a pair-wise probe shape for batch process estimation
% of the electric field in the final focal plane. The probed region is a
% square centered on the star.
%
%--INPUTS
%  psi: phase shift of the sine wave [radians]
%  InormDes: desired normalized intensity of the probes in the image
%  mp: structure of model parameters
%    -mp.est.probe.whichDM:      Which DM # to use for probing. 1 or 2. Default is 1
%    -mp.est.probe.radius:  Max horizontal (or vertical) extent of probed square region [lambda/D].
%    [Equivalent to lambda_central/D if DM is fully illuminated]
%    -mp.est.probe.xOffset: offset of probe center in x [actuators]. Use
%     to avoid central obscurations.
%    -mp.est.probe.yOffset: offset of probe center in y [actuators]. Use
%     to avoid central obscurations.
%    -mp.lambda0:       central wavelength of bandpass [meters]
%    -mp.est.probe.axis:    which axis to have the phase discontinuity
%    along [x or y]
% % % % %  Nact: number of actuators across the (square-grid) actuator array
%
%--OUTPUTS
%  probeCmd: Nact x Nact array of DM actuator commands to make a probe

function probeCmd = falco_gen_pairwise_probe_square(mp, InormDes, psi, badAxis)

% Number of actuators across DM surface (independent of beam for time being)
if mp.est.probe.whichDM == 1
    dm = mp.dm1;
elseif mp.est.probe.whichDM == 2
    dm = mp.dm2;
end
Nact = dm.Nact;
NactPerBeam = mp.P2.D / dm.dm_spacing;

%--Coordinates in actuator space
xs = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.xOffset)/Nact;
ys = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.yOffset)/Nact;
[XS,YS] = meshgrid(xs,ys);

% Probed region extend in dark hole
lamDIntoAct = Nact / NactPerBeam; % Convert units from lambda/D to actuators
probeRadius = mp.est.probe.radius * lamDIntoAct;  
if(probeRadius > Nact/2)
    probeRadius = Nact/2;
end

%--Generate the DM command for the probe
surfMax = 4*pi*mp.lambda0*sqrt(InormDes); % peak surface height to get desired intensity [meters]
switch lower(badAxis)
    case 'y'
        mX = probeRadius;
        mY = 2*probeRadius;
        omegaX = probeRadius/2;        
        probeCmd = surfMax*sinc(mX*XS).*sinc(mY*YS).*cos(2*pi*omegaX*XS + psi);

    case 'x'
        mX = 2*probeRadius;
        mY = probeRadius;
        omegaY = probeRadius/2;
        probeCmd = surfMax*sinc(mX*XS).*sinc(mY*YS).*cos(2*pi*omegaY*YS + psi);

    case 'm'
        omegaX = mp.est.probe.Xloc/2;
        omegaY = mp.est.probe.Yloc/2;
        probeCmd = zeros(size(XS));
        for iFiber = 1:mp.Fend.Nfiber
            probeCmd = probeCmd + surfMax*sin(2*pi*omegaX(iFiber)*XS + 2*pi*omegaY(iFiber)*YS + psi);
        end
end

%--Option to use just the sincs for a zero phase shift. This avoids the
% phase discontinuity along one axis (for this probe only!).
if(psi==0 && ~mp.flagFiber)
    m = 2*probeRadius;
    probeCmd = surfMax * sinc(m*XS) .* sinc(m*YS);
end 

probeCmd = falco_fit_dm_surf(dm, probeCmd);

probeCmd = mp.est.probe.gainFudge * probeCmd; % Scale the probe amplitude empirically if needed

end %--END OF FUNCTION