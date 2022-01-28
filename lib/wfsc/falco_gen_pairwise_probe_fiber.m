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
%    -mp.est.probe.xOffset: offset of probe center in x [actuators]. Use to avoid central obscurations.
%    -mp.est.probe.yOffset: offset of probe center in y [actuators]. Use to avoid central obscurations.
%    -mp.lambda0:       central wavelength of bandpass [meters]
%    -mp.est.probe.axis:    which axis to have the phase discontinuity along [x or y]

function probeCmd = falco_gen_pairwise_probe_fiber(mp,InormDes,psi,badAxis)

%--Number of actuators across DM surface (independent of beam for time being)
if(mp.est.probe.whichDM==1)
    Nact = mp.dm1.Nact;
    dm = mp.dm1;
elseif(mp.est.probe.whichDM==2)
    Nact = mp.dm2.Nact;
    dm = mp.dm2;
end
%--Coordinates in actuator space
xs = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.xOffset)/Nact;
ys = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.yOffset)/Nact;
[XS,YS] = meshgrid(xs,ys);

%--Restrict the probing region if it is not possible to achieve
if(mp.est.probe.radius > Nact/2)
    mp.est.probe.radius = Nact/2;
end

%--Generate the DM command for the probe
magn = 4*pi*mp.lambda0*sqrt(InormDes);   % surface height to get desired intensity [meters]
switch lower(badAxis)
    case 'y'
        omegaX = mp.est.probe.Xloc(1)/2;
        probeCmd = magn*sin(2*pi*omegaX*XS + psi);

    case 'x'
        omegaY = mp.est.probe.Yloc(1)/2;
        probeCmd = magn*sin(2*pi*omegaY*YS + psi);
        
    case 'm'
        omegaX = mp.est.probe.Xloc/2;
        omegaY = mp.est.probe.Yloc/2;
        probeCmd = zeros(size(XS));
        for i = 1:mp.Fend.Nfiber
            probeCmd = probeCmd + magn*sin(2*pi*omegaX(i)*XS + psi).*sin(2*pi*omegaY(i)*YS + psi);
        end
end

probeCmd = falco_fit_dm_surf(dm,probeCmd);

probeCmd = mp.est.probe.gainFudge*probeCmd; % Scale the probe amplitude empirically if needed

end %--END OF FUNCTION
