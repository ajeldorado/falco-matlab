% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute a pair-wise probe shape for batch process estimation
% of the electric field in the final focal plane. The rectangular dark
% hole region is specified by its position on one half of the focal plane.
%
% INPUTS
% ------
% mp : structure of model parameters
% InormDes: desired normalized intensity of the probes in the image
% phaseShift: phase shift of the cosine wave [radians]
% starIndex : index of star, in range from 1 to mp.compact.star.count
% 
% INPUTS IN mp
% ------------
% mp.lambda0 : central wavelength of bandpass [meters]
% mp.est.probe.whichDM : Which DM # to use for probing. 1 or 2. Default is 1
% mp.est.probe.radius : Max x/y extent of probed region [actuators].
%    [Equivalent to lambda_central/D if DM is fully illuminated]
% mp.est.probe.xOffset : offset of probe center in x at DM [actuators]. Use
%     to avoid central obscurations.
% mp.est.probe.yOffset : offset of probe center in y at DM [actuators]. Use
%     to avoid central obscurations.
% mp.est.probe.xiOffset : xi (horizontal) offset of probed region's center
% in focal plane. Units of lambda/D.
% mp.est.probe.etaOffset : eta (horizontal) offset of probed region's
% center in focal plane. Units of lambda/D.
% mp.probe.est.width : Width of probed rectangular region. Units of lambda/D.
% mp.probe.est.height : Height of probed rectangular region. Units of lambda/D. 
%
% RETURNS
% -------
% probeCmd: Nact x Nact array of DM actuator commands to make a probe [volts]

function probeCmd = falco_gen_pairwise_probe(mp, InormDes, phaseShift, starIndex)

if mp.est.probe.xiOffset(starIndex) == 0 && mp.est.probe.etaOffset(starIndex) == 0
    error("Probed region's center must be offset from the star location.")
end

% Number of actuators across DM surface (independent of beam for time being)
if mp.est.probe.whichDM == 1
    dm = mp.dm1;
elseif mp.est.probe.whichDM == 2
    dm = mp.dm2;
end
Nact = dm.Nact;
NactPerBeam = mp.P2.D / dm.dm_spacing;

% Coordinates in actuator space
xs = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.xOffset)/Nact;
ys = (-(Nact-1)/2:(Nact-1)/2)/Nact - round(mp.est.probe.yOffset)/Nact;
[XS, YS] = meshgrid(xs, ys);

% Convert units from lambda/D to actuators
lamDIntoAct = Nact / NactPerBeam;
xiOffset = mp.est.probe.xiOffset(starIndex) * lamDIntoAct; 
etaOffset = mp.est.probe.etaOffset(starIndex) * lamDIntoAct;
width = mp.est.probe.width(starIndex) * lamDIntoAct;
height = mp.est.probe.height(starIndex) * lamDIntoAct;

maxSpatialFreq = Nact / 2;
if (xiOffset + width/2) > maxSpatialFreq || (etaOffset + height/2) > maxSpatialFreq
    disp('*** WARNING: SPECIFIED PROBING REGION IN DARK HOLE IS NOT FULLY CONTROLLABLE. ***')
end

%--Generate the DM command for the probe
surfMax = 4*pi*mp.lambda0*sqrt(InormDes); % [meters]
probeHeight = surfMax * sinc(width*XS) .* sinc(height*YS) .* cos(2*pi*(xiOffset*XS + etaOffset*YS) + phaseShift);
probeCmd = falco_fit_dm_surf(dm, probeHeight);
probeCmd = mp.est.probe.gainFudge(starIndex) * probeCmd; % Scale the probe amplitude empirically if needed

end %--END OF FUNCTION