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
% rotation: CCW rotation angle of the probe shape at the DM [degrees]
%
% RETURNS
% -------
% probeCmd: Nact x Nact array of DM actuator commands to make a probe [volts]

function probeCmd = falco_gen_pairwise_probe(mp, InormDes, phaseShift, starIndex, rotation)

% Input checks
if ~isa(mp.est.probe, 'Probe')
    error('mp.est.probe must be an instance of class Probe')
end

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
xs = (-(Nact-1)/2:(Nact-1)/2)/Nact - mp.est.probe.xOffset/Nact;
ys = (-(Nact-1)/2:(Nact-1)/2)/Nact - mp.est.probe.yOffset/Nact;
[XS, YS] = meshgrid(xs, ys);

%--Rotate the coordinates
if rotation ~= 0
    RS = sqrt(XS.^2 + YS.^2);
    THETAS = atan2(YS, XS);
    XS = RS.*cos(THETAS - rotation*(pi/180));
    YS = RS.*sin(THETAS - rotation*(pi/180));
end

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

end
