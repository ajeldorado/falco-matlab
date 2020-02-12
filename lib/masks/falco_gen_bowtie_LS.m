% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate a pupil mask with a bow-tie opening.
%
% Modified on 2019-10-03 by A.J. Riggs to allow shear, clocking and
% magnification.
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Modified by A.J. Riggs on 2017-02-09 to be a function with specifyable
%   amounts of magnification, clocking, and translation. 
%   Now keeping the pupil edges non-binary.
% This script was first written by A.J. Riggs on 2017-01-11.
%
% Required Inputs:
% inputs.Dbeam: diameter of the incoming beam
% inputs.Nbeam: number of points across the incoming beam 
% inputs.ID: inner diameter of mask (in pupil diameters)
% inputs.OD: outer diameter of mask (in pupil diameters)
% inputs.ang: opening angle of the bowtie wedges (degrees)
% inputs.centering: centering of the pupil on the array ('pixel' or 'interpixel')

function mask = falco_gen_bowtie_LS(inputs)

% % %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear all
% inputs.Nbeam = 100; 
% inputs.ID = 0.25; % (pupil diameters)
% inputs.OD = 0.82; % (pupil diameters)
% inputs.ang = 115; % (degrees)
% inputs.centering = 'pixel'; % 'interpixel' or 'pixel'
% addpath(genpath('~/Repos/falco-matlab/'))

Nbeam   = inputs.Nbeam;     % number of points across the incoming beam           
ID = inputs.ID; % inner diameter of mask (in pupil diameters)
OD = inputs.OD; % outer diameter of mask (in pupil diameters)
ang = inputs.ang; %opening angle of the upper and lower bowtie wedges (degrees)

Dbeam = 1; %inputs.Dbeam; %--diameter of the beam at the mask (pupil diameters)
dx = Dbeam/Nbeam;
Dmask = Dbeam; % % width of the beam (so can have zero padding if LS is undersized) (meters)

%--Defaults for optional values
centering = 'pixel';  %--Default to pixel centering
xShear = 0; %--x-axis shear of mask [pupil diameters]
yShear = 0; %--y-axis shear of mask [pupil diameters]
clocking = 0; %--Clocking of the mask [degrees]
magfac = 1; %--magnification factor of the pupil diameter

%--Optional inputs
if(isfield(inputs,'centering')); centering = inputs.centering; end
if(isfield(inputs,'xShear')); xShear = inputs.xShear; end
if(isfield(inputs,'yShear')); yShear = inputs.yShear; end
if(isfield(inputs,'clocking')); clocking = inputs.clocking; end
if(isfield(inputs,'magfac')); magfac = inputs.magfac; end

if(strcmpi(centering,'pixel'))
    Narray = ceil_even(magfac*Nbeam+1 + 2*max(Nbeam*abs([xShear,yShear]))); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(magfac*Nbeam + 2*max(Nbeam*abs([xShear,yShear]))); %--number of points across output array. Same size as width when interpixel centered.
end

Darray = Narray*dx; %--width of the output array (meters)
bdf = Dmask/Darray; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--INITIALIZE PROPER
bm = prop_begin(Dmask, wl_dummy, Narray,'beam_diam_fraction',bdf);

%--PRIMARY MIRROR (OUTER DIAMETER)
ra_OD = (Dbeam*OD/2)*magfac;
cx_OD = 0 + cshift + xShear;
cy_OD = 0 + cshift + yShear;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);

%--SECONDARY MIRROR (INNER DIAMETER)
ra_ID = (Dbeam*ID/2)*magfac; %((2.46/2)*Dconv + pad_idod)*magfac;
cx_ID = 0 + cshift + xShear;
cy_ID = 0 + cshift + yShear;
bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);

mask = ifftshift(abs(bm.wf));

%--Create the bowtie region
if(ang<180)
    ang2 = 90-ang/2;
    
    %--Right part
    Lside = 1.1*ra_OD;
    xvert0 = [0, Lside*cosd(ang2), Lside,             Lside,             Lside*cosd(ang2), 0];
    yvert0 = [0, Lside*sind(ang2), Lside*sind(ang2), -Lside*sind(ang2), -Lside*sind(ang2), 0];
    xvert = xvert0;
    yvert = yvert0;
    for ii = 1:length(xvert0)
        xy = [cosd(clocking),sind(clocking); -sind(clocking),cosd(clocking)]*[xvert0(ii); yvert0(ii)];
        xvert(ii) = xy(1);
        yvert(ii) = xy(2);
    end
    bowtieRight = prop_irregular_polygon( bm, cshift+xShear+xvert, cshift+yShear+yvert,'DARK');
    
    %--Left part
    xvert0 = -[0, Lside*cosd(ang2), Lside,             Lside,             Lside*cosd(ang2), 0];
    xvert = xvert0;
    yvert = yvert0;
    for ii = 1:length(xvert0)
        xy = [cosd(clocking),sind(clocking); -sind(clocking),cosd(clocking)]*[xvert0(ii); yvert0(ii)];
        xvert(ii) = xy(1);
        yvert(ii) = xy(2);
    end
    bowtieLeft = prop_irregular_polygon( bm, cshift+xShear+xvert, cshift+yShear+yvert,'DARK');
    
    mask = mask.*bowtieRight.*bowtieLeft;
end

end %--END OF FUNCTION