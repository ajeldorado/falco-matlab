% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the Cycle 6 WFIRST CGI pupil (on-axis, though) in
% Matlab using PROPER. This version allows padding of the ID, OD, and
% struts to use the mask as an under-sized Lyot stop.
% Coordinates and dimensions of the struts, primary, and secondary are from
% a Visio file from Kent Wallace of the March 2, 2016 Cycle 6 WFIRST pupil.
% The pupil used is on-axis.
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.
% Modified by A.J. Riggs on 2017-02-09 to be a function with specifyable
%   amounts of magnification, clocking, and translation. 
%   Now keeping the pupil edges non-binary.
% This script was first written by A.J. Riggs on 2017-01-11.
%
%
%--Inputs:
% inputs.Dbeam: diameter of the incoming beam
% inputs.Nbeam: number of points across the incoming beam 
% inputs.ID: inner diameter of mask (in pupil diameters)
% inputs.OD: outer diameter of mask (in pupil diameters)
% inputs.ang: opening angle of the bowtie wedges (degrees)
% inputs.centering: centering of the pupil on the array ('pixel' or 'interpixel')
%
%--To allow magnification, clocking, and translation:
% Center of each feature needs to be magnified and clocked first.
%   Then translate the centers.
% Circles can be magnified and translated in either order.
% Strut width also needs to scale with magnification

function mask = falco_gen_bowtie_LS(inputs,varargin)

% Set default values of input parameters
flagRot180deg = false;
%--Look for Optional Keywords
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'rot180'}
            flagRot180deg = true; % For even arrays, beam center is in between pixels.
        otherwise
            error('falco_gen_bowtie_LS: Unknown keyword: %s\n', ...
            varargin{icav});
    end
end

% % %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear all
% inputs.Dbeam = 48e-3; % meters;
% inputs.Nbeam = 100; 
% inputs.ID = 0.25; % (pupil diameters)
% inputs.OD = 0.82; % (pupil diameters)
% inputs.ang = 115; % (degrees)
% inputs.centering = 'pixel'; % 'interpixel' or 'pixel'
% addpath(genpath('~/Repos/FALCO/'))
% flagRot180deg = false;




Dbeam = inputs.Dbeam; %--diameter of the beam at the mask (meters)
Nbeam   = inputs.Nbeam;     % number of points across the incoming beam           
ID = inputs.ID; % inner diameter of mask (in pupil diameters)
OD = inputs.OD; % outer diameter of mask (in pupil diameters)
ang = inputs.ang; %opening angle of the upper and lower bowtie wedges (degrees)

dx = Dbeam/Nbeam;
Dmask = Dbeam; % % width of the beam (so can have zero padding if LS is undersized) (meters)

% Centering of array: 'pixel' or 'interpixel'
if(isfield(inputs,'centering'))
    centering = inputs.centering;
else %--Default to pixel centering
    centering = 'pixel';
end

%--Other INPUTS
clock_deg = 0;%inputs.clock_deg; % clocking angle of the pupil (in degrees)
magfacD = 1;%inputs.magfacD; %magnification factor of the pupil diameter
xshift = 0;%inputs.xshift; % translation in x of pupil (in diameters)
yshift = 0;%inputs.yshift; % translation in y of pupil (in diameters)
pad_strut = 0; %2*pad_strut_pct/100*diam; %--Convert to meters. Factor of 2x needed since strut is padded on both sides

if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam+1/2); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even(Nbeam); %--number of points across output array. Same size as width when interpixel centered.
end

Darray = Narray*dx; %--width of the output array (meters)
bdf = Dmask/Darray; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
        if(flagRot180deg)
            cshift = -dx;
        end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%--INITIALIZE PROPER
bm = prop_begin(Dmask, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

%--PRIMARY MIRROR (OUTER DIAMETER)
% ra_OD = (diam/2 - pad_idod)*magfacD;
ra_OD = (Dbeam*OD/2)*magfacD;
cx_OD = 0 + cshift + xshift;
cy_OD = 0 + cshift + yshift;
% norm = false;
bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

%--SECONDARY MIRROR (INNER DIAMETER)
ra_ID = (Dbeam*ID/2)*magfacD; %((2.46/2)*Dconv + pad_idod)*magfacD;
cx_ID = 0 + cshift + xshift;
cy_ID = 0 + cshift + yshift;
bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);%, cx, cy, norm)
% figure(4); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;
% figure(104); imagesc(ifftshift(abs(bm.wf))-rot90(ifftshift(abs(bm.wf)),2)); axis xy equal tight; colorbar;

mask = ifftshift(abs(bm.wf));

%--Create the bowtie region
if(ang<180)
    ang2 = 90-ang/2;
    
    %--Right part
    Lside = 1.1*ra_OD;
    xvert = cshift+[0, Lside*cosd(ang2), Lside,             Lside,             Lside*cosd(ang2), 0];
    yvert = cshift+[0, Lside*sind(ang2), Lside*sind(ang2), -Lside*sind(ang2), -Lside*sind(ang2), 0];     
    %figure(30); plot(xvert,yvert,'-bo');
    bowtieRight = prop_irregular_polygon( bm, xvert, yvert,'DARK');
    % figure(5); imagesc((abs(bowtieRight))); axis xy equal tight; colorbar

    
    %--Left part
    xvert = cshift - [0, Lside*cosd(ang2), Lside,             Lside,             Lside*cosd(ang2), 0];
    bowtieLeft = prop_irregular_polygon( bm, xvert, yvert,'DARK');
    %figure(6); imagesc((abs(bowtieLeft))); axis xy equal tight; colorbar
    
    mask = mask.*bowtieRight.*bowtieLeft;
end

if(flagRot180deg)
    mask = rot90(mask,2);
end
% figure(18); imagesc(mask); axis xy equal tight; colorbar;

end %--END OF FUNCTION


% %--DEBUGGING: Visually verify that mask is centered correctly
% figure(11); imagesc(mask); axis xy equal tight; colorbar; drawnow;
% switch centering 
%     case {'pixel'}
%         maskTemp = mask(2:end,2:end);
%     otherwise
%         maskTemp = mask;
% end
% figure(12); imagesc(maskTemp-rot90(maskTemp,2)); axis xy equal tight; colorbar; 
% title('Centering Check','Fontsize',20); set(gca,'Fontsize',20);
% drawnow;
% 
% sum(sum(mask))












