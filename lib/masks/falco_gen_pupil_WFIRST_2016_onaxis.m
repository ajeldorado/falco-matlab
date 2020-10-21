% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate the on-axis version of the Cycle 6 WFIRST CGI pupil
% in Matlab using PROPER.
% Coordinates and dimensions of the struts, primary, and secondary are from
% a Visio file from Kent Wallace of the March 2, 2016 Cycle 6 WFIRST pupil.
%
%--Optional Inputs:
% changes.flagLyot: whether the mask is a Lyot stop or not (true/false)
% changes.Nbeam: Number of points across the illuminated pupil diameter
% changes.Narray: Number of points across the output array
% changes.magFac: magnification factor of the pupil diameter
% changes.clock_deg: clocking angle of the pupil (in degrees)
% changes.xShear: x-translation in (original) pupil diameters
% changes.yShear: y-translation in (original) pupil diameters
%
%--To allow magnification, clocking, and translation:
% Center of each feature needs to be magnified and clocked first.
%   Then translate the centers.
% Circles can be magnified and translated in either order.
% Strut width also needs to scale with magnification

function mask = falco_gen_pupil_WFIRST_2016_onaxis(Nbeam, centering, changes)

Dbeam = 1; % [pupil diameters]

if(isfield(changes,'flagRot180')); flagRot180 = changes.flagRot180; else; flagRot180 = false; end
if(isfield(changes,'clock_deg')); clock_deg = changes.clock_deg; else; clock_deg = 0; end
if(isfield(changes,'mag')); magFac = changes.magFac; else; magFac = 1; end
if(isfield(changes,'xShear')); xShear = changes.xShear; else; xShear = 0; end
if(isfield(changes,'yShear')); yShear = changes.yShear; else; yShear = 0; end
if(isfield(changes,'flagLyot')); flagLyot = changes.flagLyot; else; flagLyot = false; end
if(isfield(changes,'pad_strut')); pad_strut = changes.pad_strut; else; pad_strut = 0; end % Factor of 2x needed for both sides
if(isfield(changes,'pad_idod')); pad_idod = changes.pad_idod; else; pad_idod = 0; end

% Optional Lyot stop mode (concentric, circular ID and OD)
if(flagLyot == true)
    if(isfield(changes,'ID'))
        ID = changes.ID;
    else
        error('changes.ID must be defined for Lyot stop generation mode.')
    end

    if(isfield(changes,'OD'))
        OD = changes.OD;
    else
        error('changes.OD must be defined for Lyot stop generation mode.')
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--Rotation matrix used on center coordinates.
rotMat = [cosd(clock_deg), -sind(clock_deg); sind(clock_deg), cosd(clock_deg)];

dx = Dbeam/Nbeam; % pixel width [pupil diameters]
% if(strcmpi(centering,'pixel'))
%     Narray = ceil_even(Dbeam/dx+1); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
% else
%     Narray = ceil_even(Dbeam/dx); %--number of points across output array. Same size as width when interpixel centered.
% end

%--number of points across output array.
if(strcmpi(centering,'pixel'))
    Narray = ceil_even(Nbeam*max([magFac, 1.0]) + 1 + 2*magFac*Nbeam*max(abs([xShear, yShear]))); % Sometimes requires two more pixels when pixel centered.
elseif(strcmpi(centering,'interpixel'))
    Narray = ceil_even(Nbeam*max([magFac, 1.0]) + 2*magFac*Nbeam*max(abs([xShear, yShear]))); %--No zero-padding needed if beam is centered between pixels
end
if(isfield(changes,'Narray')); Narray = changes.Narray; end

Darray = Narray*dx; %--width of the output array [pupil diameters]
bdf = 1; % beam diameter factor in output array
wl_dummy = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -Darray/2/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -dx/2; 
    case {'pixel'}
        cshift = 0;
        if(flagRot180)
            cshift = -dx;
        end
end

%--DATA FROM THE VISIO FILE
D0 = 8; % inches, pupil diameter in Visio file
x0 = -26; % inches, pupil center in x in Visio file
y0 = 20.25; % inches, pupil center in y in Visio file
Dconv = Dbeam/D0; % conversion factor from inches (Visio units) to pupil diameters 

%--INITIALIZE PROPER
bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);
prop_set_antialiasing(33);

if(flagLyot == false)
    %--PRIMARY MIRROR (OUTER DIAMETER)
    ra_OD = (Dbeam/2 - pad_idod)*magFac;
    cx_OD = 0 + cshift + xShear;
    cy_OD = 0 + cshift + yShear;
    bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);

    %--SECONDARY MIRROR (INNER DIAMETER)
    ra_ID = ((2.46/2)*Dconv + pad_idod)*magFac;
    cx_ID = 0 + cshift + xShear;
    cy_ID = 0 + cshift + yShear;
    bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);%, cx, cy, norm)
else
    %--PRIMARY MIRROR (OUTER DIAMETER)
    ra_OD = OD/2*magFac;
    cx_OD = 0 + cshift + xShear;
    cy_OD = 0 + cshift + yShear;
    bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);

    %--SECONDARY MIRROR (INNER DIAMETER)
    ra_ID = ID/2*magFac;
    cx_ID = 0 + cshift + xShear;
    cy_ID = 0 + cshift + yShear;
    bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);%, cx, cy, norm)
end

if(isfield(changes,'wStrut'))
    wStrut = changes.wStrut; %--strut width [pupil diameters]
else
    wStrut = (0.209/D0)*Dbeam;
end

%--STRUTS
if wStrut > 0
    Nstrut = 6;
    angStrutVec = [77.56, -17.56, -42.44, 42.44, 17.56, 102.44];
    wStrutVec = (wStrut + 2*pad_strut)*ones(1, Nstrut);
    lStrut = 3.6*Dconv;
    xcStrutVec = Dconv*[(-24.8566 - x0), (-23.7187 - x0), (-24.8566 - x0), (-27.1434 - x0), (-28.2813 - x0), (-27.1434 - x0)];
    ycStrutVec = Dconv*[(22.2242 - y0), (20.2742 - y0), (18.2758 - y0), (18.2758 - y0), (20.2742 - y0), (22.2242 - y0)];

    for iStrut = 1:Nstrut
        angDeg = angStrutVec(iStrut) + clock_deg; % degrees
        wStrut = magFac*(wStrutVec(iStrut) + 2*pad_strut);
        lStrutIn = magFac*lStrut;
        xc = magFac*(xcStrutVec(iStrut));
        yc = magFac*(ycStrutVec(iStrut));
        cxy = rotMat*[xc; yc];
        xc = cxy(1)+xShear;
        yc = cxy(2)+yShear;
        bm = prop_rectangular_obscuration(bm, lStrutIn, wStrut, 'XC',xc+cshift, 'YC',yc+cshift, 'ROTATION',angDeg);%, norm)
    end
end

mask = ifftshift(abs(bm.wf));
if(flagRot180)
    mask = rot90(mask,2);
end

end %--END OF FUNCTION
