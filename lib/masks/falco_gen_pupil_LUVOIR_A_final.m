% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate the LUVOIR Design A (Final) telescope pupil from 
% 2018 in Matlab using PROPER.
%
% Coordinates and dimensions of the primary, secondary, and hex segments
%   are from Matthew Bolcar (NASA GSFC).
% Coordinates and dimensions of the secondary mirror support struts were a
%   best-fit match by A.J. Riggs by matching PROPER-made rectangles to the 
%   pupil file from Matthew Bolcar (NASA GSFC).
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126

function pupil = falco_gen_pupil_LUVOIR_A_final(inputs)

if(isfield(inputs,'centering')); centering = inputs.centering; else; centering = 'pixel'; end %--Centering of pupil
if(isfield(inputs,'magFac')); mag = inputs.magFac; else; mag = 1; end % Pupil Magnification
if(isfield(inputs,'wGap_m')); hexGap0 = inputs.wGap_m; else;  hexGap0 = 6e-3; end % Gap between primary mirror segments [meters]. Default is 6mm
if(isfield(inputs,'clock_deg')); clockDeg = inputs.clock_deg; else; clockDeg = 0; end
if(isfield(inputs,'flagLyot')); flagLyot = inputs.flagLyot; else; flagLyot = false; end
if(isfield(inputs,'flagRot180')); flagRot180 = inputs.flagRot180; else; flagRot180 = false; end
if(isfield(inputs,'xShear')); xShear = inputs.xShear; else; xShear = 0; end % [pupil diameters]
if(isfield(inputs,'yShear')); yShear = inputs.yShear; else; yShear = 0; end % [pupil diameters]
shearMax = max(abs([xShear, yShear])); % [pupil diameters]

%%--(Optional) Lyot stop mode (concentric, circular ID and OD)
if(flagLyot == true)
    if(isfield(inputs,'ID'))
        ID = inputs.ID;
    else
        error('inputs.ID must be defined for Lyot stop generation mode.')
    end

    if(isfield(inputs,'OD'))
        OD = inputs.OD;
    else
        error('inputs.OD must be defined for Lyot stop generation mode.')
    end
end


%--Rotation matrix used on center coordinates.
rotMat = [cosd(clockDeg), -sind(clockDeg); sind(clockDeg), cosd(clockDeg)];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%--USER INPUTS
Nbeam   = inputs.Nbeam;     % number of points across the incoming beam  
Nap   = Nbeam; % number of points across FULL usable pupil
width_hex0 = 1.2225; %-- flat-to-flat (m)
wFlatToFlat = (12*width_hex0 + 11*hexGap0); % flat-to-flat width
Dap = 1.01824*wFlatToFlat; % minimum diameter of circumscribing circle.
% fprintf('Diameter = %.5f\n meters', Dap);
dx = Dap/Nap;
xShearM = Dap*xShear; % meters
yShearM = Dap*yShear; % meters
% dxDrawing = 1.2225/158; % (m) %--In actual drawing, 158 pixels across the 1.2225m, so dx_pixel =
%dx_drawing. strut in drawing is 19 pixels across, so dx_strut =
%19*dx_drawing = 0.1470m
if(isfield(inputs,'wStrut')); wStrutM = inputs.wStrut * Dap; else; wStrutM = 0.15*mag; end % load as pupil diameters and convert to meters

if(strcmpi(centering,'pixel'))
    Narray = ceil_even((1+2*shearMax)*Nbeam*max([mag, 1])+1); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = ceil_even((1+2*shearMax)*Nbeam*max([mag, 1])); %--number of points across output array. Same size as width when interpixel centered.
end
% Darray = Narray*dx;

%--For PROPER 
wl_dummy = 1e-6; % wavelength (m)
bdf = Nbeam/Narray; %--beam diameter factor in output array
      
hexWidth = mag*width_hex0; %-- flat-to-flat (m)
nrings = 6;
hexRadius = 2/sqrt(3)*hexWidth/2;
hexGap = mag*hexGap0; % (m)
hexSep = hexWidth + hexGap; % distance from center to center of neighboring segments
switch centering % 0 for pixel-centered pupil, or -diam/np for inter-pixel centering
    case {'interpixel','even'}
        cShift = -dx/2; 
    case {'pixel','odd'}
        cShift = 0;
        if(flagRot180)
            cShift = -dx;
        end
end

% Use PROPER to generate the hexagonal mirrors and rectangular struts.
bm = prop_begin(Dap, wl_dummy, Narray,'beam_diam_fraction',bdf);
prop_set_antialiasing(33)
if ~flagLyot % (Entrance pupil mode)
    
    ap = falco_hex_aperture_LUVOIR_A_5(bm,nrings,hexRadius,hexSep,'XC',cShift+xShearM,'YC',cShift+yShearM,'DARKCENTER', 'ROTATION', clockDeg);

else % (Lyot stop mode)
    ap = 1;
    
    %--OUTER DIAMETER
    ra_OD = mag*(OD/2.)*Dap;
    bm = prop_circular_aperture(bm, ra_OD, 'XC', cShift + xShearM, 'YC', cShift + yShearM);

    %--INNER DIAMETER
    ra_ID = mag*(ID/2.)*Dap;
    bm = prop_circular_obscuration(bm, ra_ID, 'XC', cShift + xShearM, 'YC', cShift + yShearM);
end

% Struts
if wStrutM > 0
    xc = 0;
    yc = mag*Dap/4;
    xyc = rotMat*[xc; yc]; xc = xyc(1); yc = xyc(2);
    bm = prop_rectangular_obscuration(bm, wStrutM, 7*hexWidth, 'XC', xc+xShearM+cShift, 'YC', yc+yShearM+cShift, 'ROTATION', clockDeg);
    len_1b = (sqrt(93)+0.5)*hexRadius;

    xc = 1.5*hexRadius;
    yc = -11*sqrt(3)/4*hexRadius;
    xyc = rotMat*[xc; yc]; xc = xyc(1); yc = xyc(2);
    bm = prop_rectangular_obscuration(bm, wStrutM, len_1b, 'XC', xc+xShearM+cShift, 'YC', yc+yShearM+cShift, 'ROT', 12.7+clockDeg);

    xc = -1.5*hexRadius;
    yc = -11*sqrt(3)/4*hexRadius;
    xyc = rotMat*[xc; yc]; xc = xyc(1); yc = xyc(2);
    bm = prop_rectangular_obscuration(bm, wStrutM, len_1b, 'XC', xc+xShearM+cShift, 'YC', yc+yShearM+cShift, 'ROT', -12.7+clockDeg);
end

pupil = ifftshift(abs(bm.wf)).*ap;
pupil(pupil > 1) = 1; %--Get rid of overlapping segment edges at low resolution if the gap size is zero.
if(flagRot180)
    pupil = rot90(pupil, 2);
end

end %--END OF FUNCTION
