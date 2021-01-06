% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate a bowtie FPM in Matlab using PROPER.
%  -Outside the outer ring is opaque.
%   -If rhoOuter = infinity, then the outer ring and azimuthal parts are omitted 
%     and the mask is cropped down to the size of the inner spot.
%  -The output array is the smallest size that fully contains the mask.
%
% Modified by A.J. Riggs on 2018-04-26 to create a bowtie instead of just an annulus. 
% Modified by A.J. Riggs on 2017-10-25 from a file to generate a pupil 
%   to a file that generates an FPM. 
% Modified by A.J. Riggs on 2017-02-09 to be a function with specifyable
%   amounts of magnification, clocking, and translation. 
%   Now keeping the pupil edges non-binary.
% Created by A.J. Riggs on 2017-01-11.
%
%
%--INPUTS:
% pixresFPM: resolution in pixels per lambda_c/D
% rhoInner:   radius of inner FPM amplitude spot (in lambda_c/D)
% rhoOuter:   radius of outer opaque FPM ring (in lambda_c/D). Set to
%             infinity for an occulting-spot only FPM
% centering:  centering of the mask on the array. 'pixel' or 'interpixel'
% ang:        opening angle on each side (left and right) of the FPM (degrees)
%
%--OUTPUT:
% mask:    cropped-down, 2-D FPM representation. amplitude only 

function mask = falco_gen_bowtie_FPM(inputs)

% OPTIONAL INPUTS
centering = 'pixel';  %--Default to pixel centering
xOffset = 0;
yOffset = 0;
Rfillet = 0; %--fillet radius of curvature [lambda/D]
upsampleFactor = 101; %--upsample factor for grayscale edges of mask with fillets.
clocking = 0; % [degrees]

if(isfield(inputs,'centering')); centering = inputs.centering; end
if(isfield(inputs, 'xOffset')); xOffset = inputs.xOffset; end  % [lambda0/D]
if(isfield(inputs, 'yOffset')); yOffset = inputs.yOffset; end  % [lambda0/D]
if(isfield(inputs,'Rfillet')); Rfillet = inputs.Rfillet; end
if(isfield(inputs,'upsampleFactor')); upsampleFactor = inputs.upsampleFactor; end
if(isfield(inputs,'clocking')); clocking = inputs.clocking; end

% %--DEBUGGING ONLY: HARD-CODED INPUTS
% inputs.pixresFPM = 6; %--pixels per lambda_c/D
% inputs.rhoInner = 2.8; % radius of inner FPM amplitude spot (in lambda_c/D)
% inputs.rhoOuter = 10.1; % radius of outer opaque FPM ring (in lambda_c/D)
% inputs.centering = 'pixel';
% inputs.ang = 65; % (degrees)

FPMampFac = 0; %--Amplitude transmission of the center spot
pixresFPM = inputs.pixresFPM; %--pixels per lambda_c/D
rhoInner = inputs.rhoInner; % radius of inner FPM amplitude spot (in lambda_c/D)
rhoOuter = inputs.rhoOuter; % radius of outer opaque FPM ring (in lambda_c/D)
ang = inputs.ang; %--Opening angle on each side of the bowtie (degrees)


if Rfillet > 0
    
    mask = falco_gen_rounded_bowtie_FPM(rhoInner, rhoOuter, Rfillet, pixresFPM, ang, clocking, upsampleFactor, centering, 'xc', xOffset, 'yc', yOffset);
    
    if(isfield(inputs, 'Narray'))
        mask = pad_crop(mask, inputs.Narray);
    end

else

    dx = 1/pixresFPM; %--lambda_c/D per pixel. "UL" for unitless
    maxAbsOffset = max([abs(xOffset), abs(yOffset)]);
    if(isinf(rhoOuter))
        if(strcmpi(centering, 'interpixel'))
            Narray = ceil_even(2*rhoInner/dx + 2*pixresFPM*maxAbsOffset); % number of points across the inner diameter of the FPM.
        elseif(strcmpi(centering, 'pixel'))
            Narray = ceil_even(2*rhoInner/dx + 2*pixresFPM*maxAbsOffset + 1); % number of points across the inner diameter of the FPM. Another half pixel added for pixel-centered masks.
        end

        Dmask = 2*pixresFPM*rhoInner; %--Diameter of the mask

    else
        if(strcmpi(centering, 'interpixel'))
            Narray = ceil_even(2*rhoOuter/dx + 2*pixresFPM*maxAbsOffset); % number of points across the outer diameter of the FPM. 
        elseif(strcmpi(centering, 'pixel'))
            Narray = ceil_even(2*rhoOuter/dx + 2*pixresFPM*maxAbsOffset + 1); % number of points across the outer diameter of the FPM. Another half pixel added for pixel-centered masks.
        end

        Dmask = 2*pixresFPM*rhoOuter; %--Diameter of the mask

    end

    if(isfield(inputs, 'Narray'))
        Narray = inputs.Narray;
    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    %--GENERAL SPECS
    Darray = Narray*dx; % width of array in lambda_c/D
    wl_dummy = 1e-6; % wavelength (m); Dummy value--no propagation here, so not used.

    switch centering % 0 for pixel-centered FPM, or -diam/Narray for inter-pixel centering
        case {'interpixel'}
            cshift = -dx/2; 
        case {'pixel'}
            cshift = 0;
    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    %--INITIALIZE PROPER
    bdf = Dmask/Darray; %--beam diameter factor in output array
    bm = prop_begin(Dmask, wl_dummy, Narray, 'beam_diam_fraction', bdf);
    prop_set_antialiasing(ceil_odd(upsampleFactor));

    if(isinf(rhoOuter) == false)
        %--Outer opaque ring of FPM
        ra_OD = (rhoOuter);
        cx_OD = 0 + cshift + xOffset;
        cy_OD = 0 + cshift + yOffset;
        bm = prop_circular_aperture(bm, ra_OD, 'cx', cx_OD, 'cy', cy_OD);%, cx, cy, norm);

         %--Create the bowtie region
        if(ang < 180)
            %--Top part
            Lside = 2* ra_OD;
            xvert = cshift + xOffset + [0, Lside*cosd(ang/2), Lside*cosd(ang/2), -Lside*cosd(ang/2), -Lside*cosd(ang/2)];
            yvert = cshift + yOffset + [0, Lside*sind(ang/2), Lside,            Lside,            Lside*sind(ang/2)];
            bowtieTop = prop_irregular_polygon(bm, xvert, yvert, 'DARK');

            %--Bottom part
            xvert = cshift + xOffset + [0, Lside*cosd(ang/2), Lside*cosd(ang/2), -Lside*cosd(ang/2), -Lside*cosd(ang/2)];
            yvert = cshift + yOffset + -1*[0, Lside*sind(ang/2), Lside,            Lside,            Lside*sind(ang/2)];
            bowtieBottom = prop_irregular_polygon(bm, xvert, yvert, 'DARK');
        else
            bowtieTop = 1;
            bowtieBottom = 1;
        end

    else
        bowtieTop = 1;
        bowtieBottom = 1;

    end

    %--Inner spot of FPM (Amplitude transmission can be nonzero)
    ra_ID = rhoInner; 
    cx_ID = 0 + cshift + xOffset;
    cy_ID = 0 + cshift + yOffset;
    innerSpot = prop_ellipse(bm, ra_ID, ra_ID, 'cx', cx_ID, 'cy', cy_ID, 'DARK')*(1-FPMampFac) + FPMampFac;

    mask = ifftshift(abs(bm.wf)); %--undo PROPER's fftshift
    mask = mask.*innerSpot; %--Include the inner FPM spot
    mask = mask.*bowtieTop.*bowtieBottom; %--Include the azimuthal part
end

end %--END OF FUNCTION
