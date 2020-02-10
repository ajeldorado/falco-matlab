% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ arrayOut ] = hexSegMirror_addHexSegment( cenrow, cencol, numRings, ...
    apDia, wGap, piston, tiltx, tilty, loworder_struct, arrayIn)
%hexSegMirror_addHexSegment Adds hexagonal mirror segment to arrayIn, 
% centered at (cenrow, cencol). The full mirror have numRings rings of 
% hexagonal segments, flat-to-flat diameter (in samples) of apDia, 
% wGap (in samples) between segments, and piston, tiltx, and tilty
% phase offsets. Piston is in units of waves. tiltx and tilty are waves
% across the full flat-to-flat pupil diameter apDia. 
%
%   Inputs:
%   cenrow - row of hexagon center (samples)
%   cencol - column of hexagon center (samples)
%   numRings - number of rings in the segmented mirror (samples)
%   apDia - flat to flat aperture diameter (samples)
%   wGap - width of the gap between segments (samples)
%   piston - Segment piston in waves
%   tiltx - Tilt on segment in horizontal direction (waves/apDia)
%   tilty - Tilt on segment in vertical direction (waves/apDia)
%   loworder_struct - Structure that defines segment-level low order aberrations. 
%                      loworder_struct.noll_indices - list of noll indices
%                      loworder_struct.waves_rms - Zernike coefficients in
%                                                   units of waves rms. 
%   arrayIn - Input array
%   
%   Coordinate system origin: (rows/2+1, cols/2+1)

    hexFlatDiam = (apDia-numRings*2*wGap)/(2*numRings+1);
    % hexRad = hexFlatDiam/sqrt(3);% center to vertex
    hexSep = hexFlatDiam + wGap;

    [rows,cols]=size(arrayIn);

    [X,Y] = meshgrid(-cols/2:cols/2-1,-rows/2:rows/2-1); % Grids with Cartesian (x,y) coordinates 

    RHOprime = sqrt((X-cencol).^2+(Y-cenrow).^2);
    THETA = atan2(Y-cenrow,X-cencol); 

    HEXamp = exp(-(RHOprime.*sin(THETA)/(hexFlatDiam/2)).^1000)...
            .*exp(-(RHOprime.*cos(THETA-pi/6)/(hexFlatDiam/2)).^1000)...
            .*exp(-(RHOprime.*cos(THETA+pi/6)/(hexFlatDiam/2)).^1000);
    
    HEXphz = and(and(RHOprime.*sin(THETA)<=(hexSep/2),...
             RHOprime.*cos(THETA-pi/6)<=(hexSep/2)),...
             RHOprime.*cos(THETA+pi/6)<=(hexSep/2))...
            .*exp(1i*2*pi*piston)...
            .*exp(1i*2*pi*tiltx/apDia*(X-cencol))...
            .*exp(1i*2*pi*tilty/apDia*(Y-cenrow));
    
    % Add Zernike polynomials according to loworder_struct
    if(isfield(loworder_struct,'noll_indices') && numel(loworder_struct.noll_indices) > 0)
        
        Zcount = 1;
        for noll_index = loworder_struct.noll_indices 
            
            hexCircumRad = hexFlatDiam/sqrt(3); % Circumscribed radius of the hex.
            Zcoeff = loworder_struct.waves_rms(Zcount); % Zernike coefficient in waves RMS. 
            mask = logical(HEXamp); % Logical mask over which the RMS is calculated (a hex segment).
            
            Z = propcustom_gen_zernike( noll_index, hexCircumRad, RHOprime, THETA);
            
            Znorm = Z./sqrt(mean((Z(mask)).^2));% normalize by RMS across hex 
            Zwaves = Zcoeff*Znorm;% Convert to waves RMS
            HEXphz = HEXphz.*exp(1i*2*pi*Zwaves); % Apply wavefront 
            
            Zcount = Zcount + 1; 
        end
    end
        
    arrayOut = arrayIn + HEXamp.*HEXphz;
    
end

