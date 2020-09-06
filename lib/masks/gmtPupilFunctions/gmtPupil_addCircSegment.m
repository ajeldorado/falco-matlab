
function [ arrayOut ] = gmtPupil_addCircSegment( cenrow, cencol, apDia, ... 
    segDiam, segID, piston, tiltx, tilty, loworder_struct, arrayIn)
%gmtPupil_addCircSegment Adds segment to arrayIn, 
% centered at (cenrow, cencol). The full mirror has diameter (in samples) 
% of apDia, and piston, tiltx, and tilty phase offsets. 
% Piston is in units of waves. tiltx and tilty are waves
% across the full flat-to-flat pupil diameter apDia. 
%
%   Inputs:
%   cenrow - row of hexagon center (samples)
%   cencol - column of hexagon center (samples)
%   segDiam - diameter of the segment (samples)
%   segID - segment inner diameter (central obscuration) fraction 
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


    [rows,cols]=size(arrayIn);

    [X,Y] = meshgrid(-cols/2:cols/2-1,-rows/2:rows/2-1); % Grids with Cartesian (x,y) coordinates 

    RHOprime = sqrt((X-cencol).^2+(Y-cenrow).^2);
    THETA = atan2(Y-cenrow,X-cencol); 
    
    segRad = segDiam/2; % Radius of a segment
    
    CIRCamp = exp(-(RHOprime/(segRad)).^1000) - exp(-(RHOprime/(segID*segRad)).^1000);
    
    CIRCphz = and(RHOprime <= segRad,RHOprime >= segID*segRad)...
            .*exp(1i*2*pi*piston)...
            .*exp(1i*2*pi*tiltx/apDia*(X-cencol))...
            .*exp(1i*2*pi*tilty/apDia*(Y-cenrow));
    
    % Add Zernike polynomials according to loworder_struct
    if(isfield(loworder_struct,'noll_indices') && numel(loworder_struct.noll_indices) > 0)
        
        Zcount = 1;
        for noll_index = loworder_struct.noll_indices 
            
            
            Zcoeff = loworder_struct.waves_rms(Zcount); % Zernike coefficient in waves RMS. 
            mask = logical(CIRCamp); % Logical mask over which the RMS is calculated (a hex segment).
            
            Z = propcustom_gen_zernike( noll_index, segRad, RHOprime, THETA);
            
            Znorm = Z./sqrt(mean((Z(mask)).^2));% normalize by RMS across hex 
            Zwaves = Zcoeff*Znorm;% Convert to waves RMS
            CIRCphz = CIRCphz.*exp(1i*2*pi*Zwaves); % Apply wavefront 
            
            Zcount = Zcount + 1; 
        end
    end
        
    arrayOut = arrayIn + CIRCamp.*CIRCphz;
    
end

