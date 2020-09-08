
function [ arrayOut ] = gmtPupil_addCirc( cenrow, cencol, segDiam, segID, arrayIn)
%gmtPupil_addCirc Adds circular segment to arrayIn, centered at (cenrow,
%cencol), with diametet 'segDiam'.
%   Inputs:
%   cenrow - row of hexagon center (samples)
%   cencol - column of hexagon center (samples)
%   segDiam - diameter of the segment (samples)
%   segID - segment inner diameter (central obscuration) fraction 
%   arrayIn - Input array
%   
%   Coordinate system origin: (rows/2+1, cols/2+1)

    [rows,cols]=size(arrayIn);

    [X,Y] = meshgrid(-cols/2:cols/2-1,-rows/2:rows/2-1); % Grids with Cartesian (x,y) coordinates 

    RHOprime = sqrt((X-cencol).^2+(Y-cenrow).^2);

    CIRC = exp(-(RHOprime/(segDiam/2)).^1000) - exp(-(RHOprime/(segID*segDiam/2)).^1000);
    
    arrayOut = arrayIn + CIRC;
    
end

