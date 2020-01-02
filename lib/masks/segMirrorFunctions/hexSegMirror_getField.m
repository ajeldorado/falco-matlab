% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ OUT ] = hexSegMirror_getField( hexMirror_struct )
%hexSegMirror_getField Returns the complex reflectance of the pupil function defined
%by a hexagonally segmented mirror 
%   Input: hexMirror_struct - Structure with the following variables 
%   apDia - flat to flat aperture diameter (samples)
%   wGap - width of the gap between segments (samples)
%   numRings - number of rings in the segmented mirror (samples)
%   N - size of NxN computational grid 
%   pistons - Segment pistons in waves
%   tiltxs - Tilts on segment in horizontal direction (waves/apDia)
%   tiltys - Tilts on segment in vertical direction (waves/apDia)
%   loworder_struct - Structure that defines segment-level low order aberrations. 
%                      loworder_struct(i).noll_indices - list of noll indices for the ith segment. 
%                      loworder_struct(i).waves_rms - Zernike coefficients for the ith segment in units of waves rms. 

apDia = hexMirror_struct.apDia; % flat to flat aperture diameter (samples)
wGap = hexMirror_struct.wGap; % samples
numRings = hexMirror_struct.numRings;% Number of rings in hexagonally segmented mirror 
N = hexMirror_struct.Npad;
pistons = hexMirror_struct.pistons;
tiltxs = hexMirror_struct.tiltxs; 
tiltys = hexMirror_struct.tiltys; 
if(isfield(hexMirror_struct,'loworder_struct'))
    loworder_struct = hexMirror_struct.loworder_struct;
else
    loworder_struct = nan(1,hexSegMirror_numSegments( numRings ));
end

if(isfield(hexMirror_struct,'missingSegments'))
    missingSegments = hexMirror_struct.missingSegments;
else
    missingSegments = ones(1,hexSegMirror_numSegments( numRings ));
end

N1 = 2^nextpow2(apDia);
OUT = zeros(N1);

hexFlatDiam = (apDia-numRings*2*wGap)/(2*numRings+1);
hexSep = hexFlatDiam + wGap;

segNum = 1;
for ringNum = 0:numRings

    cenrow = ringNum*hexSep;
    cencol = 0;
    
    if(missingSegments(segNum)==1)
        [ OUT ] = hexSegMirror_addHexSegment( cenrow, cencol, numRings, apDia, ...
                    wGap, pistons(segNum), tiltxs(segNum), tiltys(segNum), loworder_struct(segNum), OUT);
    end
    segNum = segNum + 1;
    
    for face = 1:6
        
        step_dir = pi/6*(2*face+5);
        steprow = hexSep*sin(step_dir);
        stepcol = hexSep*cos(step_dir);
        
        stepnum = 1;

        while(stepnum<=ringNum)
            cenrow = cenrow + steprow;
            cencol = cencol + stepcol;
            if(face==6 && stepnum==ringNum)
                %disp(['Finished ring ',num2str(ringNum)]);
            else
                if(missingSegments(segNum)==1)
                    [ OUT ] = hexSegMirror_addHexSegment( cenrow, cencol, numRings, apDia, ...
                                wGap, pistons(segNum), tiltxs(segNum), tiltys(segNum), loworder_struct(segNum), OUT);
                end
                segNum = segNum + 1;
            end
            stepnum = stepnum + 1;
            
        end
    end
end
OUT = padOrCropEven(OUT,N);
end

