% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ OUT ] = hexSegMirror_getSurfHeight( hexMirror_struct, lambda0 )
%hexSegMirror_getSurfHeight Returns the surface height of a 
% hexagonally segmented mirror assuming a central wavelength lambda0
%   Input: hexMirror_struct - Structure with the following variables 
%   apDia - flat to flat aperture diameter (samples)
%   wGap - width of the gap between segments (samples)
%   numRings - number of rings in the segmented mirror (samples)
%   N - size of NxN computational grid 
%   pistons - Segment pistons in waves
%   tiltxs - Tilts on segment in horizontal direction (waves/apDia)
%   tiltys - Tilts on segment in vertical direction (waves/apDia)
%   lambda0 - Central wavelength 
%
%   Returns:
%   OUT - the surface height of the hex. segmented mirror (same units as
%   lambda0)

apDia = hexMirror_struct.apDia; % flat to flat aperture diameter (samples)
wGap = hexMirror_struct.wGap; % samples
numRings = hexMirror_struct.numRings;% Number of rings in hexagonally segmented mirror 
N = hexMirror_struct.Npad;
pistons = hexMirror_struct.pistons;
tiltxs = hexMirror_struct.tiltxs; 
tiltys = hexMirror_struct.tiltys; 

OUT = zeros(N);

hexFlatDiam = (apDia-numRings*2*wGap)/(2*numRings+1);
hexSep = hexFlatDiam + wGap;

segNum = 1;
for ringNum = 0:numRings

    cenrow = ringNum*hexSep;
    cencol = 0;
    
    [ OUT ] = hexSegMirror_addHexSegment( cenrow, cencol, numRings, apDia, ...
                wGap, pistons(segNum), tiltxs(segNum), tiltys(segNum), OUT);
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
                disp(['Finished ring ',num2str(ringNum)]);
            else
                [ OUT ] = hexSegMirror_addHexSegment( cenrow, cencol, numRings, apDia, ...
                            wGap, pistons(segNum), tiltxs(segNum), tiltys(segNum), OUT);
                segNum = segNum + 1;
            end
            stepnum = stepnum + 1;
            
        end
    end
end

OUT = angle(OUT)/(4*pi)*lambda0;

end

