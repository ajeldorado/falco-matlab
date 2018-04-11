% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ numOfSegments ] = hexSegMirror_numSegments( numRings )
%hexSegMirror_numSegments Returns the number of segments in a hexagonal
%mirror with a given number of rings, 'numRings'
%   numRings - number of rings in the mirror

% loop through rings and add up the number of segments
numOfSegments = 0;
for ringNum = 0:numRings
    numOfSegments = numOfSegments + 1; 
    for face = 1:6

        stepnum = 1;

        while(stepnum<=ringNum)
            if(face==6 && stepnum==ringNum)
                %disp(['Finished ring ',num2str(ringNum)]);
            else
                numOfSegments = numOfSegments + 1;
            end
            stepnum = stepnum + 1;
            
        end
    end
end


end

