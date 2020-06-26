
function [ OUT ] = gmtPupil_getSupport( gmtPupil_struct )
%hexSegMirror_getSupport Returns the support of the pupil function defined
%by a hexagonally segmented mirror 
%   Input: hexMirror_struct - Structure with the following variables 
%   apDia - flat to flat aperture diameter (samples)
%   Npad - size of NxN computational grid 
%   missingSegments - list of zeros and ones indicating if each segment is present 
%   offset - Array with lateral offsets in samples (NumSegments x 2)

apDia = gmtPupil_struct.apDia; % flat to flat aperture diameter (samples)
Npad = gmtPupil_struct.Npad;
segDia = gmtPupil_struct.segDia;
segID = gmtPupil_struct.segID; % central obscuration fraction of outer segments 

if(isfield(gmtPupil_struct,'offset'))
    offset = gmtPupil_struct.offset;
end

if(isfield(gmtPupil_struct,'missingSegments'))
    missingSegments = gmtPupil_struct.missingSegments;
else
    missingSegments = ones(1,7);
end

N1 = 2^nextpow2(apDia);
OUT = zeros(N1);

seg_num = 1; % segment number

%%-- Add central segment 
cenrow = 0;
cencol = 0;
if(isfield(gmtPupil_struct,'offset'))
    cenrow = cenrow + offset(seg_num,1);
    cencol = cencol + offset(seg_num,2);
end
if(missingSegments(seg_num)==1)
    [ OUT ] = gmtPupil_addCirc( cenrow,cencol, segDia, segID, OUT );
end

%%-- Add outer segments 
dist_to_center = apDia/2 - segDia/2; % Distance to center of outer segments
for seg_num = 2:7
    
    angle_to_center = (seg_num-1)*pi/3 + pi/6;
    cenrow = dist_to_center.*sin(-1*angle_to_center);
    cencol = dist_to_center.*cos(-1*angle_to_center);
   
    if(isfield(gmtPupil_struct,'offset'))
        cenrow = cenrow + offset(seg_num,1);
        cencol = cencol + offset(seg_num,2);
    end

    if(missingSegments(seg_num)==1)
        [ OUT ] = gmtPupil_addCirc( cenrow,cencol, segDia, segID, OUT );
    end

end

OUT = padOrCropEven(OUT,Npad);

end

