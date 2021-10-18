
function [ OUT ] = gmtPupil_getField( gmtPupil_struct )
%gmtPupil_getField Returns the complex reflectance of the pupil function 
%   Input: gmtPupil_struct - Structure with the following variables: 
%   apDia - full aperture diameter (samples)
%   Npad - size of the computational grid (Npad x Npad)
%   segDia - individual segment diameter (samples)
%   segID - central obscuration fraction of outer segments 
%   pistons - Segment pistons in waves
%   tiltxs - Tilts on segment in horizontal direction (waves/apDia)
%   tiltys - Tilts on segment in vertical direction (waves/apDia)
%   offset - Array with lateral offsets in samples (NumSegments x 2)
%   loworder_struct - Structure that defines segment-level low order aberrations. 
%                      loworder_struct(i).noll_indices - list of noll indices for the ith segment. 
%                      loworder_struct(i).waves_rms - Zernike coefficients for the ith segment (waves rms) 

apDia = gmtPupil_struct.apDia; % flat to flat aperture diameter (samples)
Npad = gmtPupil_struct.Npad;
segDia = gmtPupil_struct.segDia;
segID = gmtPupil_struct.segID; % central obscuration fraction of outer segments 

pistons = gmtPupil_struct.pistons;
tiltxs = gmtPupil_struct.tiltxs; 
tiltys = gmtPupil_struct.tiltys; 

if(isfield(gmtPupil_struct,'offset'))
    offset = gmtPupil_struct.offset;
end

if(isfield(gmtPupil_struct,'loworder_struct'))
    loworder_struct = gmtPupil_struct.loworder_struct;
else
    loworder_struct = nan(1,7);
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
	[ OUT ] = gmtPupil_addCircSegment( cenrow, cencol, apDia, segDia, segID, ...
                    pistons(seg_num), tiltxs(seg_num), tiltys(seg_num), loworder_struct(seg_num), OUT);
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
        [ OUT ] = gmtPupil_addCircSegment( cenrow, cencol, apDia, segDia, segID, ...
                    pistons(seg_num), tiltxs(seg_num), tiltys(seg_num), loworder_struct(seg_num), OUT);
    end

end

OUT = padOrCropEven(OUT,Npad);


end

