%
% Function to generate GMT pupil
%
% Inputs: 
% inputs.Nbeam - Number of points across the pupil diameter
% input.Npad - Number of points in the computational grid (Npad x Npad)
% input.centralID - Central obscuration fraction (wrt Nbeam)
% input.segDia - Segment diameter (units of samples)
% input.segID - Central obscuration fraction of individual segments
% input.OD - Pupil outer diameter, can be < 1
% input.offset - Array with lateral offsets in samples (NumSegments x 2)
% input.pistons - Vector of pistons per segment (waves)
% input.tiltxs - Vector of x-tilt per segment (waves/apDia)
% input.tiltys - Vector of y-tilt per segment (waves/apDia)
% input.loworder_struct - Structure to define segment-level low order aberrations. 
%                      loworder_struct(i).noll_indices - list of noll indices for the ith segment
%                      loworder_struct(i).waves_rms - Zernike coefficients for the ith segment (waves rms)

function PUPIL = falco_gen_pupil_gmtPupil( input )

    hg_expon = 1000; % hyper-gaussian exponent for anti-aliasing 

    N = input.Npad;%Number of samples in NxN grid 
    OD = input.OD; % pupil outer diameter, can be < 1
    centralID = input.centralID; % central obscuration radius 
    apRad = input.Nbeam/2; % aperture radius in samples 
    
    %Create coordinates
    [X,Y] = meshgrid(-N/2:N/2-1);
    [~,RHO] = cart2pol(X,Y); 
   
    % Create the pupil (full field of just the support) 
    input.apDia = input.Nbeam;
    if(isfield(input,'pistons') || isfield(input,'tiltxs') || ...
                isfield(input,'tiltys') || isfield(input,'loworder_struct'))
        PUPIL0 = gmtPupil_getField( input );
    else
        PUPIL0 = gmtPupil_getSupport( input );
    end
    
    % Create inner and outer circles
    if(centralID > 0)
        PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon) - exp(-(RHO/(apRad*centralID)).^hg_expon);
    else
        PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon);
    end
    
    PUPIL = PUPIL.*PUPIL0;
    
end