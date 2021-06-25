% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Inputs: 
% inputs.Nbeam
% inputs.OD
% inputs.ID
% inputs.Nstrut
% inputs.angStrut
% inputs.wStrut
%
% Corrected on 2018-08-16 by A.J. Riggs to compute 'beam_diam_fraction' correctly.

function PUPIL = falco_gen_pupil_Simple_PROPER( input )
%falco_gen_pupil_SCDA Generates a simple pupil.
%   Function may be used to generate circular, annular, and simple on-axis 
%   telescopes with radial struts. 

    hg_expon = 1000; % hyper-gaussian exponent for anti-aliasing 
    hg_expon_spider = 100; % hyper-gaussian exponent for anti-aliasing 

    if(isfield(input,'centering'))
        centering = input.centering; 
    else
        centering = 'pixel';
    end
    
    N = input.Npad;%Number of samples in NxN grid 
    OD = input.OD; % pupil outer diameter, can be < 1
    ID = input.ID; % central obscuration radius 
    apRad = input.Nbeam/2; % aperture radius in samples 
    
    %Create coordinates
    switch centering
        case 'interpixel'
            [X,Y] = meshgrid(-(N-1)/2:(N-1)/2);
        otherwise
            [X,Y] = meshgrid(-N/2:N/2-1);
    end
    [THETA,RHO] = cart2pol(X,Y); 
    
    % Make sure the inputs make sense
    if(ID > OD)
        error('Pupil generation error: Inner diameter larger than outer diameter.');
    end
    
    % Create inner and outer circles in PROPER
    
    %--INITIALIZE PROPER
    Dbeam = 1; %--diameter of beam (normalized to itself)
    dx = Dbeam/input.Nbeam;
    Narray = N;
    Darray = Narray*dx;
    wl_dummy = 1e-6; %--dummy value
    bdf = Dbeam/Darray; %--beam diameter fraction
    xshift = 0; %--x-shear of pupil
    yshift = 0; %--y-shear of pupil
    bm = prop_begin(Dbeam, wl_dummy, Narray,'beam_diam_fraction',bdf);
    
    switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
        case {'interpixel'}
            cshift = -dx/2; 
        case {'pixel'}
            cshift = 0;
    end
    
    %--PRIMARY MIRROR (OUTER DIAMETER)
    ra_OD = OD/2;
    cx_OD = 0 + cshift + xshift;
    cy_OD = 0 + cshift + yshift;
    bm = prop_circular_aperture(bm, ra_OD,'cx',cx_OD,'cy',cy_OD);%, cx, cy, norm);

    if(ID > 0)
        %--SECONDARY MIRROR (INNER DIAMETER)
        ra_ID = ID/2;
        cx_ID = 0 + cshift + xshift;
        cy_ID = 0 + cshift + yshift;
        bm = prop_circular_obscuration(bm, ra_ID,'cx',cx_ID,'cy',cy_ID);%, cx, cy, norm)
    end
    PUPIL = ifftshift(abs(bm.wf));
    
    % Create spiders 
    if(input.wStrut > 0)
        
        numSpiders = input.Nstrut;
        angs = input.angStrut;
        
        if(numSpiders~=numel(angs))
            error('Pupil generation error: ''angStrut'' should be an array of length ''Nstrut''');
        end
        
        halfwidth = input.wStrut*2*apRad;
        for ang = angs
           PUPIL = PUPIL.*(1-exp(-(RHO.*sin(THETA-ang*pi/180)/halfwidth).^hg_expon_spider).*...
               (RHO.*cos(THETA-ang*pi/180)>0));
        end
    end
    
end