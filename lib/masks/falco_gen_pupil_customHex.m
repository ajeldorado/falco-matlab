% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Inputs: 
% inputs.Nbeam - Beam diameter
% inputs.OD - Flat-to-flat diameter (fraction of Nbeam)
% inputs.ID - Inner diameter (fraction of Nbeam)
% inputs.Nstrut - Number of struts
% inputs.angStrut - Array of angles for the radial struts (deg)
% inputs.wStrut - strut width (fraction of Nbeam)


function PUPIL = falco_gen_pupil_customHex( input )
%gen_pupil_SCDA Generates a simple pupil.
%   Function may be used to generate circular, annular, and simple on-axis 
%   telescopes with radial struts. 

    hg_expon = 1000; % hyper-gaussian exponent for anti-aliasing 
    hg_expon_spider = 100; % hyper-gaussian exponent for anti-aliasing 

    N = input.Npad;%Number of samples in NxN grid 
    OD = input.OD; % pupil outer diameter, can be < 1
    ID = input.ID; % central obscuration radius 
    apRad = input.Nbeam/2; % aperture radius in samples 
    
    %Create coordinates
    [X,Y] = meshgrid(-N/2:N/2-1);
    [THETA,RHO] = cart2pol(X,Y); 
   
    input.apDia = input.Nbeam;
    if(isfield(input,'pistons'))
        PUPIL0 = hexSegMirror_getField( input );
    else
        PUPIL0 = hexSegMirror_getSupport( input );
    end
    
    % Create inner and outer circles
    if(ID > 0)
        PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon) - exp(-(RHO/(apRad*ID)).^hg_expon);
    else
        PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon);
    end
    
    PUPIL = PUPIL.*PUPIL0;
    
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

