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
% inputs.num_strut
% inputs.strut_angs
% inputs.strut_width

function PUPIL = falco_gen_pupil_Simple( input )
%falco_gen_pupil_SCDA Generates a simple pupil.
%   Function may be used to generate circular, annular, and simple on-axis 
%   telescopes with radial struts. 

    hg_expon = 1000; % hyper-gaussian exponent for anti-aliasing 
    hg_expon_spider = 100; % hyper-gaussian exponent for anti-aliasing 

    if(isfield(input,'centering'))
        centering = input.centering; disp('Checkpoint 1');
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
            [X,Y] = meshgrid(-(N-1)/2:(N-1)/2); disp('hi')
        otherwise
            [X,Y] = meshgrid(-N/2:N/2-1);
    end
    [THETA,RHO] = cart2pol(X,Y); 
    
    % Make sure the inputs make sense
    if(ID > OD)
        error('Pupil generation error: Inner diameter larger than outer diameter.');
    end
    
    % Create inner and outer circles
    if(ID > 0)
        PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon) - exp(-(RHO/(apRad*ID)).^hg_expon);
    else
        PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon);
    end
    
    % Create spiders 
    if(input.strut_width > 0)
        
        numSpiders = input.num_strut;
        angs = input.strut_angs;
        
        if(numSpiders~=numel(angs))
            error('Pupil generation error: ''strut_angs'' should be an array of length ''num_strut''');
        end
        
        halfwidth = input.strut_width*2*apRad;
        for ang = angs
           PUPIL = PUPIL.*(1-exp(-(RHO.*sin(THETA-ang*pi/180)/halfwidth).^hg_expon_spider).*...
               (RHO.*cos(THETA-ang*pi/180)>0));
        end
    end
    
    
end

