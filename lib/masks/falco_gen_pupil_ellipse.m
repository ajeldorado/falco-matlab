% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Inputs: 
% inputs.Nbeam - Number of samples across the beam 
% inputs.OD - Outer diameter (fraction of Nbeam)
% inputs.ID - Inner diameter (fraction of Nbeam)
% inputs.Nstrut - Number of struts
% inputs.angStrut - Array of struct angles (deg)
% inputs.wStrut - Strut widths (fraction of Nbeam)


function pupil = falco_gen_pupil_ellipse(inputs)
%falco_gen_pupil_SCDA Generates a simple pupil.
%   Function may be used to generate elliptical, annular, and simple on-axis 
%   telescopes with radial struts. 

%     hg_expon = 1000; % hyper-gaussian exponent for anti-aliasing 
    hg_expon_spider = 100; % hyper-gaussian exponent for anti-aliasing 

    if(isfield(inputs,'centering'))
        centering = inputs.centering;
    else
        centering = 'pixel';
    end
        
    % Create coordinates
    N = inputs.Narray;%Number of samples in NxN grid 
    switch centering
        case 'interpixel'
            [X,Y] = meshgrid(-(N-1)/2:(N-1)/2);
        otherwise
            [X,Y] = meshgrid(-N/2:N/2-1);
    end
    [THETA,RHO] = cart2pol(X,Y); 
    
    % Create ellipse for primary mirror
    inputsPrimary.Nbeam = inputs.Nbeam; % max aperture radius in samples
    inputsPrimary.Narray = inputs.Narray; % Number of samples across in square output array
    inputsPrimary.radiusX = inputs.primaryRadiusX; % x-radius of ellipse [pupil diameters]
    inputsPrimary.radiusY = inputs.primaryRadiusY; % y-radius of ellipse [pupil diameters]
    inputsPrimary.clockingDegrees = inputs.primaryClockingDegrees; % clocking of the primary [degrees]
    primary = falco_gen_ellipse(inputsPrimary);
    
    if(inputs.secondaryRadiusX > 0 && inputs.secondaryRadiusY > 0)
        inputsSecondary.Nbeam = inputs.Nbeam; % max aperture radius in samples
        inputsSecondary.Narray = inputs.Narray; % Number of samples across in square output array
        inputsSecondary.radiusX = inputs.secondaryRadiusX; % x-radius of ellipse [pupil diameters]
        inputsSecondary.radiusY = inputs.secondaryRadiusY; % y-radius of ellipse [pupil diameters]
        inputsSecondary.clockingDegrees = inputs.secondaryClockingDegrees; % clocking of the primary [degrees]
        secondary = 1 - falco_gen_ellipse(inputsSecondary);        
    else
        secondary = 1;
    end

    pupil = primary.*secondary;
    
%     if(ID > 0)
%         PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon) - exp(-(RHO/(apRad*ID)).^hg_expon);
%     else
%         PUPIL = exp(-(RHO/(apRad*OD)).^hg_expon);
%     end

    % Create spiders 
    if(inputs.wStrut > 0)
        anglesRadians = inputs.angStrut*(pi/180.);

        halfwidth = inputs.wStrut*inputs.Nbeam/2.0;
        for ang = anglesRadians
           pupil = pupil.*(1-exp(-(RHO.*sin(THETA-ang)/halfwidth).^hg_expon_spider).*...
               (RHO.*cos(THETA-ang)>0));
        end
    end
    
end