% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Inputs structure:
% inputs.Nbeam - Number of samples across the beam 
% inputs.OD - Outer diameter (fraction of Nbeam)
% inputs.ID - Inner diameter (fraction of Nbeam)
% inputs.Nstrut - Number of struts
% inputs.angStrut - Array of struct angles (deg)
% inputs.wStrut - Strut widths (fraction of Nbeam)
% inputs.stretch - Create an elliptical aperture by changing Nbeam along
%                   the horizontal direction by a factor of stretch (PROPER
%                   version isn't implemented as of March 2019).

function pupil = falco_gen_pupil_Simple(inputs)
%   Generate a simple pupil.
%   Can be used to generate circular, annular, and simple on-axis 
%   telescopes with radial struts. 

    Nbeam = inputs.Nbeam;
    Narray = inputs.Npad; %Number of samples in NxN grid 
    OD = inputs.OD; % pupil outer diameter, can be < 1
    if(isfield(inputs, 'ID')); ID = inputs.ID; else; ID = 0; end % central obscuration radius
    if(isfield(inputs, 'wStrut')); wStrut = inputs.wStrut; else; wStrut = 0; end % strut width [pupil diameters]
    if(isfield(inputs, 'angStrut')); angStrutVec = inputs.angStrut; else; angStrutVec = []; wStrut = 0; end % vector of strut orientation angles
    if(isfield(inputs,'centering')); centering = inputs.centering; else; centering = 'pixel'; end
    if(isfield(inputs,'stretch')); xStretch = inputs.stretch; else; xStretch = 1; end
    if(isfield(inputs, 'clocking')); clocking = inputs.clocking; else; clocking = 0; end % clocking [degrees]
    if(isfield(inputs, 'xShear')); xShear = inputs.xShear; else; xShear = 0; end % x-shear [pupil diameters]
    if(isfield(inputs, 'yShear')); yShear = inputs.yShear; else; yShear = 0; end % y-shear [pupil diameters]
    if(isfield(inputs, 'flagHG')); flagHG = inputs.flagHG; else; flagHG = false; end % whether to use hyper-Gaussians to generate the pupil instead. (Old behavior)

    
    if ~flagHG % By default, don't use hyger-gaussians for anti-aliasing the edges.
    
        % Create outer aperture
        inpOuter.Nbeam = Nbeam;
        inpOuter.Narray = Narray;
        inpOuter.radiusX = xStretch*0.5*OD;
        inpOuter.radiusY = 0.5*OD;
        inpOuter.centering = centering;
        inpOuter.clockingDegrees = clocking;
        inpOuter.xShear = xShear;
        inpOuter.yShear = yShear;
        apOuter = falco_gen_ellipse(inpOuter);

        % Create inner aperture
        if ID > 0
            inpInner.Nbeam = Nbeam;
            inpInner.Narray = Narray;
            inpInner.radiusX = xStretch*0.5*ID;
            inpInner.radiusY = 0.5*ID;
            inpInner.centering = centering;
            inpInner.clockingDegrees = clocking;
            inpInner.xShear = xShear;
            inpInner.yShear = yShear;
            apInner = 1 - falco_gen_ellipse(inpInner);
        elseif ID > OD
            error('Pupil generation error: Inner diameter cannot be larger than outer diameter.');
        else
            apInner = 1;
        end
        
        % Create struts in PROPER
        if isempty(angStrutVec)
            
            apStruts = 1;
            
        else
            %--INITIALIZE PROPER
            Dbeam = 1; %--diameter of beam (normalized to itself)
            dx = Dbeam/inputs.Nbeam;
            Darray = Narray*dx;
            wl_dummy = 1e-6; %--dummy value
            bdf = Dbeam/Darray; %--beam diameter fraction
            bm = prop_begin(Dbeam, wl_dummy, Narray, 'beam_diam_fraction', bdf);

            switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
                case {'interpixel'}
                    cshift = -dx/2; % [pupil diameters]  
                case {'pixel'}
                    cshift = 0;
            end

            %--STRUTS
            lStrut = 0.6; % [pupil diameters]
            rcStrut0 = lStrut / 2.0;
            for iStrut = 1:length(angStrutVec)
                ang =  angStrutVec(iStrut) + clocking;
                bm = prop_rectangular_obscuration(bm, lStrut, wStrut,...
                    'XC', rcStrut0*cosd(ang)+cshift+xShear, 'YC', rcStrut0*sind(ang)+cshift+yShear, 'ROT', ang);
            end
            apStruts = ifftshift(abs(bm.wf));
            
        end
        
        % Combine all features
        pupil = apOuter.*apInner.*apStruts;      
        
        
    else % Use hyper-Gaussians
        
        hg_expon = 1000; % hyper-gaussian exponent for anti-aliasing 
        hg_expon_spider = 100; % hyper-gaussian exponent for anti-aliasing 
        apRad = inputs.Nbeam/2; % aperture radius in samples 
    
        %Create coordinates
        switch centering
            case 'interpixel'
                [X,Y] = meshgrid(-(Narray-1)/2:(Narray-1)/2);
            otherwise
                [X,Y] = meshgrid(-Narray/2:Narray/2-1);
        end
        [THETA,RHO] = cart2pol(X/xStretch,Y); 
        
        % Create inner and outer circles
        if(ID > 0 && ID < OD)
            pupil = exp(-(RHO/(apRad*OD)).^hg_expon) - exp(-(RHO/(apRad*ID)).^hg_expon);
        else
            pupil = exp(-(RHO/(apRad*OD)).^hg_expon);
        end

        % Create spiders 
        if inputs.wStrut > 0

            halfwidth = inputs.wStrut*apRad;
            for ang = inputs.angStrut
               pupil = pupil.*(1-exp(-(RHO.*sin(THETA-ang*pi/180)/halfwidth).^hg_expon_spider).*...
                   (RHO.*cos(THETA-ang*pi/180)>0));
            end

        end
        

    end
    
end