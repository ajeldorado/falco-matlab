

function mask = falco_gen_rounded_bowtie_FPM(rhoInner, rhoOuter, rocFillet, pixresFPM, angDeg, clockDeg, upsampleFactor, centering, varargin)

    Nbeam = 2*rhoOuter*pixresFPM;
    ID = 2*rhoInner/(2*rhoOuter);
    OD = 1;%2*rhoOuter;
    Rfillet = rocFillet/(2*rhoOuter); % [aperture diameters]

    xShear = 0.;
    yShear = 0.;

    icav = 0;             % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
            case {'xc'}
                icav = icav + 1;
                xShear = varargin{icav}; % x offset (lambda/D)
            case {'yc'}
                icav = icav + 1;
                yShear(1) = varargin{icav}; % y offset (lambda/D)
            otherwise
                error('falco_gen_rounded_bowtie_FPM: Unknown keyword: %s\n', varargin{icav});
        end
    end

    DbeamUM = 1e3; % [microns]
    stepSize = 1; % [microns]
%     xcMicron = xShear*DbeamUM; % [microns]
%     ycMicron = yShear*DbeamUM; % [microns]
    xcMicron = xShear/(2*rhoOuter)*DbeamUM; % [microns]
    ycMicron = yShear/(2*rhoOuter)*DbeamUM; % [microns]


% maxOffset = max(abs([xShear, yShear]));
% % Nbeam = rhoOuter*2*pixresFPM;
% switch centering
%     case 'pixel'
%         Narray = ceil_even(2*pixresFPM*maxOffset + Nbeam + 1);
%     case 'interpixel'
%         Narray = ceil_even(2*pixresFPM*maxOffset + Nbeam);
%     otherwise
%         error('centering must be pixel or interpixel')
% end



%     xShear = 0.;
%     yShear = 0.;
%     icav = 0;  % index in cell array varargin
%     while icav < size(varargin, 2)
%         icav = icav + 1;
%         switch lower(varargin{icav})
%             case {'xc'}
%                 icav = icav + 1;
%                 xShear = varargin{icav};  % x offset [pupil diameters]
%             case {'yc'}
%                 icav = icav + 1;
%                 yShear(1) = varargin{icav};  % y offset [pupil diameters]
%             otherwise
%                 error('falco_gen_rounded_bowtie_FPM: Unknown keyword: %s\n', varargin{icav});
%         end
%     end
    


    % Masks have 180-degree rotational symmetry
    clockDeg = mod(clockDeg, 180);
    
    % Depending on the clocking, call a different function because of
    % the tangent function
    if clockDeg < (90 - angDeg/2.0)
        [mask, ~] = falco_gen_rounded_bowtie_horizontal_vector(Nbeam, ID, OD, Rfillet, angDeg, clockDeg, upsampleFactor, DbeamUM, stepSize, xcMicron, ycMicron, centering);
    elseif clockDeg > (90 - angDeg/2.0) && clockDeg < (90 + angDeg/2.0)
        [mask, ~] = falco_gen_rounded_bowtie_vertical_vector(Nbeam, ID, OD, Rfillet, angDeg, clockDeg-90, upsampleFactor, DbeamUM, stepSize, xcMicron, ycMicron, centering);
    elseif clockDeg > (90 + angDeg/2.0)
        [mask, ~] = falco_gen_rounded_bowtie_horizontal_vector(Nbeam, ID, OD, Rfillet, angDeg, clockDeg-180, upsampleFactor, DbeamUM, stepSize, xcMicron, ycMicron, centering);
    else
        error('Sorry, not set up yet to handle vertical edges of the bowtie.')
    end

end % END OF FUNCTION

