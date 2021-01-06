
function mask = falco_gen_rounded_bowtie_LS(Nbeam, ID, OD, Rfillet, upsampleFactor, angDeg, clockDeg, centering, varargin)

    xShear = 0.;
    yShear = 0.;
    icav = 0;  % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
            case {'xc'}
                icav = icav + 1;
                xShear = varargin{icav};  % x offset [pupil diameters]
            case {'yc'}
                icav = icav + 1;
                yShear(1) = varargin{icav};  % y offset [pupil diameters]
            otherwise
                error('falco_gen_rounded_bowtie_FPM: Unknown keyword: %s\n', varargin{icav});
        end
    end
    
    DbeamUM = 10e3; % [microns]
    stepSize = 10; % [microns]
    xcMicron = xShear*DbeamUM; % [microns]
    ycMicron = yShear*DbeamUM; % [microns]

    % Masks have 180-degree rotational symmetry
    clockDeg = clockDeg + 90; %--keep original convention of openings being along vertical axis. 
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