% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the LUVOIR B architecture pupil mask in Matlab using PROPER.
%
% REQUIRED INPUTS: 
% bmi  = beam structure input (used only to get sampling info)
% nr   = number of rings of hexagons in aperture
%        (e.g., 1 = a central hexagon surrounded by a ring of hexagons)
% hr   = distance from the center of a hexagonal segment to a vertex (m)
% hs   = distance between centers of adjacent hexagonal segments (m)
%
% OPTIONAL INPUTS (same as in PROPER):
% 'xc'                = aperture center to wavefront center x (m)
% 'yc'                = aperture center to wavefront center y (m)
%        By default, the aperture is centered within the wavefront.
% 'rotation'          = counter-clockwise rotation of the aperture about
%                       its center (degrees)
% 'darkcenter'        : if set, the central hexagonal segment
%                        transmission will be set to 0.0
%
% OUTPUTS:
%  ap:     2-D square array of the amplitude for the specified hex aperture 

function ap = falco_hex_aperture_LUVOIR_B(bmi, nr, hr, hs, varargin)
%        ap = falco_hex_aperture_LUVOIR_B(bmi, nr, hr, hs, varargin)
% Return an image containing a hexagonal mask consisting of multiple
% hexagons. The hexagons have antialiased edges. This routine does
% not modify the wavefront.
%
% Outputs:
% ap   = aperture transmission map
% % Required inputs:
% bmi  = beam structure input (used only to get sampling info)
% nr   = number of rings of hexagons in aperture
%        (e.g., 1 = a central hexagon surrounded by a ring of hexagons)
% hr   = distance from the center of a hexagonal segment to a vertex (m)
% hs   = distance between centers of adjacent hexagonal segments (m)
%
% Optional inputs:
% 'xc'                = aperture center to wavefront center x (m)
% 'yc'                = aperture center to wavefront center y (m)
%        By default, the aperture is centered within the wavefront.
% 'rotation'          = counter-clockwise rotation of the aperture about
%                       its center (degrees)
% 'darkcenter'        : if set, the central hexagonal segment
%                        transmission will be set to 0.0
%
% 2007 Jul     jek  created idl routine
% 2016 Aug 31  gmg  Matlab translation
% 2017 Mar 16  gmg  Revised for keyword/value for optional inputs
% 2018 Mar 27  ar   Changed to be specific to the LUVOIR A 5 aperture.
% 2020 Oct 08  ar   Changed to be specific to the LUVOIR B aperture.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
cx   = 0.0;           % pixel coordinates of aperture center X
cy   = 0.0;           % pixel coordinates of aperture center Y
useDC   = false;             % no dark center
rot  = 0.0;           % counter-clockwise rotation of aperture (deg)

icav = 0;             % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % pixel coordinates of aperture center X
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % pixel coordinates of aperture center Y
      case {'rot', 'rotation'}
        icav = icav + 1;
        rot  = varargin{icav};  % counter-clockwise rotation of aperture (deg)
      case {'dc', 'darkcenter'}
        useDC = true;               % dark center
      otherwise
        error('falco_hex_aperture_LUVOIR_B: Unknown keyword: %s\n', varargin{icav});
    end
end

ap   = zeros(size(bmi.wf));  % aperture array

counter = 0;
  
for iring = 0 : nr

    hx   = hs *  iring * cosd(30d0);       % hexagonal segment center x (m)
    hy   = hs * (iring * cosd(60d0) - nr); % hexagonal segment center y (m)

    for iseg = 0 : 2 * nr - iring

        % create hexagonal segment on one side
        if (iring ~= 0) || ~((iseg == nr) && useDC)
            hexx = cx + hx * cosd(rot) - hy * sind(rot);
            hexy = cy + hx * sind(rot) + hy * cosd(rot);
            counter = counter + 1;
            if ~any(counter == [1, 9, 52, 60])
                ap = ap + prop_polygon(bmi,6,hr,'cx',hexx,'cy',hexy,'rot',rot);
            end

        end

        % create hexagonal segment on opposite side
        if iring ~= 0
            hexx = cx - hx * cosd(rot) - hy * sind(rot);
            hexy = cy - hx * sind(rot) + hy * cosd(rot);
            counter = counter + 1;
            if ~any(counter==[1, 9, 53, 61])
                ap = ap + prop_polygon(bmi,6,hr,'cx',hexx,'cy',hexy,'rot',rot);
            end
        end

        hy   = hy + hs;

    end
end
end % EOF
