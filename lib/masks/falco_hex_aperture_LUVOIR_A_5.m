% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the LUVOIR A architecture 5/6 pupil mask in Matlab using PROPER.
%
% Created in March 2018 by A.J. Riggs.
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

function ap = falco_hex_aperture_LUVOIR_A_5(bmi, nr, hr, hs, varargin)
%        ap = prop_hex_aperture(bmi, nr, hr, hs, varargin)
% Return an image containing a hexagonal mask consisting of multiple
% hexagons.  This is useful for modeling systems with multisegmented
% mirrors, such as the Keck or JWST telescopes.  The hexagons have
% antialiased edges.  This routine does not modify the wavefront.
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

% 2007 Jul     jek  created idl routine
% 2016 Aug 31  gmg  Matlab translation
% 2017 Mar 16  gmg  Revised for keyword/value for optional inputs
% 2018 Mar 27  ar   Changed to be specific to the LUVOIR A 5 aperture.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0;           % pixel coordinates of aperture center X
  cy   = 0.0;           % pixel coordinates of aperture center Y
  dc   = 0;             % no dark center
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
        dc   = 1;               % dark center
      otherwise
        error('prop_hex_aperture: Unknown keyword: %s\n', varargin{icav});
    end
  end

  [ny, nx] = size(bmi.wf);      % number of array points
  ap   = zeros(ny, nx);         % aperture array

  counter = 0;
  Xs = 2*nr*hs*(-(nx-1)/2:(nx-1)/2)/nx;
  Ys = 2*nr*hs*(-(ny-1)/2:(ny-1)/2)/ny;
  
  for ir = 0 : nr
      
    hx   = hs *  ir * cosd(30d0);       % hexagonal segment center x (m)
    hy   = hs * (ir * cosd(60d0) - nr); % hexagonal segment center y (m)

    for iseg = 0 : 2 * nr - ir
      if(0) 
        fprintf('Now skipping segment %d \n',counter);
      else
    % create hexagonal segment on one side
          if (ir ~= 0) || ~((iseg == nr) && dc)
            hexx = cx + hx * cosd(rot) - hy * sind(rot);
            hexy = cy + hx * sind(rot) + hy * cosd(rot);
            counter = counter + 1;
            if(any(counter==[1 12 113 114 125 126]))
            else
                ap   = ap + prop_polygon(bmi,6,hr,'cx',hexx,'cy',hexy,'rot',rot);
            end

          end
          
    % create hexagonal segment on opposite side
          if (ir ~= 0)
            hexx = cx - hx * cosd(rot) - hy * sind(rot);
            hexy = cy - hx * sind(rot) + hy * cosd(rot);
            counter = counter + 1;
            if(any(counter==[1 12 113 114 125 126]))
            else
                ap   = ap + prop_polygon(bmi,6,hr,'cx',hexx,'cy',hexy,'rot',rot);
            end
          end

          hy   = hy + hs;
      end
        
    end                 % for iseg = 0 : 2 * nr - ir
  end                   % for ir = 0 : nr
end                     % prop_hex_aperture
