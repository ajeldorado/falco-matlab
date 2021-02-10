%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [zca, maxr, maxt] = prop_noll_zernikes(maxz, varargin)
%        [zca, maxr, maxt] = prop_noll_zernikes(maxz, varargin)
% Return a cell array in which each element contains the Zernike
% polynomial equation corresponding to the index of that element.
% The polynomials are orthonormal for an unobscured circular aperture.
% They follow the ordering convention of Noll (J Opt Soc Am, 66, 207 (1976)).
% The equations contain the variables "r" (normalized radius) and
% "t" (azimuth angle in radians).
% The polynomials have an RMS of 1.0 relative to a mean of 0.0.
%
% Outputs:
% zca  = A cell array with each element containing a Zernike polynomial.
%        For example:
%        zca = prop_noll_zernikes(6)
%        will display:
%        '1'
%        '2 * (r) .* cos(t)'
%        '2 * (r) .* sin(t)'
%        'sqrt(3) * (2*r.^2 - 1)'
%        'sqrt(6) * (r.^2) .* sin(2*t)'
%        'sqrt(6) * (r.^2) .* cos(2*t)'
% maxr = Maximum radial exponent
% maxt = Maximum multiplier of the azimuth angle.
%
% Required inputs:
% maxz = Maximum number of Zernike polynomials to return.
%
% Optional inputs:
% 'compact'           : if set, the equations are returned using the
%                       naming convention for terms assumed by prop_zernikes
%
% Note: to convert a cell array element to a character string, use char:
% zstr = char(zca(iz));

% Revision history:
% 2005 Feb     jek
% 2016 Mar 30  gmg  Translated to Matlab
% 2017 Feb 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cmpt = 0;
  if size(varargin, 2) > 0
    switch lower(varargin{1})
      case {'cmpt', 'compact'}
        cmpt = 1;
      otherwise
        error('prop_noll_zernikes: Unknown keyword: %s\n', varargin{1});
    end
  end

  if cmpt == 1
    pcs  = '_pow_';                     % power character string
  else
    pcs  = '.^';                        % power character string
  end

  zca(1 : maxz, 1) = cellstr('');       % Create output cell array
  maxr = 0;                             % maximum radial exponent
  maxt = 0;                             % maximum azimuth multiplier
  iz   = 1;                             % index of Zernike polynomial
  izn  = 0;                             % Zernike polynomial parameter n

  while iz <= maxz
    for izm = mod(izn, 2) : 2 : izn     % Zernike polynomial parameter m
      for izp = 0 : fix(izm ~= 0)
        if izn ~= 0

          if izm ~= 0
            val  = 2 * (izn + 1);       % normalization value squared
          else
            val = izn + 1;              % normalization value squared
          end

          srv  = fix(sqrt(val));        % integer normalization value
          if val == srv^2
            rfs  = sprintf('%d * (', srv);             % radial function string
          else
            rfs  = sprintf('sqrt(%d) * (', val);       % radial function string
          end

        else
          zca(iz, 1) = cellstr('1');
          iz   = iz + 1;
          continue
        end                     % if izn ~= 0

        for izs = 0 : ((izn - izm) / 2)
          top  = (-1)^izs * factorial(izn - izs);
          bot  = factorial(izs) ...
               * factorial((izn + izm) / 2 - izs) ...
               * factorial((izn - izm) / 2 - izs);
          rc   = top / bot;                 % radial coefficient
          rcs  = sprintf('%d.0', abs(rc));  % radial coefficient string
          rp   = izn - 2 * izs;             % radial power
          rps  = sprintf('%d', rp);         % radial power string

          maxr = max(maxr, rp);

          if top ~= 0
            if izs == 0
              if rc < 0
                sgns = '-';
              else
                sgns = '';
              end
            else
              if rc < 0
                sgns = ' - ';
              else
                sgns = ' + ';
              end
            end                 % if izs == 0

            if rp == 0
              rfs  = [rfs sgns rcs];
            elseif rp == 1
              if rc ~= 1
                rfs  = [rfs sgns rcs '*r'];
              else
                rfs  = [rfs sgns 'r'];
              end
            else
              if rc ~= 1
                rfs  = [rfs sgns rcs '*r' pcs rps];
              else
                rfs  = [rfs sgns 'r' pcs rps];
              end
            end                 % if rp == 0
          end                   % if top ~= 0
        end                     % for izs = 0 : ((izn - izm) / 2)

        maxt = max(maxt, izm);

        if izm ~= 1
          aas  = sprintf('%d*t', izm);          % azimuth angle string
        else
          aas  = 't';                           % azimuth angle string
        end
        if mod(iz, 2) == 0
          if cmpt == 1
            afs  = sprintf(' .* cos%dt', izm);  % azimuth function string
          else
            afs  = [' .* cos(' aas ')'];        % azimuth function string
          end
        else
          if cmpt == 1
            afs  = sprintf(' .* sin%dt', izm);  % azimuth function string
          else
            afs  = [' .* sin(' aas ')'];        % azimuth function string
          end
        end                     % if mod(iz, 2) == 0

        if izm == 0
          zca(iz, 1) = cellstr([rfs ')']);
        else
          zca(iz, 1) = cellstr([rfs ')' afs]);
        end

        iz   = iz + 1;
        if iz > maxz
          break
        end
      end                       % for izp = 0 : izm
      if iz > maxz
        break
      end
    end                         % for izm = mod(izn, 2) : 2 : izn
    izn  = izn + 1;
  end                           % while iz <= maxz
end                     % function prop_noll_zernikes
