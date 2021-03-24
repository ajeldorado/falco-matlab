%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function [img, pixx] = prop_run(flnm, wlu, nx, varargin)
%        [img, pixx] = prop_run(flnm, wlu, nx, varargin)
% Execute one or more instances of a Proper prescription.
% Accepts a wavelength or a file of wavelengths and weights.
%
% Outputs:
% img  = array of output images returned by the Proper prescription.
% pixx = sampling of img (m / pixel)
%        It is the responsibility of the prescription to return this
%        value (which is returned from prop_end).
%
% Required inputs:
% flnm = filename (excluding extension) of Matlab Proper prescription
% wlu  = Either the wavelength (um) or the name of a text file
%        containing a list of wavelength and weight pairs.
%        In the latter case, the prescription is run for each wavelength
%        and the results added together with the respective weights.
% nx   = size of computational grid.  Most efficient if a power of 2.
%
% Optional inputs:
% 'passvalue'         = structure containing input parameters
%                       (passed to Matlab Proper prescription).
% 'phase_offset'      : if set, then a phase offset is added as the
%                       wavefront is propagated.  For instance, if a
%                       wavefront is propagated over a distance of 1/4
%                       wavelength, a phase offset of pi/2 radians will
%                       be added.  This is useful in cases where the
%                       offset between separate beams that may be
%                       combined later may be important (e.g., the
%                       individual arms of an interferometer).
%                       By default, a phase offset is not applied.
% 'print_intensity'   : if set, then print total intensity
% 'quiet'             : if set, then intermediate messages and surface
%                       labels will not be printed
% 'table'             : if set, then prints out a table of sampling and
%                       beam size for each surface
% 'verbose'           : if set, then informational messages will be sent

% 2005 Feb     jek  created idl routine
% 2015 Jun 22  gmg  Matlab translation
% 2017 Apr 10  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon

  if nargin < 1
    error('Proper:PROP_RUN', 'Need filename of Proper prescription');
  end
  if nargin < 2
    error('Proper:PROP_RUN', 'Need wavelength(s) in micrometers');
  end
  if nargin < 3
    error('Proper:PROP_RUN', 'Need size of computational grid');
  end

% Set default values of input parameters
  do_table = 0;                 % do not print table
  iprm = 0;                     % no passvalues
  print_it = 1;                 % print intermediate messages, surface labels
  print_total_intensity = 0;    % do not print total intensity
  prop_phase_offset = 0;        % do not apply phase offset
  prop_verbose = 0;             % do not print informational messages

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'prm', 'passvalue'}
        icav = icav + 1;
        iprm = 1;
        prm  = varargin{icav};  % structure containing input parameters
      case {'qt', 'quiet'}
        print_it = 0;           % do not print intermediate messages
      case {'poff', 'phase_offset'}
        prop_phase_offset = 1;  % apply phase offset
      case {'pint', 'print_intensity'}
        print_total_intensity = 1;      % print total intensity
      case {'tbl', 'table'}
        do_table = 1;           % print table
      case {'vrbs', 'verbose'}
        prop_verbose = 1;       % print informational messages
      otherwise
        error('prop_run: Unknown keyword: %s\n', varargin{icav});
    end
  end

  funh = str2func(flnm);        % function handle for Proper prescription
  tic;
  if ischar(wlu)                % then wlu is the name of a text file
    dlmt = ' ';                 % delimiter used in file
    nhdr = 1;                   % number of header lines in file
    wla  = importdata(wlu, dlmt, nhdr); % read wavelengths and weights
    wla.data(:, 1) = wla.data(:, 1) * 1.0e-6;
    nwl  = size(wla.data, 1);   % number of wavelengths in file
    img = zeros(nx, nx);
    for iwl = 1 : nwl
      if print_it               % print wavelength and throughput
        fprintf(1, 'Lambda = %16.7e  Throughput = %13.5f\n', ...
          wla.data(iwl, 1), wla.data(iwl, 2));
      end
      if iprm                   % prm is given
        [imgt, pixx] = funh(wla.data(iwl), nx, prm);
      else
        [imgt, pixx] = funh(wla.data(iwl), nx);
      end
      img  = img + imgt * wla.data(iwl, 2);
    end
  else                          % wlu is a single wavelength
    wlm  = wlu * 1.0e-6;        % convert wlu to meters
    if print_it                 % print wavelength and throughput
      fprintf(1, 'Lambda = %16.7e  Throughput = %13.5f\n', wlm, 1.0d0);
    end
    if iprm                     % prm is given
      [img, pixx] = funh(wlm, nx, prm);
    else
      [img, pixx] = funh(wlm, nx);
    end
  end
  if do_table
    prop_table;
  end
  tdet = toc;

  if print_it                   % print elapsed time
    fprintf(1, 'Total elapsed time (seconds) = %16.9f\n', tdet);
  end
end                             % function prop_run
