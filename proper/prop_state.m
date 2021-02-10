%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_state(bm)
%        bm = prop_state(bm)
% Save the current state for the current wavelength, if one doesn't
% already exist.  If one does, read it in and use it to define the
% current wavefront.  The current contents of the wavefront array
% structure and some other info necessary to represent the current state
% of propagation are saved to a file.
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)

% 2005 Feb     jek  created idl routine
% 2016 Apr 11  gmg  Matlab translation
% 2017 Apr 13  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  propcommon
  if (numel(save_state) > 0) & (save_state ~= 0)

% If the state for the current wavelength doesn't exist, write it out.

    if ~prop_is_statesaved(bm)
      prop_savestate(bm);

% If it does exist, read it in.
    else
      [ny, nx] = size(bm.wf);
      nv   = 2 * ny * nx;                       % number of values in sarr array
      sarr = zeros(nv, 1);
      fno  = [num2str(bm.wl * 1.0e9) statefile] ;
      fid  = fopen(fno, 'r');                   % file ID
      if fid == -1
        error('Proper:PROP_STATE', ...
          'Cannot open file %s.\n', fno);
      end
      prop_total_intensity  = fread(fid, 1, '*double');          % wavefront power
      sarr = fread(fid, nv, '*double');         % complex wavefront array
      iz   = 0;
      for iy = 1 : ny
        for ix = 1 : nx
          iz   = iz + 2;
          bm.wf(iy, ix) = complex(sarr(iz - 1), sarr(iz));
        end
      end
      bm.wl = fread(fid, 1, '*double');         % wavelength (m)
      bm.dx = fread(fid, 1, '*double');         % grid sampling (m / pixel)
      tpo  = fread(fid, 7, '*char')';           % TypeOld
      if tpo == 'OUTSIDE'
        bm.TypeOld = 'OUTSIDE';
      else
        bm.TypeOld = 'INSIDE_';
      end
      rfs  = fread(fid, 6, '*char')';           % RefSurf
      if rfs == 'SPHERI'
        bm.RefSurf = 'SPHERI';
      else
        bm.RefSurf = 'PLANAR';
      end
      bm.Rbeam = fread(fid, 1, '*double');      % beam radius of curvature (m)
      bm.RbeamInf = fread(fid, 1, 'int16=>double');
      bm.pz = fread(fid, 1, '*double');         % beam position z (m)
      bm.w0_pz = fread(fid, 1, '*double');      % beam waist position z (m)
      bm.w0 = fread(fid, 1, '*double');         % beam waist radius (m)
      bm.zRay = fread(fid, 1, '*double');       % Rayleigh distance (m)
      ptyp = fread(fid, 18, '*char')';          % propagation type
      if     ptyp == 'INSIDE__TO_OUTSIDE'
        bm.PropType = 'INSIDE__to_OUTSIDE';
      elseif ptyp == 'OUTSIDE_TO_INSIDE_'
        bm.PropType = 'OUTSIDE_to_INSIDE_';
      else
        bm.PropType = 'INSIDE__to_INSIDE_';
      end
      bm.fr = fread(fid, 1, '*double');         % focal ratio
      bm.diam = fread(fid, 1, '*double');       % beam diameter (m)
      fclose(fid);

      npa  = size(bm.wf, 2);
    end                 % if ~prop_is_statesaved(bm)
  end                   % if (numel(save_state) > 0) & (save_state ~= 0)
end                     % function prop_state
