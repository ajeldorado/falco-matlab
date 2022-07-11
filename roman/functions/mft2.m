%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

%function out = mft2(in, dout, D, nout, direction, xoffset, yoffset, xc, yc)
%----------------------------------------------------------------------
% mtf:
%   Compute a matrix fourier transform.  Based on Soummer et al. 2007.
%   original format:
%   function out = CG_mft(in, nb, m, direction, xoffset, yoffset, xc, yc)
%   If the beam size is "D" pixels, then the output is sampled by "dx" cycles/D 
%   when m = dx * n_in * n_out / D
%   "xoffset" = offset_pix * dx
%
% Input parameters:
%    in : 2-D wavefront to transform
%   dout: sampling in lambda/D of output (if pupil-to-focus) or input (if focus-to-pupil)  
%      D: pupil size in pixels
%   nout: dimensions of output array (nout by nout)
%  direction : direction of transform (-1 or +1) 
%
% Optional input parameters:
%   xoffset, yoffset : offsets in output field in cycles/D
%
%  Returns:
%    2-D Fourier transform of input array.
%
%  Written by Dimitri Mawet (JPL)
%  March 2010
%  Copyright 2012 California Institute of Technology
%
%  Based on JEK version, 2/28/2019
%  Fixed 31/10/2020 by JEK: array centers were not begin computed with integer division
%------------------------------------------------------------------------
function out = mft2(in, dout, D, nout, direction, xoffset, yoffset, xc, yc)

if nargin<6  xoffset = 0.0; end
if nargin<7  yoffset = 0.0; end
if nargin<8  xc = 0.0; end
if nargin<9  yc = 0.0;end

nout = fix(nout); 
nin =  fix(size(in,1));

out= complex(zeros(nout,nout));

x=(((0:nin-1)-fix(nin/2)-xc))';
y=(((0:nin-1)-fix(nin/2)-yc))';

u=(((0:nout-1)-fix(nout/2)-xoffset/dout)*(dout/D))'; 
v=(((0:nout-1)-fix(nout/2)-yoffset/dout)*(dout/D))'; 

if direction == -1
	out= dout/D*exp((-1i*2*pi*u)*x')* in * exp((-1i*2*pi*y)*v'); 
elseif direction == 1
    	out= dout/D*exp((+1i*2*pi*u)*x')* in * exp((+1i*2*pi*y)*v');
end


return

