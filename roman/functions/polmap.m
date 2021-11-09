%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

% Function to compute the polarization aberrations by calling polab.m. 
% The mean polarization aberrations can be returned to mimic the effect of 
% wavefront sensing. 
%
% Translated by A.J. Riggs on 2017-06-26 from IDL code polmap.pro by John Krist.


% pro polmap, wavefront, polfile, pupil_diam_pix, condition, lambda_m_pol
% 
% ;-- wavefront: current wavefront structure
% ;-- polfile: rootname of file containing polarization coefficients
% ;-- pupil_diam_pix: diameter of pupil in pixels
% ;-- condition: polarization circumstance:
% ;--		-2: -45 deg in, Y out
% ;--		-1: -45 deg in, X out
% ;--		 1: +45 deg in, X out
% ;--		 2: +45 deg in, Y out
% ;--		 5: X polarization (mean of +/-45 deg in, X out)
% ;--		 6: Y polarization (mean of +/-45 deg in, X out)
% ;--		10: All polarizations (mean of +/-45 deg in, X&Y out)
% ;--	NOTE: the mean conditions (5,6,10) should only be used for sensing;
% ;--	contrast evaluation must be done by computing each in/out condition separately
% ;-- lambda_m_pol: (optional) wavelength in meters at which to compute polarization aberrations
% ;--	NOTE: this may be different from the wavelength at which one is running.
% ;--	For example, the HLC design is for a 523-578 nm bandpass, but one wants
% ;--	to know how the effect of polarization is at a shorter bandpass, assuming
% ;--	the HLC aberration sensitivity is the same in the shorter bandpass. The
% ;--	polarization aberrations IN WAVES at the shorter wavelength (e.g., 
% ;--	lambda_m_pol=460 nm) will be used at the current wavelength (e.g., 550 nm)   

function [wavefront] = polmap(wavefront, polfile, pupil_diam_pix, condition, lambda_m_pol)


n = prop_get_gridsize();
lambda_m = wavefront.wl;%prop_get_wavelength(wavefront);
% n = prop_get_gridsize()
% lambda_m = prop_get_wavelength(wavefront)

if ( numel(lambda_m_pol) == 0 ); lambda_m_pol = lambda_m; end
% if ( n_elements(lambda_m_pol) eq 0 ) then lambda_m_pol = lambda_m

if ( (condition == -2) || (condition == -1) || (condition == 1) || (condition == 2) )
	[amp, pha] = polab(polfile, lambda_m_pol, pupil_diam_pix, condition);
elseif ( condition == 5 ) %--For sensing of X polarization only
	[amp_m45_x, pha_m45_x] = polab(polfile, lambda_m_pol, pupil_diam_pix, -1);
	[amp_p45_x, pha_p45_x] = polab(polfile, lambda_m_pol, pupil_diam_pix, +1);
	amp = (amp_m45_x + amp_p45_x) / 2;
	pha = (pha_m45_x + pha_p45_x) / 2;
elseif ( condition == 6 ) %--For sensing of Y polarization only
	[amp_m45_y, pha_m45_y] = polab(polfile, lambda_m_pol, pupil_diam_pix, -2);
	[amp_p45_y, pha_p45_y] = polab(polfile, lambda_m_pol, pupil_diam_pix, +2); 
	amp = (amp_m45_y + amp_p45_y) / 2;
	pha = (pha_m45_y + pha_p45_y) / 2;
elseif ( condition == 10 ) %--For sensing of all polarizations only
	[amp_m45_x, pha_m45_x] = polab(polfile, lambda_m_pol, pupil_diam_pix, -1);
	[amp_p45_x, pha_p45_x] = polab(polfile, lambda_m_pol, pupil_diam_pix, +1);
	[amp_m45_y, pha_m45_y] = polab(polfile, lambda_m_pol, pupil_diam_pix, -2); 
	[amp_p45_y, pha_p45_y] = polab(polfile, lambda_m_pol, pupil_diam_pix, +2); 
	amp = (amp_m45_x + amp_p45_x + amp_m45_y + amp_p45_y) / 4;
	pha = (pha_m45_x + pha_p45_x + pha_m45_y + pha_p45_y) / 4;
else
	disp(['POLMAP: unmatched condition = ', num2str(condition)]);
	return
end

% if ( condition le 2 ) then begin
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, condition, amp, pha
% endif else if ( condition eq 5 ) then begin
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, -1, amp_m45_x, pha_m45_x
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, +1, amp_p45_x, pha_p45_x
% 	amp = (amp_m45_x + amp_p45_x) / 2
% 	pha = (pha_m45_x + pha_p45_x) / 2
% endif else if ( condition eq 6 ) then begin
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, -2, amp_m45_y, pha_m45_y
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, +2, amp_p45_y, pha_p45_y
% 	amp = (amp_m45_y + amp_p45_y) / 2
% 	pha = (pha_m45_y + pha_p45_y) / 2
% endif else if ( condition eq 10 ) then begin
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, -1, amp_m45_x, pha_m45_x
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, +1, amp_p45_x, pha_p45_x
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, -2, amp_m45_y, pha_m45_y
% 	polab, polfile, lambda_m_pol, pupil_diam_pix, +2, amp_p45_y, pha_p45_y
% 	amp = (amp_m45_x + amp_p45_x + amp_m45_y + amp_p45_y) / 4
% 	pha = (pha_m45_x + pha_p45_x + pha_m45_y + pha_p45_y) / 4
% endif else begin
% 	print, 'POLMAP: unmatched condition = ', condition
% 	stop
% endelse

%size(amp), size(wavefront.wf)
wavefront = prop_multiply(wavefront, custom_pad(amp,size(wavefront.wf,1)));

if ( lambda_m_pol ~= lambda_m ); pha = pha / lambda_m_pol * lambda_m; end

wavefront = prop_add_phase(wavefront, custom_pad(pha,size(wavefront.wf,1)));

end
% prop_multiply, wavefront, trim(amp,n) 
% 
% if ( lambda_pol_m ne lambda_m ) then pha = pha / lambda_pol_m * lambda_m
% 
% prop_add_phase, wavefront, trim(pha,n) 
% 
% return
% end