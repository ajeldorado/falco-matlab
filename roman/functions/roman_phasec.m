%  Copyright 2020 California Institute of Technology
% ------------------------------------------------------------------

%function [ wavefront, sampling_m ] =  roman_phasec(lambda_m, output_dim0, optval)
%
% PUBLIC RELEASE VERSION (no ITAR error maps)
%
% Phase C, Version 0.9, 15 June 2020, John Krist; beta-test version for mask evaluation
% Phase C, Version 0.9.1, 29 June 2020, John Krist; rotated SPC pupils & Lyot stops
% Phase C, Version 0.9.2, 7 July 2020, John Krist; updated HLC occulter
% Phase C, Version 1.0, 1.0.1, 15 Sept 2020, John Krist; updated HLC, replaced some error maps with measured ones, added thick lenses 
% Phase C, Version 1.0.2, 28 Oct 2020, John Krist: added optional mask array parameters (for JPL internal uses)
% Phase C, Version 1.1, 2 Dec 2020, John Krist: increased distance from TT Fold to FSM to move the
%     pupil 6 mm from the SPC mask (and Lyot stop) to represent telescope pupil defocus caused by
%     the TCA; increased the pupil imaging lens doublet separation by 0.2 mm to produce a good pupil image
%     on the detector without having to fudge the detector distance
% Version 1.2, 19 Feb 2021, John Krist: added code to subtract quadratic phase term in final image due
%     to conjugate pupil not at front focus of imaging lens (only applied when NOT using PIL or defocus
%     lenses); switched to MFT-based propagation to/from HLC FPM at high sampling to avoid interpolation
%     errors; added spc-mswc mode (provided by A.J. Riggs) for multi-star wavefront control modeling;
%     updated error maps; removed user-specfied lam0 and lambda0_m parameters (use only hard-coded values); 
%     changed default grid sizes, including getting rid of intermediate grid size changes; added
%     secondary mirror despace parameter; read in list of available HLC FPM files; removed 
%     propagation to fold 4 (not significant); added TO_PLANE option to thick lenses; added options
%     for using Band 1 with the spc-wide and spc-mswc; removed the output_field_rootname option;
%     added Phase C polarization aberrations from Jim McGuire; added HLCs for bands 2-3 and rotated SPC
% Version 1.2.1, 14 May 2021, John Krist: added option to use field_stop_array with hlc
% Version 1.2.2, 20 May 2021, John Krist: updated aperture sizes; use 2Kx2K grids for HLC when phase
%     retrieval lenses used; use OAP6 aperture when pinhole is specified
% Version 1.2.3, 15 July 2021, John Krist: added dm_sampling_m to optional parameters
% Version 1.2.4, 3 Aug 2021, John Krist: added HLC FPM files for band 3g
% Based on initial translation by Hanying Zhou of John Krist's Phase B IDL model; converted to Phase C by John Krist  


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wavefront,  sampling_m] = roman_phasec( lambda_m, output_dim0, optval )

% "output_dim" is used to specify the output dimension in pixels at the final image plane.
% The computational grid sizes are hardcoded for each coronagraph.

% mask design data path

data_dir = '/home/krist/afta/phasec/roman_phasec_v1.2.2/phasec_data';  % no trailing '/'

if ~exist('optval'); optval = 0; end
    
if ( isfield(optval,'data_dir') )
    data_dir = optval.data_dir;
    if(strcmp(data_dir(end),'/') || strcmp(data_dir(end),'\'))
        data_dir = data_dir(1:end-1); %--Remove the trailing slash for this function.
    end
end

map_dir = [data_dir '/maps/'];            % directory for surface maps
polfile = [data_dir '/pol/phasec_pol'];		% polarization aberration table rootname

pupil_array = 0;                 % 2D array containing pupil pattern (overrides default)
pupil_mask_array = 0;              % 2D array containing SPC pupil mask pattern (overrides default)
fpm_array = 0;                   % 2D array containing FPM mask pattern (overrides default)
fpm_mask = 0;			 % 2D array where 1=FPM pattern defined, 0=substrate
lyot_stop_array = 0;             % 2D array containing Lyot stop mask pattern (overrides default)
field_stop_array = 0;            % 2D array containing field stop mask pattern (overrides default)

cor_type = 'hlc';           	% 'hlc', 'spc-spec_band2', 'spc-spec_band3', 'spc-wide', 'none'
source_x_offset_mas = 0;	% source offset in milliarcsec (tilt applied at primary)
source_y_offset_mas = 0;
source_x_offset = 0;		% source offset in lambda0_m/D radians (tilt applied at primary)
source_y_offset = 0;
polaxis = 0;                	% polarization axis aberrations:
%     -2 = -45; in, Y out
%     -1 = -45; in, X out
%      1 = +45; in, X out
%      2 = +45; in, Y out
%      5 = mean of modes -1 & +1 (X channel polarizer)
%      6 = mean of modes -2 & +2 (Y channel polarizer)
%     10 = mean of all modes (no polarization filtering)
use_errors = 1;             	% use optical surface phase errors? 1 or 0
zindex = 0;                 	% array of Zernike polynomial indices
zval_m = 0;                 	% array of Zernike coefficients (meters RMS WFE)
sm_despace_m = 0;		% secondary mirror despace (meters)
use_pupil_defocus = 1;		% include pupil defocus
use_aperture = 0;           	% use apertures on all optics? 1 or 0
cgi_x_shift_pupdiam = 0;	% X,Y shift of wavefront at FSM (bulk displacement of CGI); normalized relative to pupil diameter
cgi_y_shift_pupdiam = 0;
cgi_x_shift_m = 0;		% X,Y shift of wavefront at FSM (bulk displacement of CGI) in meters
cgi_y_shift_m = 0;
end_at_fsm = 0;             	% end propagation after propagating to FSM (no FSM errors)
fsm_x_offset_mas = 0;		% offset source using FSM (milliarcsec)
fsm_y_offset_mas = 0;
fsm_x_offset = 0;           	% offset source using FSM (lambda0/D radians)
fsm_y_offset = 0;
focm_z_shift_m = 0;         	% offset (meters) of focus correction mirror (+ increases path length)
use_hlc_dm_patterns = 0;	% use Dwight-generated HLC default DM wavefront patterns? 1 or 0
use_dm1 = 0;                	% use DM1? 1 or 0
use_dm2 = 0;                	% use DM2? 1 or 0
dm_sampling_m = 0.9906e-3;  	% actuator spacing in meters;
dm1_m = zeros(48,48);
dm1_xc_act = 23.5;          	% for 48x48 DM(wavefront centered at actuator intersections: (0,0) = 1st actuator center
dm1_yc_act = 23.5;
dm1_xtilt_deg = 0;  		% tilt around X axis
dm1_ytilt_deg = 9.65;		% DM tilt in deg 
dm1_ztilt_deg = 0;
dm2_m = zeros(48,48);
dm2_xc_act = 23.5;
dm2_yc_act = 23.5;
dm2_xtilt_deg = 0;
dm2_ytilt_deg = 9.65;
dm2_ztilt_deg = 0;
hlc_dm1_file = '';
hlc_dm2_file = '';
use_pupil_mask = 1;		% (SPC only) Use SPC pupil mask (0 or 1)
mask_x_shift_pupdiam = 0;	% X,Y shift of shaped pupil mask; normalized relative to pupil diameter
mask_y_shift_pupdiam = 0;
mask_x_shift_m = 0;         	% X,Y shift of shaped pupil mask; 
mask_y_shift_m = 0;
use_fpm = 1;                	% use occulter? 1 or 0
fpm_x_offset = 0;           	% FPM x,y offset in lambda0/D
fpm_y_offset = 0;
fpm_x_offset_m = 0;          	% FPM x,y offset in meters
fpm_y_offset_m = 0;
fpm_z_shift_m = 0;          	% FPM offset in meters along optical axis (+ = away from prior optics)
pinhole_diam_m = 0;         	% pinhole diameter in meters at FPM
end_at_fpm_exit_pupil = 0;	% return field at FPM exit pupil?
use_lyot_stop = 1;          	% use Lyot stop? 1 or 0
lyot_x_shift_pupdiam = 0;	% X,Y shift of Lyot stop mask; normalized relative to pupil diameter
lyot_y_shift_pupdiam = 0;
lyot_x_shift_m = 0;		% X,Y shift of Lyot stop mask in meters
lyot_y_shift_m = 0;
use_field_stop = 1;         	% use field stop (HLC)? 1 or 0
field_stop_radius_lam0 = 0;	% field stop radius in lambda0/D (HLC || SPC-wide mask only)
field_stop_x_offset = 0;        % field stop offset in lambda0/D
field_stop_y_offset = 0;
field_stop_x_offset_m = 0;      % field stop offset in meters
field_stop_y_offset_m = 0;
use_pupil_lens = 0;         	% use pupil imaging lens?
use_defocus_lens = 0;		% use defocusing lens? Options are 1, 2, 3, 4, corresponding to +18.0, +9.0, -4.0, -8.0 waves P-V @ 550 nm
end_at_exit_pupil = 0;		% return exit pupil corresponding to final image plane
final_sampling_lam0 = 0;	% final sampling in lambda0/D
final_sampling_m = 0;       	% final sampling in meters (overrides final_sampling_lam0)
output_dim = output_dim0;	% dimension of output grid in pixels (overrides output_dim0)

if exist('optval','var')==1 
    if ( isfield(optval,'cor_type') );       cor_type = optval.cor_type; end
    if ( isfield(optval,'use_fpm') );        use_fpm = optval.use_fpm; end
    if ( isfield(optval,'pupil_array') );    pupil_array = optval.pupil_array; end
    if ( isfield(optval,'pupil_mask_array') ); pupil_mask_array = optval.pupil_mask_array; end
    if ( isfield(optval,'fpm_array') )
	fpm_array = optval.fpm_array;
	fpm_array_sampling_m = optval.fpm_array_sampling_m;
    end
    if ( isfield(optval,'lyot_stop_array') ); lyot_stop_array = optval.lyot_stop_array; end
    if ( isfield(optval,'field_stop_array') ) 
	field_stop_array = optval.field_stop_array;
	field_stop_array_sampling_m = optval.field_stop_array_sampling_m;
    end
    if ( isfield(optval,'use_pupil_lens') );         use_pupil_lens = optval.use_pupil_lens; end
    if ( isfield(optval,'use_defocus_lens') );       use_defocus_lens = optval.use_defocus_lens; end
end

is_hlc = 0;
is_spc = 0;

if  ( strfind(cor_type,'hlc') )
    is_hlc = 1;
    if ( strcmp(cor_type,'hlc') || strcmp(cor_type,'hlc_band1') )
        file_directory = [data_dir '/hlc_20190210b/'];         % must have trailing "/"
        lambda0_m = 0.575e-6;
        hlc_dm1_file = [file_directory 'hlc_dm1.fits'];
        hlc_dm2_file = [file_directory 'hlc_dm2.fits'];
    elseif ( strcmp(cor_type,'hlc_band2') )
        file_directory = [data_dir '/hlc_20200617c_band2/'];
        lambda0_m = 0.660e-6;
    elseif ( strcmp(cor_type,'hlc_band3') )
        file_directory = [data_dir '/hlc_20200614b_band3/'];
        lambda0_m = 0.730e-6;
    elseif ( strcmp(cor_type,'hlc_band4') )
        file_directory = [data_dir '/hlc_20200609b_band4/'];
        lambda0_m = 0.825e-6;
    else 
        disp('Unsupported HLC mode');
        return
    end
    pupil_diam_pix = 309.0;	% Y pupil diameter in pixels
    pupil_file = [file_directory  'pupil.fits'];
    if use_fpm ~= 0
    		lam_um = lambda_m * 1e6;
                % find nearest available FPM wavelength that matches specified wavelength
		file = fopen( [file_directory 'fpm_files.txt'] );
		fpm_nlam = fscanf( file, '%d', 1 );
		fpm_lam_um = fscanf( file,'%f', fpm_nlam );
		fpm_lams = textscan( file,'%s', fpm_nlam );
                diff = abs(fpm_lam_um - lam_um);
                [mindiff, w] = min(diff);
                if mindiff > 0.0001
                        disp('Error in roman_phasec: requested wavelength not within 0.1 nm of nearest available FPM wavelength')
                        error(['   requested (um) = ', num2str(lam_um), '  closest available (um) = ', num2str(fpm_lam_um(w))])
                end
                fpm_rootname = [file_directory fpm_lams{1}{w}];
		% read in FPM sampling and reference wavelength
                h = fitsinfo( [fpm_rootname, 'real.fits'] );
		htab = cell2table(h.PrimaryData.Keywords);
		harr = table2array(htab(:,1));
		idx = find(ismember(harr,'FPMLAM0M'));
		fpm_lam0_m = cell2mat(h.PrimaryData.Keywords(idx,2));
		idx = find(ismember(harr,'FPMDX'));
		fpm_sampling_lam0divD = cell2mat(h.PrimaryData.Keywords(idx,2));
		% read in FPM 
                r = fitsread( [fpm_rootname, 'real.fits'] );
                i = fitsread( [fpm_rootname, 'imag.fits'] );
                fpm_mask = fix(r ~= r(1,1));
                fpm_array = complex(r,i);
    end
    lyot_stop_file = [file_directory  'lyot_rotated.fits'];
    field_stop_radius_lam0 = 9.7;
    if use_pupil_lens ~= 0 || use_defocus_lens ~= 0
	n = 2048; 
    else
	n = 1024;
    end
    n_mft = 1400;
elseif ( strcmp(cor_type,'spc-spec') || strcmp(cor_type,'spc-spec_band2') || strcmp(cor_type,'spc-spec_band3') ) 
    is_spc = 1;
    file_dir = [data_dir '/spc_20200617_spec/'];       % must have trailing "/"
    pupil_diam_pix = 1000;	% Y axis pupil diameter in pixels
    pupil_file = [file_dir  'pupil_SPC-20200617_1000.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200617_1000_rounded9_rotated.fits'];
    fpm_sampling_lam0divD = 0.05; 	% sampling in lambda0/D of FPM mask
    fpm_file = [file_dir 'fpm_0.05lamD.fits'];
    if strcmp(cor_type, 'spc-spec_band2')
	fpm_lam0_m = 0.66e-6;	% central wavelength of bandpass
    else
	fpm_lam0_m = 0.73e-6;
    end 
    lambda0_m = fpm_lam0_m;
    lyot_stop_file = [file_dir  'LS_SPC-20200617_1000.fits'];
    if use_pupil_lens ~= 0 || use_defocus_lens ~= 0
	n = 4096; 
    else
	n = 2048;
    end
    n_mft = 1400;
elseif ( strcmp(cor_type,'spc-spec_rotated') )
    is_spc = 1;
    file_dir = [data_dir '/spc_20200628_specrot/'];    % must have trailing "/"
    pupil_diam_pix = 1000;	% Y axis pupil diameter in pixels
    pupil_file = [file_dir  'pupil_SPC-20200628_1000.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200628_1000_derotated.fits'];
    fpm_sampling_lam0divD = 0.05; 	% sampling in lambda0/D of FPM mask
    fpm_file = [file_dir 'FPM_SPC-20200628_res20.fits'];
    fpm_lam0_m = 0.73e-6;
    lambda0_m = fpm_lam0_m;
    lyot_stop_file = [file_dir  'LS_SPC-20200628_1000.fits'];
    if use_pupil_lens ~= 0 || use_defocus_lens ~= 0
	n = 4096; 
    else
	n = 2048;
    end
    n_mft = 1400;
elseif ( strcmp(cor_type,'spc-wide') || strcmp(cor_type,'spc-wide_band4') || strcmp(cor_type,'spc-wide_band1') )
    is_spc = 1;
    file_dir = [data_dir '/spc_20200610_wfov/'];       % must have trailing "/"
    pupil_diam_pix = 1000.0;      % Y axis pupil diameter in pixels
    pupil_file = [file_dir  'pupil_SPC-20200610_1000.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200610_1000_rounded9_gray_rotated.fits'];
    fpm_sampling_lam0divD = 0.1;       % sampling in lambda0/D of FPM mask
    fpm_file = [file_dir   'FPM_SPC-20200610_0.1_lamc_div_D.fits'];
    if ( strcmp(cor_type,'spc-wide_band1') )
    	fpm_lam0_m = 0.575e-6;	% central wavelength of bandpass
    else
    	fpm_lam0_m = 0.825e-6;
    end
    lambda0_m = fpm_lam0_m;
    lyot_stop_file = [file_dir  'LS_SPC-20200610_1000.fits'];
    if use_pupil_lens ~= 0 || use_defocus_lens ~= 0
        n = 4096; 
    else
        n = 2048;
    end
    n_mft = 1400;
elseif ( strcmp(cor_type,'spc-mswc') || strcmp(cor_type,'spc-mswc_band4') || strcmp(cor_type,'spc-mswc_band1') )
    is_spc = 1;
    file_dir = [data_dir '/spc_20200623_mswc/'];       % must have trailing "/"
    pupil_diam_pix = 982.0;      % Y axis pupil diameter in pixels
    pupil_file = [file_dir  'pupil_SPC-20200623_982_rotated.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200623_982_rounded9_gray_rotated.fits'];
    fpm_sampling_lam0divD = 0.1;       % sampling in lambda0/D of FPM mask
    fpm_file = [file_dir 'FPM_SPC-20200623_0.1_lamc_div_D.fits'];
    if ( strcmp(cor_type,'spc-mswc_band1') )
    	fpm_lam0_m = 0.575e-6;
    else
    	fpm_lam0_m = 0.825e-6;
    end
    lambda0_m = fpm_lam0_m;
    lyot_stop_file = [file_dir 'LS_SPC-20200623_982_rotated.fits'];
    if use_pupil_lens ~= 0 || use_defocus_lens ~= 0
        n = 4096; 
    else
        n = 2048;
    end
    n_mft = 1400;
elseif ( strcmp(cor_type, 'none') )
    lambda0_m = 0.575e-6;
    pupil_diam_pix = 309.0;
    use_fpm = 0;
    use_lyot_stop = 0;
    use_field_stop =0;
    use_errors = 0;
    n = 1024;
else
    disp([ 'cor_type = ' cor_type ' is not supported'])
    return
end


if (exist('optval','var')==1) 
    if ( isfield(optval,'sm_despace_m') );   sm_despace_m = optval.sm_despace_m; end
    mas_per_lamD = lambda0_m * 360.0 * 3600.0 / (2 * pi * 2.363) * 1000;	% mas per lambda0/D
    if ( isfield(optval,'source_x_offset') && isfield(optval,'source_x_offset_mas') )
	disp('ERROR: can only specify source_x_offset or source_x_offset_mas, not both');
	return
    end
    if ( isfield(optval,'source_x_offset') );   source_x_offset = optval.source_x_offset; end
    if ( isfield(optval,'source_x_offset_mas') );   source_x_offset = optval.source_x_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'source_y_offset') );   source_y_offset = optval.source_y_offset;end
    if ( isfield(optval,'source_y_offset_mas') );   source_y_offset = optval.source_y_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'use_aperture') );   use_aperture = optval.use_aperture;end
    if ( isfield(optval,'use_errors') );     use_errors = optval.use_errors;end
    if ( isfield(optval,'use_pupil_defocus') );     use_pupil_defocus = optval.use_pupil_defocus;end
    if ( isfield(optval,'polaxis') );        polaxis = optval.polaxis; end
    if ( isfield(optval,'zindex') );         zindex = optval.zindex;end
    if ( isfield(optval,'zval_m') );         zval_m = optval.zval_m; end
    if ( isfield(optval,'end_at_fsm') );     end_at_fsm = optval.end_at_fsm; end
    if ( isfield(optval,'cgi_x_shift_pupdiam') && isfield(optval,'cgi_x_shift_m') )
	disp('ERROR: can only specify cgi_x_shift_pupdiam or cgi_x_shift_m, not both');
	return
    end
    if ( isfield(optval,'cgi_x_shift_pupdiam') );   cgi_x_shift_pupdiam = optval.cgi_x_shift_pupdiam; end
    if ( isfield(optval,'cgi_x_shift_m') );  cgi_x_shift_m = optval.cgi_x_shift_m; end
    if ( isfield(optval,'cgi_y_shift_pupdiam') && isfield(optval,'cgi_y_shift_m') )
	disp('ERROR: can only specify cgi_y_shift_pupdiam or cgi_y_shift_m, not both');
	return
    end
    if ( isfield(optval,'cgi_y_shift_pupdiam') );   cgi_y_shift_pupdiam = optval.cgi_y_shift_pupdiam; end
    if ( isfield(optval,'cgi_y_shift_m') );  cgi_y_shift_m = optval.cgi_y_shift_m; end
    if ( isfield(optval,'fsm_x_offset') && isfield(optval,'fsm_x_offset_mas') )
	disp('ERROR: can only specify fsm_x_offset or fsm_x_offset_mas, not both');
	return
    end
    if ( isfield(optval,'fsm_x_offset') );   fsm_x_offset = optval.fsm_x_offset; end
    if ( isfield(optval,'fsm_x_offset_mas') );      fsm_x_offset = optval.fsm_x_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'fsm_y_offset') && isfield(optval,'fsm_y_offset_mas') )
	disp('ERROR: can only specify fsm_y_offset or fsm_y_offset_mas, not both');
	return
    end
    if ( isfield(optval,'fsm_y_offset') );   fsm_y_offset = optval.fsm_y_offset; end
    if ( isfield(optval,'fsm_y_offset_mas') );      fsm_y_offset = optval.fsm_y_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'focm_z_shift_m') );        focm_z_shift_m = optval.focm_z_shift_m; end
    if ( isfield(optval,'use_hlc_dm_patterns') );   use_hlc_dm_patterns = optval.use_hlc_dm_patterns; end
    if ( isfield(optval,'dm_sampling_m') );  dm_sampling_m = optval.dm_sampling_m; end
    if ( isfield(optval,'use_dm1') );     	use_dm1 = optval.use_dm1; end
    if ( isfield(optval,'dm1_m') );          dm1_m = optval.dm1_m; end
    if ( isfield(optval,'dm1_xc_act') );     dm1_xc_act = optval.dm1_xc_act; end
    if ( isfield(optval,'dm1_yc_act') );     dm1_yc_act = optval.dm1_yc_act; end
    if ( isfield(optval,'dm1_xtilt_deg') );  dm1_xtilt_deg = optval.dm1_xtilt_deg; end
    if ( isfield(optval,'dm1_ytilt_deg') );  dm1_ytilt_deg = optval.dm1_ytilt_deg; end
    if ( isfield(optval,'dm1_ztilt_deg') );  dm1_ztilt_deg = optval.dm1_ztilt_deg; end
    if ( isfield(optval,'use_dm2') );        use_dm2 = optval.use_dm2; end
    if ( isfield(optval,'dm2_m') );          dm2_m = optval.dm2_m; end
    if ( isfield(optval,'dm2_xc_act') );     dm2_xc_act = optval.dm2_xc_act; end
    if ( isfield(optval,'dm2_yc_act') );     dm2_yc_act = optval.dm2_yc_act; end
    if ( isfield(optval,'dm2_xtilt_deg') );  dm2_xtilt_deg = optval.dm2_xtilt_deg; end
    if ( isfield(optval,'dm2_ytilt_deg') );  dm2_ytilt_deg = optval.dm2_ytilt_deg; end
    if ( isfield(optval,'dm2_ztilt_deg') );  dm2_ztilt_deg = optval.dm2_ztilt_deg; end
    if ( isfield(optval,'use_pupil_mask') );  use_pupil_mask = optval.use_pupil_mask; end
    if ( isfield(optval,'mask_x_shift_pupdiam') && isfield(optval,'mask_x_shift_m') )
	disp('ERROR: can only specify mask_x_shift_pupdiam or mask_x_shift_m, not both');
	return
    end
    if ( isfield(optval,'mask_x_shift_pupdiam') );   mask_x_shift_pupdiam = optval.mask_x_shift_pupdiam; end
    if ( isfield(optval,'mask_x_shift_m') );     mask_x_shift_m = optval.mask_x_shift_m; end
    if ( isfield(optval,'mask_y_shift_pupdiam') && isfield(optval,'mask_y_shift_m') )
	disp('ERROR: can only specify mask_y_shift_pupdiam or mask_y_shift_m, not both');
	return
    end
    if ( isfield(optval,'mask_y_shift_pupdiam') );   mask_y_shift_pupdiam = optval.mask_y_shift_pupdiam; end
    if ( isfield(optval,'mask_y_shift_m') );     mask_y_shift_m = optval.mask_y_shift_m; end
    if ( isfield(optval,'fpm_x_offset') && isfield(optval,'fpm_x_offset_m') )
	disp('ERROR: can only specify fpm_x_offset or fpm_x_offset_m, not both');
	return
    end
    if ( isfield(optval,'fpm_x_offset') );       fpm_x_offset = optval.fpm_x_offset; end
    if ( isfield(optval,'fpm_x_offset_m') );      fpm_x_offset_m = optval.fpm_x_offset_m; end
    if ( isfield(optval,'fpm_y_offset') && isfield(optval,'fpm_y_offset_m') )
	disp('ERROR: can only specify fpm_y_offset or fpm_y_offset_m, not both');
	return
    end
    if ( isfield(optval,'fpm_y_offset') );       fpm_y_offset = optval.fpm_y_offset; end
    if ( isfield(optval,'fpm_y_offset_m') );      fpm_y_offset_m = optval.fpm_y_offset_m; end
    if ( isfield(optval,'fpm_z_shift_m') );      fpm_z_shift_m = optval.fpm_z_shift_m; end
    if ( isfield(optval,'pinhole_diam_m') );     pinhole_diam_m = optval.pinhole_diam_m; end
    if ( isfield(optval,'end_at_fpm_exit_pupil') );  end_at_fpm_exit_pupil = optval.end_at_fpm_exit_pupil; end
    if ( isfield(optval,'use_lyot_stop') );          use_lyot_stop = optval.use_lyot_stop; end
    if ( isfield(optval,'lyot_x_shift_pupdiam') && isfield(optval,'lyot_x_shift_m') )
	disp('ERROR: can only specify lyot_x_shift_pupdiam or lyot_x_shift_m, not both');
	return
    end
    if ( isfield(optval,'lyot_x_shift_pupdiam') );   lyot_x_shift_pupdiam = optval.lyot_x_shift_pupdiam; end
    if ( isfield(optval,'lyot_x_shift_m') );         lyot_x_shift_m = optval.lyot_x_shift_m; end
    if ( isfield(optval,'lyot_y_shift_pupdiam') && isfield(optval,'lyot_y_shift_m') )
	disp('ERROR: can only specify lyot_y_shift_pupdiam or lyot_y_shift_m, not both');
	return
    end
    if ( isfield(optval,'lyot_y_shift_pupdiam') );   lyot_y_shift_pupdiam = optval.lyot_y_shift_pupdiam; end
    if ( isfield(optval,'lyot_y_shift_m') );         lyot_y_shift_m = optval.lyot_y_shift_m; end
    if ( isfield(optval,'use_field_stop') );         use_field_stop = optval.use_field_stop; end
    if ( isfield(optval,'field_stop_radius_lam0') ); field_stop_radius_lam0 = optval.field_stop_radius_lam0; end
    if ( isfield(optval,'field_stop_x_offset') && isfield(optval,'field_stop_x_offset_m') )
	disp('ERROR: can only specify field_stop_x_offset or field_stop_x_offset_m, not both');
	return
    end
    if ( isfield(optval,'field_stop_x_offset') );  field_stop_x_offset = optval.field_stop_x_offset; end
    if ( isfield(optval,'field_stop_x_offset_m') );  field_stop_x_offset_m = optval.field_stop_x_offset_m; end
    if ( isfield(optval,'field_stop_y_offset') && isfield(optval,'field_stop_y_offset_m') )
	disp('ERROR: can only specify field_stop_y_offset or field_stop_y_offset_m, not both');
	return
    end
    if ( isfield(optval,'field_stop_y_offset') );  field_stop_y_offset = optval.field_stop_y_offset; end
    if ( isfield(optval,'field_stop_y_offset_m') );  field_stop_y_offset_m = optval.field_stop_y_offset_m; end
    if ( isfield(optval,'end_at_exit_pupil') );      end_at_exit_pupil = optval.end_at_exit_pupil; end
    if ( isfield(optval,'final_sampling_lam0') && isfield(optval,'final_sampling_m') )
	disp('ERROR: can only specify final_sampling_lam0 or final_sampling_m, not both');
	return
    end
    if ( isfield(optval,'final_sampling_lam0') );    final_sampling_lam0 = optval.final_sampling_lam0; end
    if ( isfield(optval,'final_sampling_m') );       final_sampling_m = optval.final_sampling_m; end
    if ( isfield(optval,'output_dim') );             output_dim = optval.output_dim; end
end

if use_hlc_dm_patterns && ~strcmp(cor_type,'hlc') && ~strcmp(cor_type,'hlc_band1')
    disp('ERROR: Can only utilize use_hlc_dm_patterns with Band 1 HLC');
    return
end

diam = 2.363114;
  fl_pri = 2.838279206904720;
d_pri_sec = 2.285150508110035 + sm_despace_m;
  fl_sec = -0.654200796568004;
  diam_sec = 0.58166;
d_sec_pomafold = 2.993753469304728 + sm_despace_m;
  diam_pomafold = 0.09;
d_pomafold_m3 = 1.680935841598811;
  fl_m3 = 0.430216463069001;
  diam_m3 = 0.2;
  d_m3_pupil = 0.469156807765176;
d_m3_m4 = 0.943514749358944;
  fl_m4 = 0.116239114833590;
  diam_m4 = 0.07;
d_m4_m5 = 0.429145636743193;
  fl_m5 = 0.198821518772608;
  d_m5_pupil = 0.716529242927776;
  diam_m5 = 0.07;
d_m5_ttfold = 0.351125431220770;
  diam_ttfold = 0.06;
d_ttfold_fsm = d_m5_pupil - d_m5_ttfold;
  if ( use_pupil_defocus )
	d_ttfold_fsm = d_ttfold_fsm + 0.033609;	% 33.6 mm to put pupil 6 mm from SPC mask
  end
  diam_fsm = 0.0508;
d_fsm_oap1 = 0.354826767220001;
  fl_oap1 = 0.503331895563883;
  diam_oap1 = 0.060;
d_oap1_focm = 0.768029932093727 + focm_z_shift_m;
  diam_focm = 0.035;
d_focm_oap2 = 0.314507535543064 + focm_z_shift_m;
  fl_oap2 = 0.579205571254990;
  diam_oap2 = 0.060;
d_oap2_dm1 = 0.775857408587825;
d_dm1_dm2 = 1.0;
d_dm2_oap3 = 0.394833855161549;
  fl_oap3 = 1.217276467668519;
  diam_oap3 = 0.06;
d_oap3_fold3 = 0.505329955078121;
  diam_fold3 = 0.06;
d_fold3_oap4 = 1.158897671642761;
  fl_oap4 = 0.446951159052363;
  diam_oap4 = 0.06;
d_oap4_pupilmask = 0.423013568764728;
d_pupilmask_oap5 = 0.408810704327559;
  fl_oap5 = 0.548189354706822;
  diam_oap5 = 0.06;
d_oap5_fpm = fl_oap5;                    % to front of FPM 
  fpm_thickness = 0.006363747896388863;     % account for FPM thickness (inclination included)
  fpm_index = glass_index('SILICA',lambda_m,data_dir);
d_fpm_oap6 = fpm_thickness / fpm_index + 0.543766629917668;     % from front of FPM
  fl_oap6 = d_fpm_oap6;
  diam_oap6 = 0.054;
d_oap6_lyotstop = 0.687476361491529;
d_oap6_exitpupil = d_oap6_lyotstop - 6e-3;
d_lyotstop_oap7 = 0.401748561745987;
  fl_oap7 = 0.708251420923810;
  diam_oap7 = 0.054;
d_oap7_fieldstop = fl_oap7; 
d_fieldstop_oap8 = 0.210985170345932 * 0.997651;
  fl_oap8 = d_fieldstop_oap8;
  diam_oap8 = 0.026;
  d_oap8_pupil = 0.237561587674008;
  d_pupil_filter = 0.130;
d_oap8_filter = d_oap8_pupil + d_pupil_filter;   % to front of filter
  diam_filter = 0.009;
  filter_thickness = 0.004016105782012525;      % account for filter thickness (inclination included)
  filter_index = glass_index('SILICA',lambda_m,data_dir);
d_filter_lens = filter_thickness / filter_index + 0.210581269256657095;   % from front of filter
  diam_lens = 0.0104;
d_lens_fold4 = 0.202432155667761;
  if use_pupil_lens ~= 0
        d_lens_fold4 = d_lens_fold4 - 0.0002;   % from back of pupil imaging lens 
  elseif use_defocus_lens ~= 0
        d_lens_fold4 = d_lens_fold4 + 0.001;    % doublet is 1 mm longer than singlet, so make up for it
  end
  diam_fold4 = 0.036;
d_fold4_image = 0.050000152941020161;


wavefront = prop_begin(diam,lambda_m, n,'beam_diam_fraction', pupil_diam_pix/n);
if ( strcmp(cor_type,'none') )
   wavefront = prop_circular_aperture(wavefront, 1.0, 'NORM');
else
   if ( isscalar(pupil_array) )
   	pupil = fitsread( pupil_file);
   	wavefront = prop_multiply(wavefront, custom_pad(pupil,n));
   	clear pupil;
   else
	wavefront = prop_multiply(wavefront, custom_pad(pupil_array,n));
   end
end
if ( polaxis ~= 0 );  wavefront =  polmap(wavefront, polfile, pupil_diam_pix, polaxis,lambda_m); end
wavefront = prop_define_entrance(wavefront);
wavefront = prop_lens(wavefront, fl_pri, 'primary');
if ( source_x_offset ~= 0 || source_y_offset ~= 0)
    % compute tilted wavefront to offset source by source_x_offset,source_y_offset lambda0_m/D
    x  = repmat( ((1:n)-n/2-1)/(pupil_diam_pix/2),n,1);
    tilt = - 1i*pi*(x*source_x_offset+x'*source_y_offset) * lambda0_m / lambda_m;%
    wavefront = prop_multiply(wavefront, exp(tilt));
    clear x tilt
end
if ( sum(zindex) ~= 0 );   wavefront = prop_zernikes(wavefront, zindex, zval_m); end
if ( use_errors ); wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_PRIMARY_synthetic_phase_error_V1.0.fits'], 'wavefront'); end

wavefront = prop_propagate(wavefront, d_pri_sec, 'surface_name','secondary');
wavefront = prop_lens(wavefront, fl_sec, 'secondary');
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_SECONDARY_synthetic_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_sec/2 );end

wavefront = prop_propagate(wavefront, d_sec_pomafold, 'surface_name','POMA FOLD');
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_POMAFOLD_measured_phase_error_V1.1.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture( wavefront, diam_pomafold/2 );end

wavefront = prop_propagate(wavefront, d_pomafold_m3, 'surface_name','M3');
wavefront = prop_lens(wavefront, fl_m3);
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_M3_measured_phase_error_V1.1.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_m3/2);end

wavefront = prop_propagate(wavefront, d_m3_m4, 'surface_name','M4');
wavefront = prop_lens(wavefront, fl_m4);
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_M4_measured_phase_error_V1.1.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_m4/2);end

wavefront = prop_propagate(wavefront, d_m4_m5, 'surface_name','M5');
wavefront = prop_lens(wavefront, fl_m5);
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_M5_measured_phase_error_V1.1.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_m5/2 ); end

wavefront =prop_propagate(wavefront, d_m5_ttfold, 'surface_name','TT FOLD');
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'roman_phasec_TTFOLD_measured_phase_error_V1.1.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_ttfold/2); end

wavefront = prop_propagate(wavefront, d_ttfold_fsm, 'surface_name','FSM');
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_LOWORDER_phase_error_V2.0.fits'], 'wavefront'); end
if ( end_at_fsm)
    [wavefront, sampling_m] = prop_end(wavefront, 'NOABS');
    return
end
if ( cgi_x_shift_pupdiam ~= 0 || cgi_y_shift_pupdiam ~= 0 || cgi_x_shift_m ~= 0 || cgi_y_shift_m ~= 0  )      % bulk coronagraph pupil shear
    % shears are normalized to pupil diameter
    % FFT the field, apply a tilt, FFT back
    x  = repmat( ((1:n)-n/2-1)/(pupil_diam_pix/2),n,1);
    if (cgi_x_shift_pupdiam ~= 0 || cgi_y_shift_pupdiam ~= 0)   % input pup beam shear in % of D
        tilt =  - 1i*pi*(x*cgi_x_shift_pupdiam + x'*cgi_y_shift_pupdiam)*pupil_diam_pix*pupil_diam_pix/n;
    elseif  (cgi_x_shift_m ~= 0 || cgi_y_shift_m ~= 0)          % input pup beam shear in abs value
        tilt =  - 1i*pi*(x*cgi_x_shift_m + x'*cgi_y_shift_m)/prop_get_sampling(wavefront)*pupil_diam_pix/n;
    end
    wavefront.wf = fft2( ifft2(wavefront.wf).*ifftshift(exp(tilt)) ); % fft the field, apply tilt, fft back
    clear x tilt 
end
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_FSM_FLIGHT_measured_coated_phase_error_V2.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_fsm/2);end
if ( fsm_x_offset ~= 0 || fsm_y_offset ~= 0)
    % compute tilted wavefront to offset source by fsm_x_offset,fsm_y_offset lambda0_m/D
    x  = repmat( ((1:n)-n/2-1)/(pupil_diam_pix/2),n,1);
    tilt =  1i*pi*(x*fsm_x_offset + x'*fsm_y_offset)* lambda0_m / lambda_m;
    wavefront = prop_multiply(wavefront, exp(tilt));
    clear x tilt 
end

wavefront =prop_propagate(wavefront, d_fsm_oap1, 'surface_name','OAP1');
wavefront = prop_lens(wavefront, fl_oap1);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP1_phase_error_V3.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap1/2);end

wavefront = prop_propagate(wavefront, d_oap1_focm, 'surface_name','FOCM');
if ( use_errors );  wavefront =  prop_errormap(wavefront,[map_dir 'roman_phasec_FCM_EDU_measured_coated_phase_error_V2.0.fits'], 'wavefront');end
if ( use_aperture );  wavefront =  prop_circular_aperture(wavefront, diam_focm/2);end

wavefront =prop_propagate(wavefront, d_focm_oap2, 'surface_name','OAP2');
wavefront = prop_lens(wavefront, fl_oap2);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP2_phase_error_V3.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap2/2 );end

wavefront = prop_propagate(wavefront, d_oap2_dm1, 'surface_name','DM1');
if ( is_hlc && use_hlc_dm_patterns && ~strcmp(hlc_dm1_file,'') ) 
    hlc_dm1 = fitsread( hlc_dm1_file );
    dm1 = dm1_m + hlc_dm1;
    use_dm1 = 1;
    hlc_dm2 = fitsread( hlc_dm2_file );
    dm2 = dm2_m + hlc_dm2;
    use_dm2 = 1;
else
    dm1 = dm1_m;
    dm2 = dm2_m;
end
if ( use_dm1 ); wavefront = prop_dm(wavefront, dm1, dm1_xc_act, dm1_yc_act, dm_sampling_m, 'xtilt',dm1_xtilt_deg, 'ytilt',dm1_ytilt_deg, 'ztilt',dm1_ztilt_deg); end
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_DM1_phase_error_V1.0.fits'], 'wavefront');end

wavefront = prop_propagate(wavefront, d_dm1_dm2, 'surface_name','DM2');
if ( use_dm2 );   wavefront = prop_dm(wavefront, dm2, dm2_xc_act, dm2_yc_act, dm_sampling_m, 'ztilt',dm2_xtilt_deg, 'ytilt',dm2_ytilt_deg, 'ztilt',dm2_ztilt_deg);end
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_DM2_phase_error_V1.0.fits'], 'wavefront');end
if ( is_hlc ) 
    dm2mask = fitsread( [file_directory 'dm2mask.fits']);
    wavefront = prop_multiply(wavefront, custom_pad(dm2mask, n));
    clear dm2mask 
end

wavefront = prop_propagate(wavefront, d_dm2_oap3, 'surface_name','OAP3');
wavefront = prop_lens(wavefront, fl_oap3);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP3_phase_error_V3.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap3/2);end

wavefront = prop_propagate(wavefront, d_oap3_fold3, 'surface_name','FOLD_3');
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_FOLD3_FLIGHT_measured_coated_phase_error_V2.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_fold3/2);end

wavefront = prop_propagate(wavefront, d_fold3_oap4, 'surface_name','OAP4');
wavefront = prop_lens(wavefront, fl_oap4);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP4_phase_error_V3.0.fits'], 'wavefront');  end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap4/2);  end

wavefront = prop_propagate(wavefront, d_oap4_pupilmask,'surface_name', 'PUPIL_MASK');	% flat/reflective shaped pupil
if ( (is_spc == 1 || isscalar(pupil_mask_array) == 0) && use_pupil_mask ~= 0 ) 
    if ( isscalar(pupil_mask_array) )
    	pupil_mask = custom_pad(fitsread(pupil_mask_file),n);
    else
	pupil_mask = custom_pad(pupil_mask_array,n);
    end
    if ( mask_x_shift_pupdiam ~=0 || mask_y_shift_pupdiam ~=0  || mask_x_shift_m ~=0 || mask_y_shift_m ~=0)
        %shift SP mask by FFTing it, applying tilt, and FFTing back
        x = repmat( ((1:n)-n/2-1) / (pupil_diam_pix/2), n, 1 );
        if ( mask_x_shift_pupdiam ~=0 || mask_y_shift_pupdiam ~=0) % shifts are in % of D
            tilt =  -1i*pi*(x*mask_x_shift_pupdiam +x'*mask_y_shift_pupdiam) * pupil_diam_pix * pupil_diam_pix/n;
        elseif (mask_x_shift_m ~=0 || mask_y_shift_m ~=0)
            tilt =  -1i*pi*(x*mask_x_shift_m +x'*mask_y_shift_m) / prop_get_sampling(wavefront) * pupil_diam_pix/n;
        end
        pupil_mask = ifft2( ifftshift(pupil_mask) ) .* ifftshift( exp(tilt) );
        pupil_mask = fftshift( fft2(pupil_mask) );
        clear x tilt 
    end
    wavefront = prop_multiply(wavefront, pupil_mask);
    clear pupil_mask 
end
if ( use_errors )
    if is_spc && use_pupil_mask
	wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_PUPILMASK_phase_error_V1.0.fits'], 'wavefront');
    else
	wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_PUPILFOLD_phase_error_V1.0.fits'], 'wavefront');
    end
end

wavefront = prop_propagate(wavefront, d_pupilmask_oap5, 'surface_name','OAP5');
wavefront = prop_lens(wavefront, fl_oap5);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP5_phase_error_V3.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap5/2); end

wavefront = prop_propagate(wavefront, d_oap5_fpm+fpm_z_shift_m,'surface_name', 'FPM', 'to_plane');
if ( use_fpm )
    if ( fpm_x_offset ~=0 || fpm_y_offset ~=0 || fpm_x_offset_m ~=0 || fpm_y_offset_m ~=0 )
        %To shift FPM, FFT field to pupil, apply tilt, FFT back to focus,
        %apply FPM, FFT to pupil, take out tilt, FFT back to focus
        x  = repmat( ((1:n)-n/2-1) / (pupil_diam_pix/2),n,1);
        if fpm_x_offset ~=0 || fpm_y_offset ~=0         % shifts are specified in lambda0/D
            tilt =  1i*pi*(x*fpm_x_offset +x'*fpm_y_offset)* lambda0_m / lambda_m;
        elseif fpm_x_offset_m ~=0 || fpm_y_offset_m ~=0 
            tilt =  1i*pi*(x*fpm_x_offset_m+x'*fpm_y_offset_m) / prop_get_sampling(wavefront) * pupil_diam_pix/n ;
        end
        wavefront.wf =  fft2(ifft2(wavefront.wf).*ifftshift(exp(tilt)) );
        clear x 
    end
    if ( is_hlc )
	wavefront0 = ifftshift( ifft2( wavefront.wf ) );	% FFT to virtual pupil
	wavefront0 = wavefront0 * fpm_array(1,1);	% apply uniform amplitude & phase from FPM clear area (assumes clear area is uniform) 
        nfpm = size(fpm_array,1);
	fpm_sampling_lamdivD = fpm_sampling_lam0divD * fpm_lam0_m / lambda_m;    	% FPM sampling at current wavelength in lambda_m/D
        wavefront_fpm = mft2(wavefront0, fpm_sampling_lamdivD, pupil_diam_pix, nfpm, -1); 	% MFT to highly-sampled focal plane
        wavefront_fpm = wavefront_fpm .* fpm_mask .* (fpm_array - 1);		% subtract field inside FPM region, add in FPM-multiplied region
        wavefront_fpm = mft2(wavefront_fpm, fpm_sampling_lamdivD, pupil_diam_pix, n, +1);  	% MFT to virtual pupil
	wavefront0 = wavefront0 + wavefront_fpm;
        wavefront.wf = fft2( fftshift(wavefront0) );	% back to normally-sampled focal plane
        clear wavefront_fpm
        clear wavefront0
    elseif ( is_spc )
        % super-sample FPM
        wavefront0 = custom_pad(ifftshift(ifft2(wavefront.wf)), n_mft);  % to virtual pupil
	if ( isscalar(fpm_array) )
        	fpm_array = fitsread( fpm_file);
		fpm_sampling_lamdivD = fpm_sampling_lam0divD * fpm_lam0_m / lambda_m;  % FPM sampling at current wavelength in lam/D
	else
		m_per_lamdivD = prop_get_sampling(wavefront) / (pupil_diam_pix / n);
		fpm_sampling_lamdivD = fpm_array_sampling_m / m_per_lamdivD;
	end
        nfpm = size(fpm_array,1);
        wavefront0 = mft2(wavefront0, fpm_sampling_lamdivD, pupil_diam_pix, nfpm, -1); % MFT to highly-sampled focal plane
        wavefront0 = wavefront0.*fpm_array;
        wavefront0 = mft2(wavefront0, fpm_sampling_lamdivD, pupil_diam_pix, n, +1);  % MFT to virtual pupil
        wavefront.wf = fft2(fftshift(wavefront0));% back to normally-sampled focal plane
        clear wavefront0
    end
    if ( fpm_x_offset ~=0 || fpm_y_offset ~=0 || fpm_x_offset_m ~= 0 || fpm_y_offset_m ~= 0)
        wavefront.wf = fft2( ifft2(wavefront.wf).*ifftshift(exp(-tilt)) );
        clear tilt
    end
end
if ( pinhole_diam_m ~=0 )
    % "pinhole_diam_m" is pinhole diameter in meters
    dx_m = prop_get_sampling(wavefront);
    dx_pinhole_diam_m = pinhole_diam_m / 101.0;		% 101 samples across pinhole
    n_out = 105;
    m_per_lamD = dx_m * n / pupil_diam_pix;         % current focal plane sampling in lambda_m/D
    dx_pinhole_lamD = dx_pinhole_diam_m /m_per_lamD;% pinhole sampling in lambda_m/D
    n_in = round(pupil_diam_pix * 1.2);
    wavefront0 = ifftshift(ifft2( wavefront.wf ));   % to virtual pupil
    wavefront0 = custom_pad( wavefront0, n_in );
    wavefront0 = mft2( wavefront0, dx_pinhole_lamD, pupil_diam_pix, n_out, -1 );		% MFT to highly-sampled focal plane
    p = (radpix(n_out,n_out) .* dx_pinhole_diam_m) <= pinhole_diam_m/2.0;
    wavefront0 = wavefront0 .* p;
    wavefront0 = mft2( wavefront0, dx_pinhole_lamD, pupil_diam_pix, n, +1 );	% MFT back to virtual pupil
    wavefront.wf = fft2(fftshift(wavefront0)); % back to normally-sampled focal plane
    clear wavefront0 
end

wavefront = prop_propagate(wavefront, d_fpm_oap6-fpm_z_shift_m, 'surface_name','OAP6');
wavefront = prop_lens(wavefront, fl_oap6);
if ( use_errors && ~end_at_fpm_exit_pupil );   wavefront =prop_errormap(wavefront, [map_dir 'roman_phasec_OAP6_phase_error_V3.0.fits'], 'wavefront'); end
if ( use_aperture || pinhole_diam_m ~= 0 );   wavefront = prop_circular_aperture(wavefront, diam_oap6/2); end

if ( end_at_fpm_exit_pupil )
    wavefront = prop_propagate(wavefront, d_oap6_exitpupil, 'surface_name','FPM EXIT PUPIL');
    [wavefront, sampling_m] = prop_end( wavefront, 'NOABS' );
    wavefront = custom_pad( wavefront, n );
    return
end

wavefront = prop_propagate(wavefront, d_oap6_lyotstop, 'surface_name','LYOT_STOP');
if ( use_lyot_stop )
    if ( isscalar(lyot_stop_array) )
    	lyot = fitsread( lyot_stop_file );
    	lyot = custom_pad( lyot, n );    
    else
	lyot = custom_pad( lyot_stop_array, n );
    end
    if ( lyot_x_shift_pupdiam ~=0 || lyot_y_shift_pupdiam ~=0 || lyot_x_shift_m ~=0 || lyot_y_shift_m ~=0 )
        % apply shift to lyot stop by FFTing the stop, applying a tilt, and FFTing back
        x  = repmat( ((1:n)-n/2-1)/(pupil_diam_pix/2),n,1);
        if lyot_x_shift_pupdiam ~=0 || lyot_y_shift_pupdiam ~=0 %offsets are normalized to pupil diameter
            tilt = -1i*pi*(x*lyot_x_shift_pupdiam +x'*lyot_y_shift_pupdiam)* pupil_diam_pix * pupil_diam_pix/n;
        elseif lyot_x_shift_m ~=0 || lyot_y_shift_m ~=0
            tilt = -1i*pi*(x*lyot_x_shift_m +x'*lyot_y_shift_m)/prop_get_sampling(wavefront) * pupil_diam_pix/n;
        end
        lyot = ifft2( ifftshift(lyot) ) .* ifftshift( exp(tilt) );
        lyot = fftshift( fft2(lyot) );
        clear x tilt
    end
    wavefront = prop_multiply( wavefront, lyot );
    clear lyot
end
if ( use_pupil_lens  || pinhole_diam_m ~=0 ); wavefront =prop_circular_aperture(wavefront, 1.079, 'norm'); end

wavefront = prop_propagate(wavefront, d_lyotstop_oap7, 'surface_name','OAP7');
wavefront = prop_lens(wavefront, fl_oap7);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP7_phase_error_V3.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap7/2); end

wavefront = prop_propagate(wavefront, d_oap7_fieldstop, 'surface_name','FIELD_STOP');
if ( use_field_stop )
	if ( isscalar(field_stop_array) && is_hlc )
    		sampling_lamD = pupil_diam_pix / n;      % sampling at focus in lambda_m/D
    		stop_radius = field_stop_radius_lam0 / sampling_lamD * (lambda0_m/lambda_m) * prop_get_sampling(wavefront);
    		if ( field_stop_x_offset ~= 0 || field_stop_y_offset ~= 0 )
          		% convert offsets in lambda0/D to meters
          		x_offset_lamD = field_stop_x_offset * lambda0_m / lambda_m;
          		y_offset_lamD = field_stop_y_offset * lambda0_m / lambda_m;
          		pupil_ratio = pupil_diam_pix / n;
          		field_stop_x_offset_m = x_offset_lamD / pupil_ratio * prop_get_sampling(wavefront);
          		field_stop_y_offset_m = y_offset_lamD / pupil_ratio * prop_get_sampling(wavefront);
    		end
    		wavefront = prop_circular_aperture(wavefront, stop_radius, 'XC', -field_stop_x_offset_m, 'YC', -field_stop_y_offset_m);
    	elseif ( isscalar(field_stop_array) == 0 )
        	wavefront0 = custom_pad(ifftshift(ifft2(wavefront.wf)), n_mft);  % to virtual pupil
		m_per_lamdivD = prop_get_sampling(wavefront) / (pupil_diam_pix / n);
		field_stop_sampling_lamdivD = field_stop_array_sampling_m / m_per_lamdivD;
        	nstop = size(field_stop_array,1);
        	wavefront0 = mft2(wavefront0, field_stop_sampling_lamdivD, pupil_diam_pix, nstop, -1); % MFT to highly-sampled focal plane
        	wavefront0 = wavefront0.*field_stop_array;
        	wavefront0 = mft2(wavefront0, field_stop_sampling_lamdivD, pupil_diam_pix, n, +1);  % MFT to virtual pupil
        	wavefront.wf = fft2(fftshift(wavefront0));% back to normally-sampled focal plane
        	clear wavefront0
	end 
end

wavefront = prop_propagate(wavefront, d_fieldstop_oap8, 'surface_name','OAP8');
wavefront = prop_lens(wavefront, fl_oap8);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_OAP8_phase_error_V3.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap8/2);end

wavefront = prop_propagate(wavefront, d_oap8_filter, 'surface_name','filter');
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'roman_phasec_FILTER_phase_error_V1.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_filter/2); end

% propagate to lens, apply lens, then propagate to detector

if ( use_pupil_lens == 0 && use_defocus_lens == 0 )
        % use imaging lens to create normal focus
        if ( use_errors )
		imaging_lens_error_file = [map_dir 'roman_phasec_LENS_phase_error_V1.0.fits']; 
	else
		imaging_lens_error_file = ' ';
	end
        wavefront = to_from_doublet( wavefront, d_filter_lens, d_lens_fold4+d_fold4_image, ...
                0.10792660718579995, -0.10792660718579995, 0.003, glass_index('S-BSL7R', lambda_m, data_dir), 0.0005, ...
                1e10, 0.10608379812011390, 0.0025, glass_index('PBM2R', lambda_m, data_dir), ...
                'IMAGING LENS', 'IMAGE', imaging_lens_error_file, 1 );
elseif ( use_pupil_lens ) 
        % use pupil imaging lens
        if ( use_errors )
		pupil_lens_error_file = [map_dir 'roman_phasec_PUPILLENS_phase_error_V1.0.fits']; 
	else
		pupil_lens_error_file = ' ';
	end
        wavefront = to_from_doublet( wavefront, d_filter_lens, d_lens_fold4+d_fold4_image, ...
                0.03449714, -0.20846556, 0.003, glass_index('S-BSL7R', lambda_m, data_dir), 0.0007, ...
                1e10, 0.05453792, 0.0025, glass_index('PBM2R', lambda_m, data_dir), ...
                'PUPIL IMAGING LENS', 'IMAGE', pupil_lens_error_file, 0 );
elseif ( use_defocus_lens ~= 0 )     
        % use one of 4 defocusing lenses
        c_lens = [0.01325908247149297, 0.01243982235933671, 0.01163480668768688, 0.005663476241717166];
        if ( use_errors )
		defocus_lens_error_file = [map_dir 'roman_phasec_DEFOCUSLENS' num2str(round(use_defocus_lens)) '_phase_error_V1.0.fits']; 
	else
		defocus_lens_error_file = ' ';
	end
        wavefront = to_from_singlet( wavefront, d_filter_lens, d_lens_fold4+d_fold4_image, ...
                 0.001 / c_lens(use_defocus_lens), 1e10, 0.005, glass_index('SILICA', lambda_m, data_dir), ...
                'DEFOCUS LENS', 'IMAGE', defocus_lens_error_file, use_defocus_lens == 3 || use_defocus_lens == 4 );
end

if ( end_at_exit_pupil )
   % create a sharp pupil image
   wavefront = prop_propagate( wavefront, 0.1 );
   wavefront = prop_lens( wavefront, 0.1 );
   wavefront = prop_propagate( wavefront, 0.0882 );
elseif ( use_defocus_lens == 0 && use_pupil_lens == 0 ) 
   rsqr = prop_radius( wavefront ).^2;
end

[wavefront, sampling_m] = prop_end(wavefront, 'NOABS');

% remove phase term due to pupil not at front focus of lens (only in direct imaging mode)

if use_defocus_lens == 0 && use_pupil_lens == 0 && end_at_exit_pupil == 0
        wavefront = wavefront .* exp((complex(0,1) * 16 * pi / lambda_m * 0.0735) * rsqr);
end

% resample as needed

if ( final_sampling_m ~=0 )
        mag = sampling_m / (final_sampling_m);
        sampling_m = final_sampling_m;
    	wavefront = prop_magnify( wavefront, mag, 'size_out',output_dim,'amp_conserve');
elseif ( final_sampling_lam0 ~=0 )
	if ( use_pupil_lens ~= 0 || use_defocus_lens ~= 0 || end_at_exit_pupil ~= 0 )
    		disp('ERROR: Cannot specify final_sampling_lam0 when using pupil or defocus lens or ending at exit pupil.')
    		return
	else
        	mag = (pupil_diam_pix/n) / final_sampling_lam0 * (lambda_m/lambda0_m);
        	sampling_m = sampling_m / mag;
    		wavefront = prop_magnify( wavefront, mag, 'size_out',output_dim,'amp_conserve');
    	end
else
	wavefront = custom_pad(wavefront, output_dim);
end

return
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function n = glass_index( glass, lambda_m, data_dir )

file = fopen( [data_dir '/glass/' glass '_index.txt'], 'r' );
val = fscanf( file, '%f %f', [2,101] );
fclose( file );

lam_um = val(1,:);
index = val(2,:);

n = spline( lam_um, index, lambda_m*1e6 );

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function wavefront = to_from_singlet( wavefront, dz_to_lens, dz_from_lens, r1, r2, thickness, nglass, ...
			surface_name, next_surface_name, error_file, to_plane )
% single thick lens
% r1, r2 = front & rear radii of curvature

f = 1 / ( (nglass - 1) * ( 1.0/r1 - 1.0/r2 + (nglass - 1)*thickness / (nglass*r1*r2) ) );      % effective focal length
h1 = -f*(nglass - 1)*thickness / (nglass*r2);    % front principal plane
h2 = -f*(nglass - 1)*thickness / (nglass*r1);    % rear principal plane

wavefront = prop_propagate( wavefront, dz_to_lens+h1, 'SURFACE_NAME', surface_name );
wavefront = prop_lens( wavefront, f, surface_name );
if error_file ~= ' '; wavefront = prop_errormap( wavefront, error_file, 'wavefront' ); end
if ( to_plane )
   wavefront = prop_propagate( wavefront, -h2+dz_from_lens, 'SURFACE_NAME', next_surface_name, 'TO_PLANE' );
else
   wavefront = prop_propagate( wavefront, -h2+dz_from_lens, 'SURFACE_NAME', next_surface_name );
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function wavefront = to_from_doublet( wavefront, dz_to_lens, dz_from_lens, r1_a, r2_a, thickness_a, nglass_a, separation_m, ... 
                        r1_b, r2_b, thickness_b, nglass_b, surface_name, next_surface_name, error_file, to_plane )

% two air-gapped thick lenses (lenses a & b)

f_a = 1 / ( (nglass_a - 1) * ( 1.0/r1_a - 1.0/r2_a + (nglass_a - 1)*thickness_a / (nglass_a*r1_a*r2_a) ) );
h1_a = -f_a*(nglass_a - 1)*thickness_a / (nglass_a*r2_a);
h2_a = -f_a*(nglass_a - 1)*thickness_a / (nglass_a*r1_a);

f_b = 1 / ( (nglass_b - 1) * ( 1.0/r1_b - 1.0/r2_b + (nglass_b - 1)*thickness_b / (nglass_b*r1_b*r2_b) ) );
h1_b = -f_b*(nglass_b - 1)*thickness_b / (nglass_b*r2_b);
h2_b = -f_b*(nglass_b - 1)*thickness_b / (nglass_b*r1_b);

wavefront = prop_propagate( wavefront, dz_to_lens+h1_a, 'SURFACE_NAME', surface_name );
wavefront = prop_lens( wavefront, f_a, [surface_name ' lens #1'] );
if error_file ~= ' '; wavefront = prop_errormap( wavefront, error_file, 'wavefront' ); end
wavefront = prop_propagate( wavefront, -h2_a+separation_m+h1_b );
wavefront = prop_lens( wavefront, f_b, [surface_name ' lens #2'] );
if ( to_plane )
   wavefront = prop_propagate( wavefront, -h2_b+dz_from_lens, 'SURFACE_NAME', next_surface_name, 'TO_PLANE' );
else
   wavefront = prop_propagate( wavefront, -h2_b+dz_from_lens, 'SURFACE_NAME', next_surface_name );
end

end
