%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

%function [ wavefront, sampling_m ] =  wfirst_phaseb(lambda_m, output_dim0, optval)

% Version 1.0, 3 January 2019, JEK
% Version 1.0a, 9 January 2019, JEK
%    Changes: 1) now append polaxis to output_field_rootname
%             2) can output field at FPM exit pupil (via output_field_rootname) without
%                having to stop there
% Version 1.0b, 16 January 2019, JEK
%    Changes: Multiplied wfirst_phaseb_GROUND_TO_ORBIT_phase_error_V1.0.fits by 4.2x to get the
%		  the RMS WFE at the CGI entrance to be 76.4 nm RMS.  The new maps is in
%		  result to wfirst_phaseb_GROUND_TO_ORBIT_4.2X_phase_error_V1.0.fits
% Version 1.0c, 13 February 2019, JEK
%    Changes: Changed HLC design to HLC-20190210; added 'fpm_axis' parameter to specify
%		  axis of HLC FPM ('s' || 'p') due to FPM tilt; changed 'spc' coronagraph
%		  option to 'cor_type' to 'spc-ifs' and changed default SPC to SPC-20190130 and
%		  updated associated parameters; added 'spc-wide' option to 'cor_type' to use
%		  SPC-20181220. Changed DM tilt to rotation around Y axis.
% Version 1.2, 13 March 2019, JEK
%    Changes: Added new parameters: cgi & mask & lyot stop & FPM shifts in meters; final
%		  sampling in meters; simplified how some masks are shifted; changed option
%		  of cor_type='spc-ifs' to 'spc-ifs_short' for 660 nm FPM, 'spc-ifs_long' for
%		  730 nm FPM.  Rotated SPC masks to match orientation of pupil in HLC.
%		  Added Erkin's HLC design.  dm_sampling set to 0.9906 mm.
% Version 1.4, 25 Sept 2019, JEK
%    Changes: Added optional use_pupil_mask parameter; fixed error that omitted pupil mask
%    if end_at_fpm_exit_pupil was set (Matlab only) 
% Version 1.5, 7 Oct 2019, JKE
%    Changes: Aliased spc-spec_* to spc-ifs_*
% Version 1.6, 18 Oct 2019, JEK
%    Changes: Added option to shift field stop (HLC)
% Version 1.7, 3 Dec 2019, JEK
%    Changes: Replaced defocus-vs-focal-length table with revised version; added logic to
%             force planar propagation when requested defocus is <=4 waves.
%
% Based on Zemax prescription "WFIRST_CGI_DI_LOWFS_Sep24_2018.zmx" by Hong Tang.
% IDL version (original), John E Krist
% Matlab translation, 3/1/2019, Hanying Zhou
% Last modified, 3/14/2019, Hanying Zhou
% Modified on 2019-04-04 by A.J. Riggs to clean up the syntax and to make
% the function work better with an external WFSC wrapper.
% Modified on 2019-05-01 by A.J. Riggs to include bug fixes in some of the 
% if statements from Hanying Zhou.
% Modified on 2019-05-06 by A.J. Riggs to remove the "/" or "\" at the end
% of the data_dir path if the slash exists.
% Modified on 2019-05-07 by A.J. Riggs to zero out the phase in the outer
% part of the HLC occulters get rid of arbitrary phase differences between
% at different wavelengths.
% Modified on 2019-05-07 by A.J. Riggs to include the case 'spc_ifs_custom'
% for custom SPC-IFS designs.
% Modified on 2019-06-06 by J. Krist to take out the code that subtracts the
% phase offset; the FPM files have instead been phase adjusted

function [wavefront,  sampling_m]= model_full_wfirst_phaseb(lambda_m, output_dim0, optval)

% "output_dim" is used to specify the output dimension in pixels at the final image plane.
% The computational grid sizes are hardcoded for each coronagraph.

%--mask design data path

data_dir = '/home/krist/cgi/phaseb/phaseb_data';  % no trailing '/'

if ~exist('optval'); optval.dummy = 0; end;

if ( isfield(optval,'data_dir') )
    data_dir = optval.data_dir;
    if(strcmp(data_dir(end),'/') || strcmp(data_dir(end),'\'))
        data_dir = data_dir(1:end-1); %--Remove the trailing slash for this function.
    end
end

map_dir = [data_dir '/maps/'];            % directory for surface maps
polfile = [data_dir '/pol/new_toma_'];		% polarization aberration table rootname

cor_type = 'hlc';           % 'hlc', 'spc-ifs', 'spc-wide', || 'none' (none = clear aperture, no coronagraph)
source_x_offset_mas = 0;	% source offset in milliarcsec (tilt applied at primary)
source_y_offset_mas = 0;
source_x_offset = 0;		% source offset in lambda0_m/D radians (tilt applied at primary)
source_y_offset = 0;
polaxis = 0;                % polarization axis aberrations:
%     -2 = -45; in, Y out
%     -1 = -45; in, X out
%      1 = +45; in, X out
%      2 = +45; in, Y out
%      5 = mean of modes -1 & +1 (X channel polarizer)
%      6 = mean of modes -2 & +2 (Y channel polarizer)
%     10 = mean of all modes (no polarization filtering)
use_errors = 1;             % use optical surface phase errors? 1 or 0
zindex = 0;                 % array of Zernike polynomial indices
zval_m = 0;                 % array of Zernike coefficients (meters RMS WFE)
use_aperture = 0;           % use apertures on all optics? 1 or 0
cgi_x_shift_pupdiam = 0;	% X,Y shift of wavefront at FSM (bulk displacement of CGI); normalized relative to pupil diameter
cgi_y_shift_pupdiam = 0;
cgi_x_shift_m = 0;			%X,Y shift of wavefront at FSM (bulk displacement of CGI) in meters
cgi_y_shift_m = 0;
end_at_fsm = 0;             % end propagation after propagating to FSM (no FSM errors)
fsm_x_offset_mas = 0;		% offset source using FSM (milliarcsec)
fsm_y_offset_mas = 0;
fsm_x_offset = 0;           % offset source using FSM (lambda0/D radians)
fsm_y_offset = 0;
focm_z_shift_m = 0;         % offset (meters) of focus correction mirror (+ increases path length)
use_hlc_dm_patterns = 0;	% use Dwight-generated HLC default DM wavefront patterns? 1 or 0
use_dm1 = 0;                % use DM1? 1 or 0
use_dm2 = 0;                % use DM2? 1 or 0
dm_sampling_m = 0.9906e-3;  % actuator spacing in meters;
dm1_xc_act = 23.5;          % for 48x48 DM(wavefront centered at actuator intersections: (0,0) = 1st actuator center
dm1_yc_act = 23.5;
dm1_xtilt_deg = 0;   		% tilt around X axis
dm1_ytilt_deg = 5.7;		% effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
dm1_ztilt_deg = 0;
dm2_xc_act = 23.5;
dm2_yc_act = 23.5;
dm2_xtilt_deg = 0;
dm2_ytilt_deg = 5.7;
dm2_ztilt_deg = 0;
use_pupil_mask = 1;		% (SPC only) Use SPC pupil mask (0 or 1)
mask_x_shift_pupdiam = 0;	% X,Y shift of shaped pupil mask; normalized relative to pupil diameter
mask_y_shift_pupdiam = 0;
mask_x_shift_m = 0;         % X,Y shift of shaped pupil mask; 
mask_y_shift_m = 0;
use_fpm = 1;                % use occulter? 1 or 0
fpm_axis = 'p';             % HLC FPM axis: '', 's', 'p'
fpm_x_offset = 0;           % FPM x,y offset in lambda0/D
fpm_y_offset = 0;
fpm_x_offset_m = 0;          % FPM x,y offset in meters
fpm_y_offset_m = 0;
fpm_z_shift_m = 0;          % FPM offset in meters along optical axis (+ = away from prior optics)
pinhole_diam_m = 0;         % pinhole diameter in meters at FPM
end_at_fpm_exit_pupil = 0;	% return field at FPM exit pupil?
getWFE = 0;                 % return wfe instead of field
output_field_rootname = '';	% rootname of FPM exit pupil field file (must set end_at_fpm_exit_pupil=1)
use_lyot_stop = 1;          % use Lyot stop? 1 or 0
lyot_x_shift_pupdiam = 0;	% X,Y shift of Lyot stop mask; normalized relative to pupil diameter
lyot_y_shift_pupdiam = 0;
lyot_x_shift_m = 0;			%X,Y shift of Lyot stop mask in meters
lyot_y_shift_m = 0;
use_field_stop = 1;         % use field stop (HLC)? 1 or 0
field_stop_radius_lam0 = 0;	% field stop radius in lambda0/D (HLC || SPC-wide mask only)
field_stop_x_offset = 0;        % field stop offset in lambda0/D
field_stop_y_offset = 0;
field_stop_x_offset_m = 0;      % field stop offset in meters
field_stop_y_offset_m = 0;
use_pupil_lens = 0;         % use pupil imaging lens?
use_defocus_lens = 0;		% use defocusing lens? Options are 1, 2, 3, 4, corresponding to +18.0, +9.0, -4.0, -8.0 waves P-V @ 550 nm
defocus = 0;                % instead of specific lens, defocus in waves P-V @ 550 nm (-8.7 to 42.0 waves)
final_sampling_lam0 = 0;	% final sampling in lambda0/D
final_sampling_m = 0;       % final sampling in meters (overrides final_sampling_lam0)
output_dim = output_dim0;	% dimension of output grid in pixels (overrides output_dim0)

if exist('optval','var')==1 %if  exist('optval')
    if ( isfield(optval,'cor_type') );   cor_type = optval.cor_type; end
    if ( isfield(optval,'use_fpm') );    use_fpm = optval.use_fpm; end
    if ( isfield(optval,'fpm_axis') );   fpm_axis = optval.fpm_axis; end
end

is_hlc = 0;
is_spc = 0;
is_fpm_grid = 0;

if  strcmp(cor_type,'hlc')
    is_hlc = 1;
    file_directory = [data_dir '/hlc_20190210/'];         % must have trailing "/"
    prefix = [file_directory  'run461_'];
    pupil_diam_pix = 309.0;
    pupil_file = [prefix  'pupil_rotated.fits'];
    lyot_stop_file = [prefix  'lyot.fits'];
    lambda0_m = 0.575e-6;
    nlams = 19 ;              % number of occ trans lambda provided
    bw = 0.1;
    lam_occ = linspace(1-bw/2, 1+bw/2, nlams)*lambda0_m; %[(1-bw/2):bw/(nlams-mod(nlams,2)):(1+bw/2)]*lambda0_m; 	% wavelengths at which occ trans provided
    wlam = find( round(1e13*lambda_m) == round(1e13*lam_occ) ); 	% find exactly matching FPM wavelength
    occulter_file_r = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'real.fits'];
    occulter_file_i = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'imag.fits'];
    n_default = 1024;	% gridsize in non-critical areas
    if  use_fpm;    n_to_fpm = 2048; else; n_to_fpm = 1024; end
    n_from_lyotstop = 1024;
    field_stop_radius_lam0 = 9.0;
elseif  strcmp(cor_type,'hlc_erkin')
    is_hlc = 1;
    file_directory = [data_dir '/hlc_20190206_v3/'];         % must have trailing "/"
    prefix = [file_directory  'dsn17d_run2_pup310_fpm2048_'];
    pupil_diam_pix = 310.0;
    pupil_file = [prefix  'pupil.fits'];
    lyot_stop_file = [prefix  'lyot.fits'];
    fpm_axis = 's';
    lambda0_m = 0.575e-6;
    nlams = 19 ;              % number of occ trans lambda provided
    bw = 0.1;
    lam_occ = [(1-bw/2):bw/(nlams-mod(nlams,2)):(1+bw/2)]*lambda0_m; % wavelengths at which occ trans provided
    [~,wlam] = min(abs(lambda_m-lam_occ)); % find nearest matching FPM wavelength
    occulter_file_r = [prefix  'occ_lam' num2str(lam_occ(wlam),5) 'theta6.69pol'   fpm_axis  '_real_rotated.fits'];
    fprintf('occulter_file_r = %s\n',occulter_file_r);
    occulter_file_i = [prefix  'occ_lam' num2str(lam_occ(wlam),5) 'theta6.69pol'   fpm_axis  '_imag_rotated.fits'];
    n_default = 1024;	% gridsize in non-critical areas
    if  use_fpm;    n_to_fpm = 2048; else; n_to_fpm = 1024; end
    n_from_lyotstop = 1024;
    field_stop_radius_lam0 = 9.0;
elseif  strcmp(cor_type,'hlc_custom')  
    is_hlc = 1;  
    if(isfield(optval,'hlc_name')==false || isfield(optval,'prefix')==false)
        error('You must define the variables hlc_name and prefix when using hlc_custom as the coronagraph.')
    end
    file_directory = [data_dir filesep 'hlc_custom' filesep optval.hlc_name filesep];         % must have trailing "/"
    prefix = [file_directory  optval.prefix];
    pupil_diam_pix = 309.0;
    pupil_file = [prefix  'pupil_rotated.fits'];
    lyot_stop_file = [prefix  'lyot.fits'];
    %--For CUSTOM HLC ONLY: Defined again because model_full_wfirst_phaseb.m has too many hard-coded values
    lambda0_m = optval.lambda0_m;
    nlams = optval.nlams;             % number of occ trans lambda provided
    bw = optval.bw;
    thetaString = optval.thetaString;
    if(nlams==1)
        lam_occ = lambda0_m;
    else
        lam_occ = linspace(1-bw/2,1+bw/2,nlams)*lambda0_m;%[(1-bw/2):bw/(nlams-mod(nlams,2)):(1+bw/2)]*lambda0_m; 	% wavelengths at which occ trans provided
    end
    wlam = find( round(1e14*lambda_m) == round(1e14*lam_occ) ); 	% find exactly matching FPM wavelength
    occulter_file_r = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta' thetaString 'pol'   fpm_axis   '_' 'real.fits'];
    occulter_file_i = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta' thetaString 'pol'   fpm_axis   '_' 'imag.fits'];

    %--FOR TILED FPM STUDY ONLY
    if(isfield(optval, 'fpm_grid_sep'))
        fpm_grid_sep = optval.fpm_grid_sep;
    else
        fpm_grid_sep = 30; % lambda0/D
    end
    if(isfield(optval, 'is_fpm_grid'))
        is_fpm_grid = optval.is_fpm_grid;
    else
        is_fpm_grid = false;
    end
    if(isfield(optval, 'is_fpm_grid_hex'))
        is_fpm_grid_hex = optval.is_fpm_grid_hex;
    else
        is_fpm_grid_hex = false;
    end
    if(is_fpm_grid)
        if(is_fpm_grid_hex)
            fpm_grid = fitsread([prefix  'occ_grid_hexSep' num2str(fpm_grid_sep) '_lam' num2str(lambda_m,12) '_' 'abs.fits']);
        else
            fpm_grid = fitsread([prefix  'occ_grid_sep' num2str(fpm_grid_sep) '_lam' num2str(lambda_m,12) '_' 'abs.fits']);
        end
%         figure(221); imagesc(fpm_grid); axis xy equal tight; colorbar; drawnow;
    end
    
    n_default = 1024;	% gridsize in non-critical areas
    if  use_fpm;    n_to_fpm = 2048; else; n_to_fpm = 1024; end
    n_from_lyotstop = 1024;
    field_stop_radius_lam0 = 9.0; 
elseif(strcmpi(cor_type, 'spc_spec_custom') || strcmpi(cor_type, 'spc_ifs_custom')) 
    is_spc = 1;
    pupil_file = optval.pupil_file;
    pupil_diam_pix = optval.pupil_diam_pix; %1000;
    pupil_mask_file = optval.pupil_mask_file;% [file_dir  'SPM_SPC-20190130.fits']; %--SPM file name
    fpm_file = optval.fpm_file; %[file_dir 'fpm_0.05lamdivD.fits'];
    fpm_sampling_lam0 = optval.fpm_sampling_lam0;% 0.05; 	% sampling in lambda0/D of FPM mask
    lyot_stop_file = optval.lyot_stop_file; %[file_dir  'LS_SPC-20190130.fits'];
    lambda0_m = optval.lambda0_m; %0.73e-6;        % FPM scaled for this central wavelength
    
    n_default = 2048;           % gridsize in non-critical areas
    n_to_fpm = 2048;            % gridsize to/from FPM
    n_mft = 1400;
    n_from_lyotstop = 4096;
elseif  sum(strfind(cor_type, 'spc-ifs')) || sum(strfind(cor_type, 'spc-spec'))
    is_spc = 1;
    file_dir = [data_dir '/spc_20190130/'];       % must have trailing "/"
    pupil_diam_pix = 1000;
    pupil_file = [file_dir  'pupil_SPC-20190130_rotated.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20190130.fits'];
    fpm_file = [file_dir 'fpm_0.05lamdivD.fits'];
    fpm_sampling_lam0 = 0.05; 	% sampling in lambda0/D of FPM mask
    lyot_stop_file = [file_dir  'LS_SPC-20190130.fits'];
    lambda0_m = 0.73e-6;        % FPM scaled for this central wavelength
    if strcmp(cor_type, 'spc-spec_short'); lambda0_m = 0.66e-6;  end      
    n_default = 2048;           % gridsize in non-critical areas
    n_to_fpm = 2048;            % gridsize to/from FPM
    n_mft = 1400;
    n_from_lyotstop = 4096;
elseif  strcmp(cor_type, 'spc-wide' )
    is_spc = 1;
    file_dir = [data_dir '/spc_20181220/'];       % must have trailing "/"
    pupil_diam_pix = 1000;
    pupil_file = [file_dir  'pupil_SPC-20181220_1k_rotated.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20181220_1000_rounded9_gray.fits'];
    fpm_file = [file_dir   'fpm_0.05lamdivD.fits'];
    fpm_sampling_lam0 = 0.05;       % sampling in lambda0/D of FPM mask
    lyot_stop_file = [file_dir  'LS_SPC-20181220_1k.fits'];
    lambda0_m = 0.825e-6;       % FPM scaled for this central wavelength
    n_default = 2048;           % gridsize in non-critical areas
    n_to_fpm = 2048;            % gridsize to/from FPM
    n_mft = 1400;
    n_from_lyotstop = 4096;
elseif strcmp(cor_type, 'none' )
    file_directory = './hlc_20190210/';         % must have trailing "/"
    prefix = file_directory + 'run461_';
    pupil_file = prefix + 'pupil_rotated.fits';
    lambda0_m = 0.575e-6;
    use_fpm = 0;
    use_lyot_stop = 0;
    use_field_stop =0;
    pupil_diam_pix = 309.0;
    n_default = 1024;
    n_to_fpm = 1024;
    n_from_lyotstop = 1024;
else
    disp([ 'cor_type = ' cor_type ' is not supported in model_full_wfirst_phaseb'])
    return
end

mas_per_lamD = lambda0_m * 360.0 * 3600.0 / (2 * pi * 2.363) * 1000;	% mas per lambda0/D

if(exist('optval','var')==1) %if exist('optval')
    if ( isfield(optval,'lam0') );           lambda0_m = optval.lam0 * 1.0e-6;end
    if ( isfield(optval,'lambda0_m') );      lambda0_m = optval.lambda0_m;end
    if ( isfield(optval,'source_x_offset') );   source_x_offset = optval.source_x_offset; end
    if ( isfield(optval,'source_y_offset') );   source_y_offset = optval.source_y_offset;end
    if ( isfield(optval,'source_x_offset_mas') );   source_x_offset = optval.source_x_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'source_y_offset_mas') );   source_y_offset = optval.source_y_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'use_errors') );     use_errors = optval.use_errors;end
    if ( isfield(optval,'polaxis') );        polaxis = optval.polaxis; end
    if ( isfield(optval,'zindex') );         zindex = optval.zindex;end
    if ( isfield(optval,'zval_m') );         zval_m = optval.zval_m; end
    if ( isfield(optval,'end_at_fsm') );     end_at_fsm = optval.end_at_fsm; end
    if ( isfield(optval,'cgi_x_shift_pupdiam') );   cgi_x_shift_pupdiam = optval.cgi_x_shift_pupdiam; end
    if ( isfield(optval,'cgi_y_shift_pupdiam') );   cgi_y_shift_pupdiam = optval.cgi_y_shift_pupdiam; end
    if ( isfield(optval,'cgi_x_shift_m') );  cgi_x_shift_m = optval.cgi_x_shift_m; end
    if ( isfield(optval,'cgi_y_shift_m') );  cgi_y_shift_m = optval.cgi_y_shift_m; end
    if ( isfield(optval,'fsm_x_offset') );   fsm_x_offset = optval.fsm_x_offset; end
    if ( isfield(optval,'fsm_y_offset') );   fsm_y_offset = optval.fsm_y_offset; end
    if ( isfield(optval,'fsm_x_offset_mas') );      fsm_x_offset = optval.fsm_x_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'fsm_y_offset_mas') );      fsm_y_offset = optval.fsm_y_offset_mas / mas_per_lamD; end
    if ( isfield(optval,'focm_z_shift_m') );        focm_z_shift_m = optval.focm_z_shift_m; end
    if ( isfield(optval,'use_hlc_dm_patterns') );   use_hlc_dm_patterns = optval.use_hlc_dm_patterns; end
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
    if ( isfield(optval,'mask_x_shift_pupdiam') );   mask_x_shift_pupdiam = optval.mask_x_shift_pupdiam; end
    if ( isfield(optval,'mask_y_shift_pupdiam') );   mask_y_shift_pupdiam = optval.mask_y_shift_pupdiam; end
    if ( isfield(optval,'mask_x_shift_m') );     mask_x_shift_m = optval.mask_x_shift_m; end
    if ( isfield(optval,'mask_y_shift_m') );     mask_y_shift_m = optval.mask_y_shift_m; end
    if ( isfield(optval,'fpm_x_offset') );       fpm_x_offset = optval.fpm_x_offset; end
    if ( isfield(optval,'fpm_y_offset') );       fpm_y_offset = optval.fpm_y_offset; end
    if ( isfield(optval,'fpm_x_offset_m') );      fpm_x_offset_m = optval.fpm_x_offset_m; end
    if ( isfield(optval,'fpm_y_offset_m') );      fpm_y_offset_m = optval.fpm_y_offset_m; end
    if ( isfield(optval,'fpm_z_shift_m') );      fpm_z_shift_m = optval.fpm_z_shift_m; end
    if ( isfield(optval,'pinhole_diam_m') );     pinhole_diam_m = optval.pinhole_diam_m; end
    if ( isfield(optval,'end_at_fpm_exit_pupil') );  end_at_fpm_exit_pupil = optval.end_at_fpm_exit_pupil; end
    if ( isfield(optval,'output_field_rootname') );  output_field_rootname = optval.output_field_rootname; end
    if ( isfield(optval,'use_lyot_stop') );          use_lyot_stop = optval.use_lyot_stop; end
    if ( isfield(optval,'lyot_x_shift_pupdiam') );   lyot_x_shift_pupdiam = optval.lyot_x_shift_pupdiam; end
    if ( isfield(optval,'lyot_y_shift_pupdiam') );   lyot_y_shift_pupdiam = optval.lyot_y_shift_pupdiam; end
    if ( isfield(optval,'lyot_x_shift_m') );         lyot_x_shift_m = optval.lyot_x_shift_m; end
    if ( isfield(optval,'lyot_y_shift_m') );         lyot_y_shift_m = optval.lyot_y_shift_m; end
    if ( isfield(optval,'use_field_stop') );         use_field_stop = optval.use_field_stop; end
    if ( isfield(optval,'field_stop_radius_lam0') ); field_stop_radius_lam0 = optval.field_stop_radius_lam0; end
    if ( isfield(optval,'field_stop_x_offset') );  field_stop_x_offset = optval.field_stop_x_offset; end
    if ( isfield(optval,'field_stop_y_offset') );  field_stop_y_offset = optval.field_stop_y_offset; end
    if ( isfield(optval,'field_stop_x_offset_m') );  field_stop_x_offset_m = optval.field_stop_x_offset_m; end
    if ( isfield(optval,'field_stop_y_offset_m') );  field_stop_y_offset_m = optval.field_stop_y_offset_m; end
    if ( isfield(optval,'use_pupil_lens') );         use_pupil_lens = optval.use_pupil_lens; end
    if ( isfield(optval,'use_defocus_lens') );       use_defocus_lens = optval.use_defocus_lens; end
    if ( isfield(optval,'defocus') );                defocus = optval.defocus; end
    if ( isfield(optval,'final_sampling_lam0') );    final_sampling_lam0 = optval.final_sampling_lam0; end
    if ( isfield(optval,'final_sampling_m') );       final_sampling_m = optval.final_sampling_m; end
    if ( isfield(optval,'output_dim') );             output_dim = optval.output_dim; end
end


diam = 2.3633372;
fl_pri = 2.83459423440 * 1.0013;
d_pri_sec = 2.285150515460035;
d_focus_sec = d_pri_sec - fl_pri;
fl_sec = -0.653933011 * 1.0004095;
d_sec_focus = 3.580188916677103; 	
diam_sec = 0.58166;
d_sec_fold1 = 2.993753476654728;
d_fold1_focus = 0.586435440022375;	
diam_fold1 = 0.09;
d_fold1_m3 = 1.680935841598811;
fl_m3 = 0.430216463069001;
d_focus_m3 = 1.094500401576436;	
d_m3_pupil = 0.469156807701977;	
d_m3_focus = 0.708841602661368;	
diam_m3 = 0.2;
d_m3_m4 = 0.943514749358944;
fl_m4 = 0.116239114833590;
d_focus_m4 = 0.234673014520402;	
d_pupil_m4 = 0.474357941656967;	
d_m4_focus = 0.230324117970585;	
diam_m4 = 0.07;
d_m4_m5 = 0.429145636743193;
d_m5_focus = 0.198821518772608;	
fl_m5 = 0.198821518772608;		
d_m5_pupil = 0.716529242882632;	
diam_m5 = 0.07;
d_m5_fold2 = 0.351125431220770;
diam_fold2 = 0.06;
d_fold2_fsm = 0.365403811661862;	
d_fsm_oap1 = 0.354826767220001;
fl_oap1 = 0.503331895563883;
diam_oap1 = 0.06;
d_oap1_focm = 0.768005607094041;
d_focm_oap2 = 0.314483210543378;
fl_oap2 = 0.579156922073536;
diam_oap2 = 0.06;
d_oap2_dm1 = 0.775775726154228;		
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
d_pupilmask_oap5 = 0.408810648253099; 	
fl_oap5 =  0.548189351937178;
diam_oap5 = 0.06;
d_oap5_fpm = 0.548189083164429;
d_fpm_oap6 = 0.548189083164429;		
fl_oap6 = 0.548189083164429;		
diam_oap6 = 0.06;
d_oap6_lyotstop = 0.687567667550736;	
d_lyotstop_oap7 = 0.401748843470518;
fl_oap7 = 0.708251083480054;
diam_oap7 = 0.06;
d_oap7_fieldstop = 0.708251083480054;	
d_fieldstop_oap8 = 0.210985967281651;
fl_oap8 = 0.210985967281651;		
diam_oap8 = 0.06;
d_oap8_pupil = 0.238185804200797;	
d_oap8_filter = 0.368452268225530;
diam_filter = 0.01;
d_filter_lens = 0.170799548215162;
fl_lens = 0.246017378417573 + 0.050001306014153;
diam_lens = 0.01;
d_lens_fold4 = 0.246017378417573;
diam_fold4 = 0.02;
d_fold4_image = 0.050001578514650;


fl_pupillens = 0.149260576823040;	

n = n_default;	% start off with less padding

wavefront = prop_begin(diam,lambda_m, n,'beam_diam_fraction', pupil_diam_pix/n);%

pupil =fitsread( pupil_file);
wavefront = prop_multiply(wavefront, custom_pad(pupil,n));
clear pupil

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
if ( use_errors )
    wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_PRIMARY_phase_error_V1.0.fits'], 'wavefront');
    wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_GROUND_TO_ORBIT_4.2X_phase_error_V1.0.fits'], 'wavefront');
end

wavefront = prop_propagate(wavefront, d_pri_sec, 'surface_name','secondary');
wavefront = prop_lens(wavefront, fl_sec, 'secondary');
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_SECONDARY_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_sec/2 );end

wavefront = prop_propagate(wavefront, d_sec_fold1, 'surface_name','FOLD_1');
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_FOLD1_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture( wavefront, diam_fold1/2 );end

wavefront = prop_propagate(wavefront, d_fold1_m3, 'surface_name','M3');
wavefront = prop_lens(wavefront, fl_m3);
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_M3_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_m3/2);end

wavefront = prop_propagate(wavefront, d_m3_m4, 'surface_name','M4');
wavefront = prop_lens(wavefront, fl_m4);
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_M4_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_m4/2);end

wavefront = prop_propagate(wavefront, d_m4_m5, 'surface_name','M5');
wavefront = prop_lens(wavefront, fl_m5);
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_M5_phase_error_V1.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_m5/2 ); end

wavefront =prop_propagate(wavefront, d_m5_fold2, 'surface_name','FOLD_2');
if ( use_errors );   wavefront = prop_errormap(wavefront, [map_dir 'wfirst_phaseb_FOLD2_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_fold2/2); end

wavefront = prop_propagate(wavefront, d_fold2_fsm, 'surface_name','FSM');

if ( end_at_fsm)
    [wavefront, sampling_m] = prop_end(wavefront, 'noabs');
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
    clear x tilt %x =0; tilt =0;
end

if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_FSM_phase_error_V1.0.fits'], 'wavefront'); end
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
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP1_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap1/2);end

wavefront = prop_propagate(wavefront, d_oap1_focm+focm_z_shift_m, 'surface_name','FOCM');
if ( use_errors );  wavefront =  prop_errormap(wavefront,[map_dir 'wfirst_phaseb_FOCM_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );  wavefront =  prop_circular_aperture(wavefront, diam_focm/2);end

wavefront =prop_propagate(wavefront, d_focm_oap2+focm_z_shift_m, 'surface_name','OAP2');
wavefront = prop_lens(wavefront, fl_oap2);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP2_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap2/2 );end

wavefront = prop_propagate(wavefront, d_oap2_dm1, 'surface_name','DM1');
if ( use_dm1 ); wavefront = prop_dm(wavefront, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m, 'xtilt',dm1_xtilt_deg, 'ytilt',dm1_ytilt_deg, 'ztilt',dm1_ztilt_deg); end
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_DM1_phase_error_V1.0.fits'], 'wavefront');end
if ( is_hlc  && use_hlc_dm_patterns )
    dm1wfe = fitsread( [prefix 'dm1wfe.fits']);
    wavefront = prop_add_phase(wavefront, custom_pad(dm1wfe, n));
    clear dm1wfe 
end

wavefront = prop_propagate(wavefront, d_dm1_dm2, 'surface_name','DM2');
if ( use_dm2 );   wavefront = prop_dm(wavefront, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m, 'ztilt',dm2_xtilt_deg, 'ytilt',dm2_ytilt_deg, 'ztilt',dm2_ztilt_deg);end
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_DM2_phase_error_V1.0.fits'], 'wavefront');end
if ( is_hlc ) 
    if ( use_hlc_dm_patterns )
        dm2wfe =fitsread( [prefix 'dm2wfe.fits']);
        wavefront = prop_add_phase(wavefront, custom_pad(dm2wfe, n));
        clear dm2wfe 
    end
    dm2mask=fitsread( [prefix 'dm2mask.fits']);
    wavefront = prop_multiply(wavefront, custom_pad(dm2mask, n));
    clear dm2mask 
end

wavefront = prop_propagate(wavefront, d_dm2_oap3, 'surface_name','OAP3');
wavefront = prop_lens(wavefront, fl_oap3);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP3_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap3/2);end

wavefront = prop_propagate(wavefront, d_oap3_fold3, 'surface_name','FOLD_3');
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_FOLD3_phase_error_V1.0.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_fold3/2);end

wavefront = prop_propagate(wavefront, d_fold3_oap4, 'surface_name','OAP4');
wavefront = prop_lens(wavefront, fl_oap4);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP4_phase_error_V1.0.fits'], 'wavefront');  end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap4/2);  end

wavefront = prop_propagate(wavefront, d_oap4_pupilmask,'surface_name', 'PUPIL_MASK');	% flat/reflective shaped pupil
if ( is_spc && use_pupil_mask ~= 0 )
    pupil_mask = custom_pad(fitsread( pupil_mask_file),n);
    if ( mask_x_shift_pupdiam ~=0 || mask_y_shift_pupdiam ~=0  || mask_x_shift_m ~=0 || mask_y_shift_m ~=0)
        %shift SP mask by FFTing it, applying tilt, and FFTing back
        x = repmat( ([1:n]-n/2-1) / (pupil_diam_pix/2), n, 1 );
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
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_PUPILMASK_phase_error_V1.0.fits'], 'wavefront'); end

% while at a pupil, use more padding to provide better sampling at FPM
diam = 2 * prop_get_beamradius(wavefront);
wavefront = prop_end(wavefront, 'noabs');

n = n_to_fpm;
wavefront0 = custom_pad(wavefront,n);
wavefront = prop_begin(diam, lambda_m, n, 'beam_diam_fraction', pupil_diam_pix/n);
wavefront.wf = prop_shift_center(wavefront0);
clear wavefront0 

wavefront = prop_propagate(wavefront, d_pupilmask_oap5, 'surface_name','OAP5');
wavefront = prop_lens(wavefront, fl_oap5);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP5_phase_error_V1.0.fits'], 'wavefront'); end
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
        occ = complex(fitsread(occulter_file_r),fitsread(occulter_file_i));
        if(is_fpm_grid)
            occ = occ.*fpm_grid;
%             figure(123); imagesc(abs(occ)); axis xy equal tight; colorbar; drawnow;
        end
        %--DEBUGGING
        if angle(occ(1,1)) ~= 0
            occ = occ.*exp(-1j*angle(occ(1,1))); %--Standardize the phase of the masks to be 0 for the outer glass part.
        end
        wavefront = prop_multiply(wavefront, custom_pad(occ,n));
        clear occ
    elseif ( is_spc )
        % super-sample FPM
        wavefront0 = custom_pad(ifftshift(fft2(wavefront.wf)), n_mft);  % to virtual pupil
        fpm = fitsread( fpm_file);
        fpm = custom_pad(fpm, size(fpm,1)+mod(size(fpm,1),2)); % make it even; otherwise the phase part is incorrect
        nfpm = size(fpm,1);
        fpm_sampling_lam = fpm_sampling_lam0 * lambda0_m / lambda_m;
        wavefront0 = mft2(wavefront0, fpm_sampling_lam, pupil_diam_pix, nfpm, -1); % MFT to highly-sampled focal plane
        wavefront0 = wavefront0.*fpm;
        clear fpm
        wavefront0 = mft2(wavefront0, fpm_sampling_lam, pupil_diam_pix, n, +1);  % MFT to virtual pupil
        wavefront.wf = ifft2(ifftshift(wavefront0));% back to normally-sampled focal plane
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
    wavefront0 = fftshift(fft2( wavefront.wf ));   % to virtual pupil
    wavefront0 = custom_pad( wavefront0, n_in );
    wavefront0 = mft2( wavefront0, dx_pinhole_lamD, pupil_diam_pix, n_out, -1 );		% MFT to highly-sampled focal plane
    p = (radpix(n_out,n_out) .* dx_pinhole_diam_m) <= pinhole_diam_m/2.0;
    wavefront0 = wavefront0 .* p;
    wavefront0 = mft2( wavefront0, dx_pinhole_lamD, pupil_diam_pix, n, +1 );	% MFT back to virtual pupil
    wavefront.wf = ifft2(ifftshift(wavefront0)); % back to normally-sampled focal plane
    clear wavefront0 
end

wavefront = prop_propagate(wavefront, d_fpm_oap6-fpm_z_shift_m, 'surface_name','OAP6');
wavefront = prop_lens(wavefront, fl_oap6);
if ( use_errors && ~end_at_fpm_exit_pupil );   wavefront =prop_errormap(wavefront, [map_dir 'wfirst_phaseb_OAP6_phase_error_V1.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap6/2); end

wavefront = prop_propagate(wavefront, d_oap6_lyotstop, 'surface_name','LYOT_STOP');
% while at a pupil, switch grid size (smaller for HLC for better speed, larger for SPC for better sampling)
diam = 2 * prop_get_beamradius(wavefront);
[wavefront, sampling_m ]= prop_end(wavefront,  'noabs');

n = n_from_lyotstop;
wavefront = custom_pad(wavefront,n);
wavefront0 = wavefront;

if  ~isempty(output_field_rootname)
    lams = num2str(lambda_m*1e6, '%6.4f');
    pols = ['polaxis'   num2str(polaxis,2)];
    wavefront = custom_pad(wavefront, output_dim); %
    fitswrite(real(wavefront), [output_field_rootname '_' lams 'um_' pols '_real.fits']);
    fitswrite(imag(wavefront), [output_field_rootname '_' lams 'um_' pols '_imag.fits']);
end

if ( end_at_fpm_exit_pupil )
    if getWFE
        pupil = fitsread( pupil_file);
        wavefront = custom_pad(wavefront.*custom_pad(pupil,n), pupil_diam_pix); % if to get wfe
    end
    return
end

wavefront = prop_begin (diam, lambda_m, n, 'beam_diam_fraction', pupil_diam_pix/n);
wavefront.wf = prop_shift_center(wavefront0);
clear wavefront0 

if ( use_lyot_stop )
    lyot = fitsread( lyot_stop_file );
    lyot = custom_pad( lyot, n );    
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

if ( use_pupil_lens  || pinhole_diam_m ~=0 ); wavefront =prop_circular_aperture(wavefront, 1.1, 'norm'); end

wavefront = prop_propagate(wavefront, d_lyotstop_oap7, 'surface_name','OAP7');
wavefront = prop_lens(wavefront, fl_oap7);
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP7_phase_error_V1.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_oap7/2); end

wavefront = prop_propagate(wavefront, d_oap7_fieldstop, 'surface_name','FIELD_STOP');

if ( use_field_stop && is_hlc )
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
end

wavefront = prop_propagate(wavefront, d_fieldstop_oap8, 'surface_name','OAP8');
wavefront =prop_lens(wavefront, fl_oap8);
if ( use_errors );   prop_errormap(wavefront,[map_dir 'wfirst_phaseb_OAP8_phase_error_V1.0.fits'], 'wavefront'); end
if ( use_aperture );   prop_circular_aperture(wavefront, diam_oap8/2);end

wavefront = prop_propagate(wavefront, d_oap8_filter, 'surface_name','filter');
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_FILTER_phase_error_V1.0.fits'], 'wavefront'); end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_filter/2); end

wavefront = prop_propagate(wavefront, d_filter_lens, 'surface_name','LENS');
if ( use_pupil_lens == 0 && use_defocus_lens == 0 && defocus == 0 )
    % use imaging lens to create normal focus
    wavefront = prop_lens(wavefront, fl_lens);
    if ( use_errors );   prop_errormap(wavefront,[map_dir 'wfirst_phaseb_LENS_phase_error_V1.0.fits'], 'wavefront'); end
else
    if ( use_pupil_lens )       % use pupil imaging lens
        wavefront = prop_lens(wavefront, fl_pupillens);
        if ( use_errors );   prop_errormap(wavefront,[map_dir 'wfirst_phaseb_PUPILLENS_phase_error_V1.0.fits'], 'wavefront'); end
    else                        % table is waves P-V @ 575 nm
        z4_pv_waves = [ -9.0545 -8.5543 -8.3550 -8.0300 -7.54500 -7.03350 -6.03300 -5.03300 -4.02000  ...
                        -2.51980 0.00000000 3.028000 4.95000 6.353600 8.030000 10.10500 12.06000 14.06000 ...
                        20.26000 28.34000 40.77500 56.65700 ]; 
        fl_defocus_lens = [ 5.09118 1.89323 1.54206 1.21198 0.914799 0.743569 0.567599 0.470213 0.406973  ...
                            0.350755 0.29601868 0.260092 0.24516 0.236606 0.228181 0.219748 0.213278 0.207816 ...
                            0.195536 0.185600 0.176629 0.169984 ];
	% subtract ad-hoc function to make z4 vs f_length more accurately spline interpolatible
	f = fl_defocus_lens / 0.005;
	f0 = 59.203738;
	z4t = z4_pv_waves - (0.005*(f0-f-40))./f.^2/0.575e-6;
        if ( use_defocus_lens ~=0 )    	% defocus lens (1-4)
            % use one of 4 defocusing lenses
            defocus = [ 18.0, 9.0, -4.0, -8.0 ];	% waves P-V @ 575 nm
            z4x = interp1( z4_pv_waves, z4t, defocus, 'spline' );
            lens_fl = interp1( z4t, fl_defocus_lens, z4x, 'spline' );
            wavefront =prop_lens(wavefront, lens_fl(use_defocus_lens));
            if ( use_errors );   wavefront =prop_errormap(wavefront,[map_dir 'wfirst_phaseb_DEFOCUSLENS' ...
                    num2str(use_defocus_lens,2)  '_phase_error_V1.0.fits'], 'wavefront');
            end
	    defocus = defocus(use_defocus_lens);
        else
            % specify amount of defocus (P-V waves @ 575 nm)
	    z4x = interp1( z4_pv_waves, z4t, defocus, 'spline' );
            lens_fl = interp1( z4t, fl_defocus_lens, z4x, 'spline' );
            wavefront =prop_lens(wavefront, lens_fl);
            if ( use_errors );   wavefront =prop_errormap(wavefront,[map_dir 'wfirst_phaseb_DEFOCUSLENS1_phase_error_V1.0.fits'], 'wavefront');end
        end
    end
end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_lens/2); end

wavefront = prop_propagate(wavefront, d_lens_fold4, 'surface_name','FOLD_4');
if ( use_errors );   wavefront = prop_errormap(wavefront,[map_dir 'wfirst_phaseb_FOLD4_phase_error_V1.1.fits'], 'wavefront');end
if ( use_aperture );   wavefront = prop_circular_aperture(wavefront, diam_fold4/2);end

if ( defocus ~= 0 || use_defocus_lens ~= 0 ) 
	if ( abs(defocus) <= 4 )
		wavefront = prop_propagate(wavefront, d_fold4_image, 'surface_name','IMAGE','TO_PLANE');
	else
		wavefront = prop_propagate(wavefront, d_fold4_image, 'surface_name','IMAGE');
	end
else
	wavefront = prop_propagate(wavefront, d_fold4_image, 'surface_name','IMAGE');
end

[wavefront, sampling_m] = prop_end(wavefront, 'noabs');

if ( final_sampling_lam0 ~=0 || final_sampling_m ~=0)
    if  final_sampling_m ~=0
        mag = sampling_m / (final_sampling_m);
        sampling_m = final_sampling_m;
    elseif final_sampling_lam0 ~=0
        mag = (pupil_diam_pix/n) / final_sampling_lam0 * (lambda_m/lambda0_m);
        sampling_m = sampling_m / mag;
    end
    wavefront = prop_magnify( wavefront, mag, 'size_out',output_dim,'amp_conserve');
else
    wavefront = custom_pad(wavefront, output_dim);
end

return
end
