%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

%   Copyright 2019, California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   Version 1.0, 3 Jan 2019, JEK
%   Version 1.0a, 9 Jan 2019, JEK
%     Changes: If reading in input field using input_field_rootname, the
%              user can specify polaxis, and the appropriate file will
%              be used.  By default, polaxis = 0.  This is the only
%              time polaxis can be specified for the compact model.
%   Version 1.0c, 13 February 2019, JEK 
%      Changes: Changed HLC design to HLC-20190210; added 'fpm_axis' parameter to specify
%                 axis of HLC FPM ('s' or 'p') due to FPM tilt; changed 'spc' coronagraph
%                 option to 'cor_type' to 'spc-ifs' and changed default SPC to SPC-20190130 and
%                 updated associated parameters; added 'spc-wide' option to 'cor_type' to use
%                 SPC-20191220.  Changed DM tilt to rotation around Y axis.
%   Version 1.2 , 13 March 2019, JEK 
%      Changes:   Changed option of cor_type='spc-ifs' to 'spc-ifs_short' for 660 nm FPM, 
%  		  'spc-ifs_long' for 730 nm FPM.  Rotated SPC masks to match orientation 
%  		  of pupil in HLC. Added Erkin's HLC design. dm_sampling set to 0.9906 mm.
%   Based on Zemax prescription "WFIRST_CGI_DI_LOWFS_Sep24_2018.zmx" by Hong Tang.
%
%  IDL version (original), John E Krist
%  Matlab translation, 3/14/2019, Hanying Zhou;
%  Minor changes to syntax, 4/3/2019, A.J. Riggs.

function [wavefront,  sampling_m]= model_compact_wfirst_phaseb(lambda_m, output_dim0, optval)

%   "output_dim" is used to specify the output dimension in pixels at the final image plane.  
%   The computational grid sizes are hardcoded for each coronagraph.

data_dir = '/home/krist/afta/phaseb/phaseb_data';	% no trailing '/'

if ( isfield(optval,'data_dir') )
    data_dir = optval.data_dir;
    if(strcmp(data_dir(end),'/') || strcmp(data_dir(end),'\'))
        data_dir = data_dir(1:end-1); %--Remove the trailing slash for this function.
    end
end

cor_type = 'hlc';           %   'hlc', 'spc', or 'none' (none = clear aperture, no coronagraph)
input_field_rootname = '';	%   rootname of files containing aberrated pupil
polaxis = 0;                %   polarization condition (only used with input_field_rootname)
source_x_offset = 0;		%   source offset in lambda0_m/D radians
source_y_offset = 0;
use_hlc_dm_patterns = 0;	%   use Dwight-generated HLC default DM wavefront patterns? 1 or 0
use_dm1 = 0;                %   use DM1? 1 or 0
use_dm2 = 0;                %   use DM2? 1 or 0
dm_sampling_m = 0.9906e-3;    %   actuator spacing in meters; default is 1 mm
dm1_xc_act = 23.5;          %   for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
dm1_yc_act = 23.5;
dm1_xtilt_deg = 0;   		%   tilt around X axis
dm1_ytilt_deg = 5.7;		%   effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
dm1_ztilt_deg = 0;
dm2_xc_act = 23.5;		
dm2_yc_act = 23.5;
dm2_xtilt_deg = 0;   
dm2_ytilt_deg = 5.7;
dm2_ztilt_deg = 0;
use_fpm  = 1;
fpm_axis = 'p';             %   HLC FPM axis: '', 's', 'p'
final_sampling_lam0 = 0;	%   final sampling in lambda0/D
output_dim = output_dim0;	%  dimensions of output in pixels (overrides output_dim0)

if(exist('optval','var')==1) %if exist('optval') 
	if ( isfield(optval,'cor_type') );  cor_type = optval.cor_type; end
	if ( isfield(optval,'fpm_axis') );  fpm_axis = optval.fpm_axis; end
end

if strcmp(cor_type,'hlc') 
    file_directory = [data_dir '/hlc_20190210/'];         %   must have trailing "/"
    prefix = [file_directory   'run461_'];
    pupil_diam_pix = 309.0;
    pupil_file = [prefix   'pupil_rotated.fits'];
    lyot_stop_file = [prefix  'lyot.fits'];
    lambda0_m = 0.575e-6;
    nlams =9 ;              % number of occ trans lambda provided
    bw = 0.1;
    lam_occ = ((1-bw/2):bw/(nlams-mod(nlams,2)):(1+bw/2))*lambda0_m; % wavelength at which occ trans provided
    %   find nearest matching FPM wavelength
    [~,wlam] = min(abs(lambda_m-lam_occ));
    occulter_file_r = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'real_rotated.fits'];
    occulter_file_i = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'imag_rotated.fits'];
    n_small = 1024;         %   gridsize in non-critical areas
    n_big = 2048;           %   gridsize to/from FPM
elseif  strcmp(cor_type,'hlc_erkin')
    file_directory = [data_dir '/hlc_20190206_v3/'];         % must have trailing "/"
    prefix = [file_directory  'dsn17d_run2_pup310_fpm2048_'];
    pupil_diam_pix = 310.0;
    pupil_file = [prefix  'pupil.fits'];
    lyot_stop_file = [prefix  'lyot.fits'];
    fpm_axis = 's';
    lambda0_m = 0.575e-6;
    nlams = 19 ;              % number of occ trans lambda provided
    bw = 0.1;
    lam_occ = ((1-bw/2):bw/(nlams-mod(nlams,2)):(1+bw/2))*lambda0_m; % wavelengths at which occ trans provideed
    
    [~,wlam] = min(abs(lambda_m-lam_occ)); % find nearest matching FPM wavelength
    occulter_file_r = [prefix  'occ_lam' num2str(lam_occ(wlam),5) 'theta6.69pol'   fpm_axis  '_real.fits'];
    occulter_file_i = [prefix  'occ_lam' num2str(lam_occ(wlam),5) 'theta6.69pol'   fpm_axis  '_imag.fits'];
    
    n_small = 1024;	% gridsize in non-critical areas
    n_big =2048;    % gridsize to/from FPM

elseif  contains(cor_type, 'spc-spec' )
    file_dir = [data_dir '/spc_20190130/'];        %   must have trailing "/"
    pupil_diam_pix  = 1000.0;
    pupil_file      = [file_dir  'pupil_SPC-20190130_rotated.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20190130_rotated.fits'];
    fpm_file        = [file_dir  'fpm_0.05lamdivD.fits'];
    lyot_stop_file  = [file_dir  'lyotstop_0.5mag.fits'];
    fpm_sampling_lam0 = 0.05;	%   sampling in lambda0/D of FPM mask
    lambda0_m = 0.73e-6;        %   FPM scaled for this central wavelength
    if ( contains(cor_type,'spc-spec_short' ));  lambda0_m = 0.66e-6; end 
    n_small = 2048;             %   gridsize in non-critical areas
    n_big = 1400;               %   gridsize to FPM (propagation to/from FPM handled by MFT)
elseif  strcmp(cor_type, 'spc-wide' )
    file_dir = [data_dir '/spc_20181220/'];        %   must have trailing "/"
    pupil_diam_pix = 1000.0;
    pupil_file      = [file_dir  'pupil_SPC-20181220_1k_rotated.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20181220_1000_rounded9_gray_rotated.fits'];
    fpm_file        = [file_dir  'fpm_0.05lamdivD.fits'];
    lyot_stop_file  = [file_dir  'LS_half_symm_CGI180718_Str3.20pct_38D91_N500_pixel.fits'];
    fpm_sampling_lam0 = 0.05;	%   sampling in lambda0/D of FPM mask
    lambda0_m = 0.825e-6;       %   FPM scaled for this central wavelength
    n_small = 2048;             %   gridsize in non-critical areas
    n_big = 1400;               %   gridsize to FPM (propagation to/from FPM handled by MFT)
elseif strcmp(cor_type, 'none' )
    file_directory = './hlc_20190210/';         % must have trailing "/"
    prefix = file_directory + 'run461_';
    pupil_file = prefix + 'pupil_rotated.fits';
    use_fpm = 0;
    pupil_diam_pix = 309.0;
%     n_default = 1024;
%     n_to_fpm = 1024;
%     n_from_lyotstop = 1024;
else
    disp([ '  ERROR: Unsupported cor_type: ' cor_type ])
    return
end

if( exist('optval','var')==1 && numel(fieldnames(optval)) ~=0 ) 
	if ( isfield(optval,'lam0') );           lambda0_m = optval.lam0 * 1.0e-6;end
	if ( isfield(optval,'lambda0_m') );      lambda0_m = optval.lambda0_m;end
	if ( isfield(optval,'input_field_rootname') );  input_field_rootname = optval.input_field_rootname;end
	if ( isfield(optval,'polaxis') );        polaxis = optval.polaxis;end
	if ( isfield(optval,'source_x_offset') );  source_x_offset = optval.source_x_offset;end 
	if ( isfield(optval,'source_y_offset') );  source_y_offset = optval.source_y_offset;end
	if ( isfield(optval,'use_hlc_dm_patterns') );  use_hlc_dm_patterns = optval.use_hlc_dm_patterns;end
	if ( isfield(optval,'use_dm1') );        use_dm1 = optval.use_dm1;end
	if ( isfield(optval,'dm1_m') );          dm1_m = optval.dm1_m;end
	if ( isfield(optval,'dm1_xc_act') );     dm1_xc_act = optval.dm1_xc_act;end
	if ( isfield(optval,'dm1_yc_act') );     dm1_yc_act = optval.dm1_yc_act;end
	if ( isfield(optval,'dm1_xtilt_deg') );  dm1_xtilt_deg = optval.dm1_xtilt_deg;end
	if ( isfield(optval,'dm1_ytilt_deg') );  dm1_ytilt_deg = optval.dm1_ytilt_deg;end
	if ( isfield(optval,'dm1_ztilt_deg') );  dm1_ztilt_deg = optval.dm1_ztilt_deg;end
	if ( isfield(optval,'use_dm2') );        use_dm2 = optval.use_dm2;end
	if ( isfield(optval,'dm2_m') );          dm2_m = optval.dm2_m;end
	if ( isfield(optval,'dm2_xc_act') );     dm2_xc_act = optval.dm2_xc_act;end
	if ( isfield(optval,'dm2_yc_act') );     dm2_yc_act = optval.dm2_yc_act;end
	if ( isfield(optval,'dm2_xtilt_deg') );  dm2_xtilt_deg = optval.dm2_xtilt_deg;end
	if ( isfield(optval,'dm2_ytilt_deg') );  dm2_ytilt_deg = optval.dm2_ytilt_deg;end
	if ( isfield(optval,'dm2_ztilt_deg') );  dm2_ztilt_deg = optval.dm2_ztilt_deg;end
	if ( isfield(optval,'final_sampling_lam0') );  final_sampling_lam0 = optval.final_sampling_lam0;end
    if ( isfield(optval,'use_fpm') );        use_fpm = optval.use_fpm;end
    if ( isfield(optval,'output_dim') );    output_dim = optval.output_dim; end
end

if( polaxis ~=0 && isempty(input_field_rootname) )  
	disp( '  ERROR: polaxis can only be defined when input_field_rootname is given')
	wavefront = 0;
	sampling_m = 0.0;
	return
end

diam_at_dm1 = 0.0463;
d_dm1_dm2 = 1.0;

n = n_small;        %   start off with less padding

wavefront = prop_begin(diam_at_dm1,lambda_m, n,'beam_diam_fraction', pupil_diam_pix/n);%

if  isempty(input_field_rootname)
    pupil = fitsread(pupil_file);
    wavefront = prop_multiply(wavefront, custom_pad(pupil,n));
    clear pupil
    %pupil = 0;
else
    lams = num2str(lambda_m*1e6, '%6.4f');
    pols = ['polaxis'  num2str(polaxis,2)];
    rval = fitsread ([input_field_rootname '_' lams 'um_' pols '_real.fits']);
    ival = fitsread ([input_field_rootname '_' lams 'um_' pols '_imag.fits']);
    wavefront	= prop_multiply(wavefront, custom_pad(complex(rval, ival),n));
    clear rval ival
end


wavefront  = prop_define_entrance(wavefront);

if(source_x_offset ~= 0.0 || source_y_offset ~= 0.0)
    %   compute tilted wavefront to offset source by source_x_offset,source_y_offset lambda0_m/D
    x  = repmat( ((1:n)-n/2-1)/(pupil_diam_pix/2),n,1);
    tilt = -1i*pi*(x*source_x_offset + x'*source_y_offset)*lambda0_m / lambda_m; % shift in lam/D at lam0,
    wavefront = prop_multiply(wavefront, exp(tilt));
    clear x tilt
end

if(use_dm1)
    wavefront = prop_dm(wavefront, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m,...
        'xtilt',dm1_xtilt_deg, 'ytilt', dm1_ytilt_deg, 'ztilt',dm1_ztilt_deg);
end

if( contains(cor_type,'hlc') && use_hlc_dm_patterns )
    dm1wfe = fitsread([prefix 'dm1wfe.fits']);
    wavefront = prop_add_phase(wavefront, custom_pad(dm1wfe, n));
    clear dm1wfe
end

wavefront = prop_propagate(wavefront, d_dm1_dm2, 'surface_name','DM2');

if(use_dm2)
    wavefront = prop_dm(wavefront, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m,...
        'xtilt',dm2_xtilt_deg, 'ytilt', dm2_ytilt_deg, 'ztilt',dm2_ztilt_deg);
end

if  contains(cor_type,'hlc' )
    if  use_hlc_dm_patterns
        dm2wfe = fitsread([prefix 'dm2wfe.fits']);
        wavefront = prop_add_phase( wavefront, custom_pad(dm2wfe, n));
        clear dm2wfe
    end
    dm2mask = fitsread([prefix 'dm2mask.fits']);
    wavefront = prop_multiply (wavefront, custom_pad(dm2mask, n));
    clear dm2wfe 
end

wavefront = prop_propagate (wavefront, -d_dm1_dm2, 'surface_name','back to DM1');

[wavefront, sampling_m]= prop_end (wavefront, 'noabs');

if contains(cor_type, 'spc')   
    pupil_mask = fitsread(pupil_mask_file);
    wavefront = wavefront .* custom_pad(pupil_mask,n);
    clear pupil_mask 
end

if contains(cor_type,'hlc' ) && use_fpm 
    n = n_big;
    wavefront = custom_pad(wavefront,n);
    wavefront = fftshift(fft2(ifftshift(wavefront)));  %   to focus
    occ = complex(fitsread(occulter_file_r),fitsread(occulter_file_i));
    wavefront = wavefront .* custom_pad(occ,n);
    clear occ 
    wavefront = fftshift(ifft2(ifftshift(wavefront)));  %   FFT to Lyot stop
elseif( contains(cor_type, 'spc') && use_fpm ) 
    n = n_big;
    wavefront = custom_pad(wavefront,n);
    fpm = fitsread(fpm_file);
    fpm = custom_pad(fpm, size(fpm,1)+mod(size(fpm,1),2));
    nfpm = size(fpm,1);
    fpm_sampling_lam = fpm_sampling_lam0 * lambda0_m / lambda_m;
    wavefront = mft2(wavefront,fpm_sampling_lam, pupil_diam_pix, nfpm,  -1); %   MFT to highly-sampled focal plane

    wavefront = wavefront .* fpm;
    clear fpm 
    pupil_diam_pix = pupil_diam_pix / 2;	%   Shrink pupil by 1/2
    wavefront = mft2(wavefront, fpm_sampling_lam, pupil_diam_pix,fix(pupil_diam_pix), +1);
end

n = n_small;
wavefront = custom_pad(wavefront,n);
if  ~strcmp(cor_type, 'none' )
    lyot = fitsread (lyot_stop_file);
    wavefront = wavefront .* custom_pad(lyot,n);
    clear lyot  
end

wavefront = fftshift(ifft2(ifftshift(wavefront)))*n;  %   FFT to final focus

if ( final_sampling_lam0 ~= 0 )
    mag = (pupil_diam_pix/n) / final_sampling_lam0 * (lambda_m/lambda0_m);
    wavefront = prop_magnify( wavefront, mag, 'size_out',output_dim, 'amp_conserve');
else
    wavefront = custom_pad(wavefront, output_dim);
end

% wavefront = circshift(rot90(wavefront,2),[1,1]);	%   rotate to same orientation as full model

sampling_m = 0.0;

return

end

