%   Copyright 2020 California Institute of Technology
% ------------------------------------------------------------------
%
% PUBLIC RELEASE VERSION (no ITAR error maps)
%
% Phase C, Version 0.9, 15 June 2020, John Krist; beta-test version for mask evaluation
% Phase C, Version 0.9.1, 29 June 2020, John Krist; rotated SPC pupils & Lyot stops
% Phase C, Version 0.9.2, 7 July 2020, John Krist; updated HLC occulter
% Phase C, Version 1.0, 1.0.1, 22 Sept 2020, John Krist
% Phase C, Version 1.1, 3 Dec 2020, John Krist: added field-dependent phase term at image
%     plane to match E-field from full model due to pupil offset
% Phase C, Version 1.2, 2 Feg 2021, John Krist: switched to using MFTs to propagate to/from HLC FPM;
%     updated the list of available HLC FPM wavelengths; removed field-dependent phase term 
%     (now subtracting it from the full model instead); read in list of available HLC FPM files;
%     added spc-mswc mode; added modes to allow spc-wide or spc-mswc to use band 1; removed the
%     input_field_rootname and polaxis options
% Phase C, Version 1.2.3, 15 July 2021, John Krist: added dm_sampling_m as optional parameter
% Based on translation by Hanying Zhou of John Krist's Phase B IDL model; converted
% to Phase C by John Krist  

function [wavefront,  sampling_m]= roman_phasec_compact(lambda_m, output_dim0, optval)

%   "output_dim" is used to specify the output dimension in pixels at the final image plane.  
%   The computational grid sizes are hardcoded for each coronagraph.

data_dir = '/home/krist/afta/phasec/roman_phasec_v1.2.2/phasec_data';	% no trailing '/'

if ( isfield(optval,'data_dir') )
    data_dir = optval.data_dir;
    if(strcmp(data_dir(end),'/') || strcmp(data_dir(end),'\'))
        data_dir = data_dir(1:end-1); %--Remove the trailing slash for this function.
    end
end

cor_type = 'hlc';           	%   'hlc', 'spc-spec_band3', 'spc-spec_band2', 'spc-wide'
source_x_offset = 0;		%   source offset in lambda0_m/D radians
source_y_offset = 0;
use_hlc_dm_patterns = 0;	%   use Dwight-generated HLC default DM wavefront patterns? 1 or 0
use_dm1 = 0;                	%   use DM1? 1 or 0
use_dm2 = 0;                	%   use DM2? 1 or 0
dm_sampling_m = 0.9906e-3;    	%   actuator spacing in meters; default is 1 mm
dm1_m = zeros(48,48);
dm1_xc_act = 23.5;          	%   for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
dm1_yc_act = 23.5;
dm1_xtilt_deg = 0;   		%   tilt around X axis
dm1_ytilt_deg = 9.65;		%   effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
dm1_ztilt_deg = 0;
dm2_m = zeros(48,48);
dm2_xc_act = 23.5;		
dm2_yc_act = 23.5;
dm2_xtilt_deg = 0;   
dm2_ytilt_deg = 9.65;
dm2_ztilt_deg = 0;
hlc_dm1_file = '';
hlc_dm2_file = '';
use_fpm = 1;
final_sampling_lam0 = 0;	%   final sampling in lambda0/D
output_dim = output_dim0;	%  dimensions of output in pixels (overrides output_dim0)

if (exist('optval','var')==1) 
	if ( isfield(optval,'cor_type') );  cor_type = optval.cor_type; end
	if ( isfield(optval,'use_fpm') );  use_fpm = optval.use_fpm; end
end

is_hlc = 0;
is_spc = 0;

if strcmp(cor_type,'hlc') 
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
    pupil_diam_pix = 309.0;   % Y pixel diameter in pixels
    pupil_file = [file_directory 'pupil.fits'];
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
                        disp('Error in roman_phasec_compact: requested wavelength not within 0.1 nm of nearest available FPM wavelength')
                        error(['   requested (um) = ', num2str(lam_um), '  closest available (um) = ', numtostr(fpm_lam_um(w))])
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
    n_small = 1024;	%   gridsize in non-critical areas
    n_big = 1024;	%   gridsize to FPM (propagation to/from FPM handled by MFT)
elseif ( strcmp(cor_type,'spc-spec') || strcmp(cor_type,'spc-spec_band3') || strcmp(cor_type,'spc-spec_band2') ) 
    is_spc = 1;
    file_dir = [data_dir '/spc_20200617_spec/'];        %   must have trailing "/"
    pupil_diam_pix  = 1000.0;
    pupil_file      = [file_dir  'pupil_SPC-20200617_1000.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200617_1000_rounded9.fits'];
    fpm_sampling_lam0divD = 0.05;	%   sampling in lambda0/D of FPM mask
    fpm_file        = [file_dir  'fpm_0.05lamD.fits'];
    if ( strcmp(cor_type,'spc-spec_band2') )
	fpm_lam0_m = 0.66e-6;
    else
	fpm_lam0_m = 0.73e-6;
    end
    lambda0_m = fpm_lam0_m; 
    lyot_stop_file  = [file_dir  'LS_SPC-20200617_1000.fits'];
    n_small = 2048;             %   gridsize in non-critical areas
    n_big = 1400;               %   gridsize to FPM (propagation to/from FPM handled by MFT)
elseif ( strcmp(cor_type,'spc-spec_rotated') )
    is_spc = 1;
    file_dir = [data_dir '/spc_20200628_specrot/'];      %   must have trailing "/"
    pupil_diam_pix  = 1000.0;
    pupil_file      = [file_dir  'pupil_SPC-20200628_1000.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200628_1000_derotated.fits'];
    fpm_sampling_lam0divD = 0.05;	%   sampling in lambda0/D of FPM mask
    fpm_file        = [file_dir  'FPM_SPC-20200628_res20.fits'];
    fpm_lam0_m = 0.73e-6;
    lambda0_m = fpm_lam0_m; 
    lyot_stop_file  = [file_dir  'LS_SPC-20200628_1000.fits'];
    n_small = 2048;             %   gridsize in non-critical areas
    n_big = 1400;               %   gridsize to FPM (propagation to/from FPM handled by MFT)
elseif ( strcmp(cor_type,'spc-wide') || strcmp(cor_type,'spc-wide_band4') || strcmp(cor_type,'spc-wide_band1') )
    is_spc = 1;
    file_dir = [data_dir '/spc_20200610_wfov/'];        %   must have trailing "/"
    pupil_diam_pix = 1000.0;      % Y axis pupil diameter in pixels
    pupil_file = [file_dir  'pupil_SPC-20200610_1000.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200610_1000_rounded9_gray.fits'];
    fpm_sampling_lam0divD = 0.1;       % sampling in lambda0/D of FPM mask
    fpm_file = [file_dir   'FPM_SPC-20200610_0.1_lamc_div_D.fits'];
    if ( strcmp(cor_type,'spc-wide_band1') )
    	fpm_lam0_m = 0.575e-6;
    else
    	fpm_lam0_m = 0.825e-6;
    end
    lambda0_m = fpm_lam0_m;       % FPM scaled for this central wavelength
    lyot_stop_file = [file_dir  'LS_SPC-20200610_1000.fits'];
    n_small = 2048;             %   gridsize in non-critical areas
    n_big = 1400;               %   gridsize to FPM (propagation to/from FPM handled by MFT)
elseif ( strcmp(cor_type,'spc-mswc') || strcmp(cor_type,'spc-mswc_band4') || strcmp(cor_type,'spc-mswc_band1') )
    is_spc = 1;
    file_dir = [data_dir '/spc_20200623_mswc/'];        %   must have trailing "/"
    pupil_diam_pix = 982.0;      % Y axis pupil diameter in pixels
    pupil_file = [file_dir  'pupil_SPC-20200623_982_rotated.fits'];
    pupil_mask_file = [file_dir  'SPM_SPC-20200623_982_rounded9_gray.fits'];
    fpm_sampling_lam0divD = 0.1;       % sampling in lambda0/D of FPM mask
    fpm_file = [file_dir   'FPM_SPC-20200623_0.1_lamc_div_D.fits'];
    if ( strcmp(cor_type,'spc-wide_band1') )
    	fpm_lam0_m = 0.575e-6;
    else
    	fpm_lam0_m = 0.825e-6;
    end
    lambda0_m = fpm_lam0_m;       % FPM scaled for this central wavelength
    lyot_stop_file = [file_dir  'LS_SPC-20200623_982_rotated.fits'];
    n_small = 2048;             %   gridsize in non-critical areas
    n_big = 1400;               %   gridsize to FPM (propagation to/from FPM handled by MFT)
else
    disp([ '  ERROR: Unsupported cor_type: ' cor_type ])
    return
end

if ( exist('optval','var')==1 && numel(fieldnames(optval)) ~=0 ) 
	if ( isfield(optval,'source_x_offset') );  source_x_offset = optval.source_x_offset;end 
	if ( isfield(optval,'source_y_offset') );  source_y_offset = optval.source_y_offset;end
	if ( isfield(optval,'use_hlc_dm_patterns') );  use_hlc_dm_patterns = optval.use_hlc_dm_patterns;end
	if ( isfield(optval,'use_dm1') );        use_dm1 = optval.use_dm1;end
	if ( isfield(optval,'dm_sampling_m') );  dm_sampling_m = optval.dm_sampling_m;end
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
	if ( isfield(optval,'output_dim') );    output_dim = optval.output_dim; end
end

if use_hlc_dm_patterns && ~strcmp(cor_type,'hlc') && ~strcmp(cor_type,'hlc_band1')
    disp('ERROR: Can only utilize use_hlc_dm_patterns with Band 1 HLC');
    return
end

diam_at_dm1 = 0.0463;
d_dm1_dm2 = 1.0;

n = n_small;        	% start off with less padding

wavefront = prop_begin(diam_at_dm1,lambda_m, n,'beam_diam_fraction', pupil_diam_pix/n);%

pupil = fitsread(pupil_file);
wavefront = prop_multiply(wavefront, custom_pad(pupil,n));
clear pupil

wavefront  = prop_define_entrance(wavefront);

if ( source_x_offset ~= 0.0 || source_y_offset ~= 0.0 )
    %   compute tilted wavefront to offset source by source_x_offset,source_y_offset lambda0_m/D
    x  = repmat( ((1:n)-n/2-1)/(pupil_diam_pix/2),n,1);
    tilt = -1i*pi*(x*source_x_offset + x'*source_y_offset)*lambda0_m / lambda_m; % shift in lam/D at lam0,
    wavefront = prop_multiply(wavefront, exp(tilt));
    clear x tilt
end

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

if ( use_dm1 )
    wavefront = prop_dm(wavefront, dm1, dm1_xc_act, dm1_yc_act, dm_sampling_m,...
        'xtilt',dm1_xtilt_deg, 'ytilt', dm1_ytilt_deg, 'ztilt',dm1_ztilt_deg);
end

wavefront = prop_propagate(wavefront, d_dm1_dm2, 'surface_name','DM2');

if ( use_dm2 )
    wavefront = prop_dm(wavefront, dm2, dm2_xc_act, dm2_yc_act, dm_sampling_m,...
        'xtilt',dm2_xtilt_deg, 'ytilt', dm2_ytilt_deg, 'ztilt',dm2_ztilt_deg);
end

if ( is_hlc )
    dm2mask = fitsread([file_directory 'dm2mask.fits']);
    wavefront = prop_multiply (wavefront, custom_pad(dm2mask, n));
    clear dm2mask
end

wavefront = prop_propagate (wavefront, -d_dm1_dm2, 'surface_name','back to DM1');

[wavefront, sampling_m]= prop_end (wavefront, 'NOABS');
n = n_big;
wavefront = custom_pad(wavefront,n);

if ( is_spc ) 
    pupil_mask = fitsread(pupil_mask_file);
    wavefront = wavefront .* custom_pad(pupil_mask,n);
    clear pupil_mask 
end

if ( use_fpm )
    if ( is_hlc ) 
    	wavefront = wavefront * fpm_array(1,1);       % apply uniform amplitude & phase from FPM clear area (assumes clear area is uniform) 
    	nfpm = size(fpm_array,1);
    	fpm_sampling_lamdivD = fpm_sampling_lam0divD * fpm_lam0_m / lambda_m;        % FPM sampling at current wavelength in lambda_m/D
    	wavefront_fpm = mft2(wavefront, fpm_sampling_lamdivD, pupil_diam_pix, nfpm, -1);       % MFT to highly-sampled focal plane
    	wavefront_fpm = wavefront_fpm .* fpm_mask .* (fpm_array - 1);           % subtract field inside FPM region, add in FPM-multiplied region
    	wavefront_fpm = mft2(wavefront_fpm, fpm_sampling_lamdivD, pupil_diam_pix, n, +1);       % MFT to pupil (Lyot stop)
    	wavefront = wavefront + wavefront_fpm;
    	clear wavefront_fpm
    elseif ( is_spc ) 
    	fpm = fitsread(fpm_file);
    	fpm = custom_pad(fpm, size(fpm,1)+mod(size(fpm,1),2));
    	nfpm = size(fpm,1);
    	fpm_sampling_lamdivD = fpm_sampling_lam0divD * fpm_lam0_m / lambda_m;
    	wavefront = mft2(wavefront, fpm_sampling_lamdivD, pupil_diam_pix, nfpm,  -1); %   MFT to highly-sampled focal plane
    	wavefront = wavefront .* fpm;
    	clear fpm 
    	wavefront = mft2(wavefront, fpm_sampling_lamdivD, pupil_diam_pix,fix(pupil_diam_pix), +1);
    end
end

n = n_small;
wavefront = custom_pad(wavefront,n);
if ( ~strcmp(cor_type, 'none' ) )
    lyot = fitsread (lyot_stop_file);
    wavefront = wavefront .* custom_pad(lyot,n);
    clear lyot  
end

wavefront = fftshift(ifft2(ifftshift(wavefront)))*n;  %   FFT to final focus

% resample, if requested

if ( final_sampling_lam0 ~= 0 )
    mag = (pupil_diam_pix/n) / final_sampling_lam0 * (lambda_m/lambda0_m);
    wavefront = prop_magnify( wavefront, mag, 'size_out',output_dim, 'amp_conserve');
else
    wavefront = custom_pad(wavefront, output_dim);
end

sampling_m = 0.0;

return

end

