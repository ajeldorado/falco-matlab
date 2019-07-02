% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Simple script to load the final result of an HLC design run and save out
% the DM maps, FPM at several wavelengths, and the PSF with and without
% jitter or stellar diameter.
%
% Written on 2019-05-29 by A.J. Riggs.

clear all

flagPlot = true; %--Plot results
flagSave = false;%true;%false; %--Save FITS files

%% Set up Paths for FALCO and PROPER

%--Library locations
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = '~/Repos/falco-matlab/data/brief/'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = '~/Repos/falco-matlab/data/ws/'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

%%
% fn_prefix = '/home/ajriggs/Repos/falco-matlab/data/brief/';
fn_prefix = '/home/ajriggs/Downloads/hlc20190523/';


%% Define file names to load
% runLabel = '~/Repos/falco-matlab/data/brief/Series0022_Trial0006_HLC_WFIRST180718_3DM48_z1_IWA2.7_OWA9_5lams575nm_BW10_plannedEFC';  %--good trials are 2,3,6,9
runLabel = '~/Repos/falco-matlab/data/brief/Series0867_Trial5309_HLC_LUVOIRAfinal_2DM64_z3_IWA2.8_OWA10_5lams500nm_BW10_gridsearchEFC';
fn_init_ws = [runLabel, '_config.mat']; %--configuration file
fn_snippet_ws = [runLabel, '_snippet.mat']; %--post-trial data in structure named 'out'

fn_brief = fn_init_ws; 

%% Load final DM settings and assign to mp structure

load(fn_init_ws)
load(fn_snippet_ws)
figure(70); semilogy(0:mp.Nitr,out.InormHist,'Linewidth',3); set(gca,'Fontsize',20); grid on;

%%
mp.dm1.V = out.DM1V;
mp.dm2.V = out.DM2V;
mp.dm8.V = out.DM8V;
mp.dm9.V = out.DM9V;

%% Increase Spatial Resolution in FP4 if desired
mp.Fend.res = 10;

mp.Nsbp = 9;
mp.Nwpsbp = 1;

%% Save to new config file that includes the final DM settings as the new starting point
fn_config_updated = [fn_init_ws(1:end-4) '_updated.mat'];
save(fn_config_updated)

%% Run initialization function to complete the workspace

[mp,~] = falco_init_ws(fn_config_updated);

%% Plot Thickness Profiles for FPM (DM8 and DM9)
modelType = 'compact';
t_diel_bias = mp.t_diel_bias_nm*1e-9; % bias thickness of dielectric [m]

%--Generate thickness profiles of each layer
DM8surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm8,modelType); %--Metal layer profile [m]
DM9surf = t_diel_bias + falco_gen_HLC_FPM_surf_from_cube(mp.dm9,modelType); %--Dielectric layer profile [m]

figure(8); imagesc(DM8surf*1e9); axis xy equal tight; colorbar;
figure(9); imagesc(DM9surf*1e9); axis xy equal tight; colorbar; drawnow;

if(flagSave)
    fitswrite(DM8surf,[fn_prefix,'data_DM8_surf.fits']);
    fitswrite(DM9surf,[fn_prefix,'data_DM9_surf.fits']);
end

%% Make the FPM at each wavelength
N = size(DM8surf,1);
FPMcube = zeros(N,N,mp.Nsbp);
wi = 1;
for si=1:mp.Nsbp
    FPMcube(:,:,si) = falco_gen_HLC_FPM_complex_trans_mat(mp,si,wi,'full');
    FPMcube(:,:,si) = FPMcube(:,:,si)/FPMcube(1,1,si);
    figure(21); imagesc(abs(FPMcube(:,:,si))); axis xy equal tight; colorbar; drawnow;
    figure(22); imagesc(angle(FPMcube(:,:,si))); axis xy equal tight; colorbar; drawnow;
end

if(flagSave)
    fitswrite(real(FPMcube),[fn_prefix, 'realFPMcube.fits' ])
    fitswrite(imag(FPMcube),[fn_prefix, 'imagFPMcube.fits' ])
end

%% Save out the pupil and Lyot stop

figure(23); imagesc(mp.P1.full.mask); axis xy equal tight;
figure(24); imagesc(mp.P4.full.mask); axis xy equal tight;


if(flagSave)
    fitswrite(mp.P1.full.mask,[fn_prefix, 'pupil.fits' ])
    fitswrite(mp.P4.full.mask,[fn_prefix, 'LyotStop.fits' ])
end

%% Compute and plot the DM1 and DM2 surfaces
lambda0 = mp.lambda0;

if(any(mp.dm_ind==1)) 
    
    if(flagSave);  fitswrite(mp.dm1.V,[fn_prefix,'data_DM1_V.fits']);  end

    DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  
    if(flagSave);  fitswrite(DM1surf,[fn_prefix,'data_DM1_surf.fits']); end
end

if(any(mp.dm_ind==2)) 
    if(flagSave);  fitswrite(mp.dm2.V,[fn_prefix,'data_DM2_V.fits']);  end

    DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  
    if(flagSave);  fitswrite(DM2surf,[fn_prefix,'data_DM2_surf.fits']);  end
end

if(flagPlot)
    figure(3); imagesc(DM1surf*1e9); axis xy equal tight; colorbar;
    figure(4); imagesc(DM2surf*1e9); axis xy equal tight; colorbar; drawnow;
end

%% Plot the PSF
% 
% I0f = falco_get_summed_image(mp);
% 
% if(flagPlot)
%     figure(502); imagesc(log10(I0f),[-10.5 -8]); axis xy equal tight; colorbar;
%     title('(Full Model: Normalization Check Using Starting PSF)'); drawnow;
% end
% 
% if(flagSave);  fitswrite(I0f,[fn_prefix,'data_Inorm_broadband.fits']); end

%% Get PSFs
%--Compute coherent-only image
tic; fprintf('Generating the PSF of only coherent light... ')
mp.full.pol_conds = 0;%10; %--Which polarization states to use when creating an image.
ImCoh = falco_get_summed_image(mp);
fprintf('done. Time = %.2f s\n',toc)

%--Compute coherent+incoherent light image
tic; fprintf('Generating the PSF with polarization aberrations and stellar size... ')
mp.full.pol_conds = 0;%[-2,-1,1,2]; %--Which polarization states to use when creating an image.
mp.full.TTrms = 0; % [mas]
mp.full.Dstar = 1; % [mas]
mp.full.Dtel = 15.0; % [meters]
mp.full.TipTiltNacross = 7;%5; 
ImBoth = falco_get_summed_image_TipTiltPol(mp);
fprintf('done. Time = %.2f s\n',toc)

%--Subtract to get the incoherent component
ImInco = ImBoth-ImCoh;
ImInco(ImInco<0) = 0;

if(mp.flagPlot)
    figure(81); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(ImCoh),[-10,-7]); axis xy equal tight; colorbar; set(gcf,'Color','w'); set(gca,'Fontsize',16)
        rectangle('Position',[-1/2,-1/2,1,1]*mp.Fend.corr.Rin*2,'Curvature',[1,1], 'Linewidth',2,'EdgeColor','w','LineStyle','--'); drawnow;% [x0 y0 w h]
        rectangle('Position',[-1/2,-1/2,1,1]*mp.Fend.corr.Rout*2,'Curvature',[1,1], 'Linewidth',2,'EdgeColor','w','LineStyle','--'); drawnow;% [x0 y0 w h]        
    figure(82); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(ImBoth),[-10,-7]); axis xy equal tight; colorbar; set(gcf,'Color','w'); set(gca,'Fontsize',16)
        rectangle('Position',[-1/2,-1/2,1,1]*mp.Fend.corr.Rin*2,'Curvature',[1,1], 'Linewidth',2,'EdgeColor','w','LineStyle','--'); drawnow;% [x0 y0 w h]
        rectangle('Position',[-1/2,-1/2,1,1]*mp.Fend.corr.Rout*2,'Curvature',[1,1], 'Linewidth',2,'EdgeColor','w','LineStyle','--'); drawnow;% [x0 y0 w h]  
    figure(83); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(ImInco),[-10,-7]); axis xy equal tight; colorbar; set(gcf,'Color','w'); set(gca,'Fontsize',16)
        rectangle('Position',[-1/2,-1/2,1,1]*mp.Fend.corr.Rin*2,'Curvature',[1,1], 'Linewidth',2,'EdgeColor','w','LineStyle','--'); drawnow;% [x0 y0 w h]
        rectangle('Position',[-1/2,-1/2,1,1]*mp.Fend.corr.Rout*2,'Curvature',[1,1], 'Linewidth',2,'EdgeColor','w','LineStyle','--'); drawnow;% [x0 y0 w h]  
end

if(flagSave);  fitswrite(ImCoh,[fn_prefix,'data_Inorm_BB.fits']); end
if(flagSave);  fitswrite(ImBoth,[fn_prefix,'data_Inorm_BBTT.fits']); end