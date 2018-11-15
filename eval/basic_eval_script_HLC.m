
flagPlot = true; %--Plot results
flagSave = false; %--Save FITS files

%% Set up Paths for FALCO and PROPER


%% Define file names to load
runLabel = '/full/path/to/data/Series0015_Trial0009_HLC_LUVOIRA5_3DM64_z3.24_IWA2.7_OWA10_11lams550nm_BW10_plannedEFC'; 

fn_init_ws = [runLabel, '_config.mat']; %--configuration file
fn_snippet_ws = [runLabel, '_snippet.mat']; %--post-trial data in structure named 'out'

fn_brief = fn_init_ws; 

%% Load final DM settings and assign to mp structure
load(fn_snippet_ws)

mp.dm1.V = out.DM1V;
mp.dm2.V = out.DM2V;
mp.dm8.V = out.DM8V;
mp.dm9.V = out.DM9V;

% % whichItr = out.Nitr+1;%54;%6;%5;
% % DM1V = out.dm1.Vall(:,:,whichItr);
% % DM2V = out.dm2.Vall(:,:,whichItr);
% % if(whichItr==out.Nitr+1);  extra = -1; end
% % DM8V = out.dm8.Vall(:,whichItr+extra);
% % DM9V = out.dm9.Vall(:,whichItr+extra); 

%% Increase Spatial Resolution in FP4 if desired
mp.F4.res = 10;

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

if(flagSave);  
    fitswrite(DM8surf,[fn_prefix,'data_DM8_surf.fits']);
    fitswrite(DM9surf,[fn_prefix,'data_DM9_surf.fits']);
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

I0f = falco_get_summed_image(mp);

if(flagPlot)
    figure(502); imagesc(log10(I0f),[-10.5 -8]); axis xy equal tight; colorbar;
    title('(Full Model: Normalization Check Using Starting PSF)'); drawnow;
end

if(flagSave);  fitswrite(I0f,[fn_prefix,'data_Inorm_broadband.fits']); end
