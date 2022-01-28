
flagPlot = true; %--Plot results
flagSave = false;%true;%false; %--Save FITS files

%% Set up Paths for FALCO and PROPER

%--Library locations
mp.path.falco = 'C:\CoronagraphSims\falco-matlab\';  %--Location of FALCO
mp.path.proper = 'C:\CoronagraphSims\FALCO\lib\PROPER\'; %--Location of the MATLAB PROPER library
% mp.path.cvx = '~/Documents/MATLAB/cvx/'; %--Location of MATLAB CVX

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
mp.path.config = 'C:\CoronagraphSims\falco-matlab\data\brief\'; %--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
mp.path.ws = 'C:\CoronagraphSims\falco-matlab\data\ws\'; % (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];

%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% addpath(genpath(mp.path.cvx)) %--Add CVX to MATLAB path

%%
fn_prefix = 'C:\CoronagraphSims\falco-matlab\data\brief\';

%% Define file names to load
runLabel = 'C:\CoronagraphSims\falco-matlab\data\brief\Series0010_Trial0003_vortex_LUVOIR_B_offaxis_1DM32_z1_IWA2_OWA15_12lams690nm_BW40_gridsearchEFC';

fn_init_ws = [runLabel, '_config.mat']; %--configuration file
fn_snippet_ws = [runLabel, '_snippet.mat']; %--post-trial data in structure named 'out'

fn_brief = fn_init_ws; 

%% Load final DM settings and assign to mp structure

load(fn_init_ws)
load(fn_snippet_ws)

%%
mp.dm1.V = out.DM1V;
mp.dm2.V = out.DM2V;

%% Change FP4 and other Properties if desired
% mp.Fend.res = 10;
% mp.Fend.FOV = 30; %For displaying the image of Fend
mp.Nsbp = 150;
mp.fracBW = 0.65;

%% Save to new config file that includes the final DM settings as the new starting point
fn_config_updated = [fn_init_ws(1:end-4) '_updated.mat'];
save(fn_config_updated)

%% Run initialization function to complete the workspace

[mp,~] = falco_init_ws(fn_config_updated);

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

% if(flagPlot)
%     figure(3); imagesc(DM1surf*1e9); axis xy equal tight; colorbar;
%     figure(4); imagesc(DM2surf*1e9); axis xy equal tight; colorbar; drawnow;
% end

%% Plot the PSF

% I0f = falco_get_summed_image(mp);
% 
% if(flagPlot)
%     figure; imagesc(log10(I0f),[-10 -3]); axis xy equal tight; colorbar;
%     title('(Full Model: Normalization Check Using Starting PSF)'); drawnow;
% end
% 
% if(flagSave);  fitswrite(I0f,[fn_prefix,'data_Inorm_broadband.fits']); end

%% Calculate and plot contrast in the lenslets

rawcontrast = [];
thput = [];

for fiberloc = [[6.1888; 0], [-3.0944; 5.3597], [-3.0944; -5.3597]]%[[5.3405; 0], [-2.6702; 4.625], [-2.6702; -4.625]]%[[6; 0], [-3; 5.1962], [-3; -5.1962]]
    for wi = 1:mp.Nsbp
        modvar.wpsbpIndex = 1;
        modvar.whichSource = 'star';
        modvar.sbpIndex = wi;
        modvar.ttIndex = 1;
        [Eout, Efiber] = model_full(mp, modvar);
        Istarmax = abs(Efiber(int8(fiberloc(2)*mp.Fend.res+length(Efiber)/2+1), int8(fiberloc(1)*mp.Fend.res+length(Efiber)/2+1))).^2;%max(max(abs(Efiber).^2));
        
        modvar.x_offset = fiberloc(1);
        modvar.y_offset = fiberloc(2);
        modvar.whichSource = 'offaxis';
        [Eout, Efiber] = model_full(mp, modvar);
        Iplanetmax = max(max(abs(Efiber).^2));
        rawcontrast = [rawcontrast Istarmax./Iplanetmax];
        thput = [thput sum(sum(abs(Efiber).^2))/mp.sumPupil];
    end
end

rawcontrast = squeeze(rawcontrast);
thput = squeeze(thput);

%% Plot raw contrast in each fiber
figure;
hold on;
% plot(mp.sbp_centers*1e9, log10(rawcontrast));
for i=1:mp.Fend.Nfiber
    plot(mp.sbp_centers*1e9, log10(rawcontrast(1+(i-1)*mp.Nsbp:i*mp.Nsbp)));
end
plot(586.5*ones(7), -12:-6, '--k');
plot(793.5*ones(7), -12:-6, '--k');
plot([550:830], -8*ones(281), '--k');

hold off;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'white');
xlabel('$\lambda$ (nm)', 'FontSize', 18, 'Interpreter', 'Latex');
ylabel('Contrast', 'FontSize', 18, 'Interpreter', 'Latex');
xlim([550 830]);
ylim([-12 -6]);
yticks(-12:1:-6);
legend({'(6.2, 0) $\lambda_0/D$', '(-3.0, 5.4) $\lambda_0/D$', '(-3.0, -5.4) $\lambda_0/D$'}, 'Location', 'North', 'Interpreter', 'Latex');
legend('boxoff');

%% Plot throughput in each fiber
figure;
hold on;
for i=1:mp.Fend.Nfiber
    plot(mp.sbp_centers*1e9, thput(1+(i-1)*mp.Nsbp:i*mp.Nsbp));
end

hold off;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'white');
xlabel('$\lambda$ (nm)', 'FontSize', 18, 'Interpreter', 'Latex');
ylabel('Throughput', 'FontSize', 18, 'Interpreter', 'Latex');
xlim([550 830]);
ylim([0 0.6]);