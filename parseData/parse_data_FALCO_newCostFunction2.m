% parse_compare2jobs.m
%
% read FALCO output data files and parse data
%
% Jorge Llop - Nov 15, 2018

clear all;close all;

addpath(genpath('/Users/jllopsay/Documents/MATLAB/WFIRST/parseFALCOData/utils'));

survey_str = 'testNewFALCO';
units = 'nm';
fontsz=16;

path = '/Users/jllopsay/Documents/GitHub/falco-matlab/data/ws/';
list_files = dir([path,'*',survey_str,'*all.mat']);
% baselinestrategy = 'Series0001_Trial0001_HLC_WFIRST180718_3DM48_z1_IWA2.7_OWA10_6lams575nm_BW10_plannedEFC_ResolutionPupil120_all.mat';
% load([path,baselinestrategy]);
% ni_os = InormHist(end);
% thput_os = thput;
% dm1rms_os = out.dm1.Srms(end);
% dm2rms_os = out.dm2.Srms(end);

numfiles = numel(list_files);
filename = 'gif_FinalImage';
fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
    len = numel(nm);
    num_str = nm(len-16:len-14);
    wightNewCF = str2num(num_str);
    load([path,nm])
    wightNewCF_arr(II) = wightNewCF;
    RCHist = out.InormHist;
    RC_arr(II) = InormHist(end);
    thput_arr(II) = thput;
    dm1rms_arr(II) = out.dm1.Srms(end);
    dm2rms_arr(II) = out.dm2.Srms(end);

    plot(1:out.Nitr+1,log10(RCHist),'LineWidth',3)
    hold on
end
xlabel('Iteration')
ylabel('log10(Normalized Intensity)')
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
hold off
figure(fig0)
% export_fig('RC_Feb22.png','-r300');

fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
    len = numel(nm);
    num_str = nm(len-16:len-14);
    wightNewCF = str2num(num_str);
    load([path,nm])
    wightNewCF_arr(II) = wightNewCF;
    betaHist = out.log10regHist;

    plot(1:out.Nitr,betaHist,'LineWidth',3)
    hold on
end
xlabel('Iteration')
ylabel('Beta')
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
figure(fig0)
hold off
% export_fig('Beta_Feb22.png','-r300');

fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
    len = numel(nm);
    num_str = nm(len-16:len-14);
    wightNewCF = str2num(num_str);
    load([path,nm])
    wightNewCF_arr(II) = wightNewCF;
    thputHist = out.thput;

    plot(1:out.Nitr+1,thputHist,'LineWidth',3)
    hold on
end
xlabel('Iteration')
ylabel('Throughput')
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
figure(fig0)
hold off
% export_fig('Thput_Feb22.png','-r300');

%%
fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
    len = numel(nm);
    num_str = nm(len-16:len-14);
    wightNewCF = str2num(num_str);
    load([path,nm])
    wightNewCF_arr(II) = wightNewCF;
    thputHist = out.thput;
    RCHist = out.InormHist;
    
    plot(log10(RCHist),thputHist,'LineWidth',3)
    hold on
end
xlabel('RC')
ylabel('Throughput')
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
figure(fig0)
hold off
% export_fig('ThputVsRC_Feb22.png','-r300');

% %% Get rid of bad elements
% % wightNewCF_arr(4) = [];
% % RC_arr(4) = [];
% % thput_arr(4) = [];
% % dm1rms_arr(4) = [];
% % dm2rms_arr(4) = [];
% %%
% figure(1)
% scatter(wightNewCF_arr(3:end),log10(RC_arr(3:end)))
% title('NI WFIRST BB planned EFC - iteration 270');
% % title(titl,'FontSize', fontsz);
% xlabel([survey_str,' [',units,']'])
% ylabel('NI (log scale)')
% set(gca,'FontSize',fontsz);
% export_fig('NIvswightNewCF_Nov19.png','-r300');
% figure(2)
% scatter(wightNewCF_arr(3:end),thput_arr(3:end))
% title('Thput WFIRST BB planned EFC - iteration 270');
% % title(titl,'FontSize', fontsz);
% xlabel([survey_str,' [',units,']'])
% ylabel('Thput')
% set(gca,'FontSize',fontsz);
% % export_fig('ThputvswightNewCF_Nov19.png','-r300');
% %%
% figure(3)
% scatter(wightNewCF_arr,dm1rms_arr)
% hold on
% scatter(wightNewCF_arr,dm2rms_arr)
% hold off
% title('DMs RMS WFIRST BB planned EFC - iteration 270');
% % title(titl,'FontSize', fontsz);
% xlabel([survey_str,' [',units,']'])
% ylabel('DMs nm RMS')
% set(gca,'FontSize',fontsz);
% legend('DM1' ,'DM2')
% % export_fig('DMsRMSvswightNewCF_Nov19.png','-r300');
% 
% %% Low order aberrations sensitivity
% %--Library locations
% mp.path.falco = '~/Documents/MATLAB/WFIRST/FALCOatHPCCaltech/falco-matlab/';  %--Location of FALCO
% mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
% 
% %%--Add to the MATLAB Path
% addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
% addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path
% 
% 
% %   ***   THE ONLY VARIABLES YOU NEED TO CHANGE ARE IN THIS CELL    ***
% %--Specify zernike modes and their RMS values to use.
% indsZnoll = 2:6; %--Noll indices of Zernikes to compute values for
% %--Annuli to compute sensitivities over. 
% % Columns are [inner radius, outer radius]. One row per annulus.
% Rsens = [3, 4;...
%          4, 8];
% list_files_config = dir([path,'*wightNewCF*config.mat']);
% list_files_snippet = dir([path,'*wightNewCF*snippet.mat']);
% 
% for II=1:numfiles
% 
%     %--Input files from the FALCO trial (the configuration data and final result snippet)
% %     fn_config = '~/Repos/falco-matlab/data/brief/Series0029_Trial0004_HLC_WFIRST180718_3DM48_z1_IWA2.7_OWA10_5lams575nm_BW10_plannedEFC_config.mat';
%     fn_config = [path,list_files_config(II).name];%
%     fn_snippet = [path,list_files_snippet(II).name];%
%     % Compute sensitivities to 1nm RMS of the specified Zernikes
%     dE2_array = func_eval_Zernike_sensitivities(indsZnoll,Rsens,fn_config,fn_snippet);
% end