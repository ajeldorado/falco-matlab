% parse_compare2jobs.m
%
% read FALCO output data files and parse data
%
% Jorge Llop - Nov 15, 2018

clear all;close all;

addpath(genpath('/Users/jllopsay/Documents/MATLAB/WFIRST/parseFALCOData/utils'));

survey_str = 'DisOne';
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
%     num_str = nm(len-16:len-14);
%     wightNewCF = str2num(num_str);
    load([path,nm])
%     wightNewCF_arr(II) = wightNewCF;
    RCHist = out.InormHist;
    RC_arr(II) = InormHist(end);
    thput_arr(II) = thput;
    dm1rms_arr(II) = out.dm1.Srms(end);
    dm2rms_arr(II) = out.dm2.Srms(end);

    plot(1:out.Nitr,log10(RCHist(1:out.Nitr)),'LineWidth',3)
    hold on
end
xlabel('Iteration')
ylabel('log10(Normalized Intensity)')
lgd = legend('Omega = 5', 'Omega = 0', 'Omega = 5, First in itr 2, Beta -1','Omega = 5, First in itr 2',...
    'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
% lgd.Location = 'southeast';
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
hold off
figure(fig0)
% export_fig('/Users/jllopsay/Documents/GitHub/falco-matlab/parseData/RC_Mar07.png','-r300');

% fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
% for II=1:numfiles
%     nm = list_files(II).name;
%     len = numel(nm);
% %     num_str = nm(len-16:len-14);
% %     wightNewCF = str2num(num_str);
%     load([path,nm])
% %     wightNewCF_arr(II) = wightNewCF;
%     betaHist = out.log10regHist;
% 
%     plot(1:out.Nitr,betaHist,'LineWidth',3)
%     hold on
% end
% xlabel('Iteration')
% ylabel('Beta')
% lgd = legend('Omega = 5', 'Omega = 0', 'Omega = 5, First in itr 2, Beta -1','Omega = 5, First in itr 2',...
%     'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
% lgd.Location = 'southeast';
% % lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% % title(lgd,'Weight of new CF')
% set(gca,'FontSize',15)
% figure(fig0)
% hold off
% export_fig('/Users/jllopsay/Documents/GitHub/falco-matlab/parseData/Beta_Mar07.png','-r300');

fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
    len = numel(nm);
%     num_str = nm(len-16:len-14);
%     wightNewCF = str2num(num_str);
    load([path,nm])
%     wightNewCF_arr(II) = wightNewCF;
    thputHist = out.thput;

    plot(1:out.Nitr,thputHist(1:out.Nitr),'LineWidth',3)
    hold on
end
xlabel('Iteration')
ylabel('Throughput')
lgd = legend('Omega = 5', 'Omega = 0', 'Omega = 5, First in itr 2, Beta -1','Omega = 5, First in itr 2',...
    'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
% lgd.Location = 'southeast';
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
figure(fig0)
hold off
% export_fig('/Users/jllopsay/Documents/GitHub/falco-matlab/parseData/Thput_Mar07.png','-r300');

%%
fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
%     len = numel(nm);
%     num_str = nm(len-16:len-14);
%     wightNewCF = str2num(num_str);
    load([path,nm])
%     wightNewCF_arr(II) = wightNewCF;
    thputHist = out.thput;
    RCHist = out.InormHist;
    
    plot(log10(RCHist(1:out.Nitr)),thputHist(1:out.Nitr),'LineWidth',3)
    hold on
end
xlabel('RC')
ylabel('Throughput')
lgd = legend('Omega = 5', 'Omega = 0', 'New Omega','Omega = 5, First in itr 2',...
    'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
% lgd = legend('Omega = 5', 'Omega = 0', 'Omega = 5, First in itr 2, Beta -1','Omega = 5, First in itr 2',...
%     'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
lgd.Location = 'southeast';
% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
figure(fig0)
hold off
export_fig('/Users/jllopsay/Documents/GitHub/falco-matlab/parseData/ThputVsRC_Mar14.png','-r300');

%%
fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
%     len = numel(nm);
%     num_str = nm(len-16:len-14);
%     wightNewCF = str2num(num_str);
    load([path,nm])
%     wightNewCF_arr(II) = wightNewCF;
    thputHist = out.thput;
    SrmsHist = out.dm1.Srms;
    Srms2Hist = out.dm2.Srms;
    
    plot(1:out.Nitr,SrmsHist*1e9,'LineWidth',3)
    hold on
%     plot(1:out.Nitr,Srms2Hist,'LineWidth',3)
%     hold on
end
xlabel('it')
ylabel('Srms DM1 [nmRMS]')
ylim([10 50])
title('DM1 Stroke Height nmRMS')
lgd = legend('Omega = 5', 'Omega = 0', 'Omega = 5, First in itr 2, Beta -1','Omega = 5, First in itr 2',...
    'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
lgd.Location = 'southeast';
set(gca,'FontSize',15)
figure(fig0)
hold off
% export_fig('/Users/jllopsay/Documents/GitHub/falco-matlab/parseData/SrmsDM1_Mar07.png','-r300');


fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
for II=1:numfiles
    nm = list_files(II).name;
%     len = numel(nm);
%     num_str = nm(len-16:len-14);
%     wightNewCF = str2num(num_str);
    load([path,nm])
%     wightNewCF_arr(II) = wightNewCF;
    thputHist = out.thput;
    SrmsHist = out.dm1.Srms;
    Srms2Hist = out.dm2.Srms;
    
    plot(1:out.Nitr,Srms2Hist*1e9,'LineWidth',3)
    hold on
%     plot(1:out.Nitr,Srms2Hist,'LineWidth',3)
%     hold on
end
xlabel('it')
ylabel('Srms DM2 [nmRMS]')
ylim([10 50])
title('DM2 Stroke Height nmRMS')
lgd = legend('Omega = 5', 'Omega = 0', 'Omega = 5, First in itr 2, Beta -1','Omega = 5, First in itr 2',...
    'Omega = 5, First in itr 5','New Grid Search NI/thput','New Grid Search NI/thput^2');
% lgd.Location = 'southeast';

% lgd = legend('0.1','0.25','0.3','0.35','0.5','0.75','0.9','0.95','1.0');%,'Location','northeastoutside');
% title(lgd,'Weight of new CF')
set(gca,'FontSize',15)
figure(fig0)
hold off
% export_fig('/Users/jllopsay/Documents/GitHub/falco-matlab/parseData/SrmsDM2_Mar07.png','-r300');

