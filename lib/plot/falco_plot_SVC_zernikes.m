% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to quickly plot SVC data.
%--Originally created to plot Zernike mode analysis.
%--Also used to generate contrast curves from saved mat files.


% REVISION HISTORY:
% --------------
% Modified on 2022-04-06 by Niyati Desai.

% load vortex8zernikes.mat %zernikes across varying bandwidths

% load frenchwrapped8zernikes.mat
% load staircase8zernikes.mat
% load sawtooth8zernikes.mat

load staircase6zernikes.mat
% load vortex6zernikes.mat
% load sawtooth6zernikes.mat
% load wrapped6zernikes.mat


%% Load Charge 6 Contrast Data

load staircasecontrasts.mat
bws(1) = 0;
% xaxis = DFTres(1:length(vals));

xaxis = bws;
stairvals = vals;
figure
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;

plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for Charge 6 Staircase SVC')
set(gca, 'YScale', 'log')


figure
load mcmccontrasts.mat
bws(1) = 0;
xaxis = bws;
mcmcvals = vals;
plot(xaxis,vals,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2)
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for Charge 6 MCMC SVC')
set(gca, 'YScale', 'log')

figure
load classicalcontrasts.mat
bws(1) = 0;
xaxis = bws;
classicvals = vals;
plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for Charge 6 Classically Wrapped SVC')
set(gca, 'YScale', 'log')

figure
load vortexcontrasts.mat
vortexvals = vals;
bws(1) = 0;
xaxis = bws;
plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for Charge 6 Vortex SVC')
set(gca, 'YScale', 'log')

%% Load Charge 8 Contrast Data

figure
load frenchcontrasts.mat
bws(1) = 0;
frenchvals = vals;
xaxis = bws;
plot(xaxis,vals,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2)
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for French Wrapped SVC')
set(gca, 'YScale', 'log')

load classical8contrasts.mat
classic8vals = vals;
load vortex8contrasts.mat
vortex8vals = vals;
load staircase8contrasts.mat
stair8vals = vals;


%% Contrast subtraction Figure
figure
hold on
xaxis = bws;
mcmcvort = mcmcvals-vortexvals;
classvort = classicvals-vortexvals;
stairvort = stairvals - vortexvals;

plot(bws,mcmcvort,'LineWidth',2)
hold on
plot(xaxis,classvort,'LineWidth',2)
hold on
plot(xaxis,stairvort,'LineWidth',2)
legend('MCMC-Vortex','Classic-Vortex', 'Staircase-Vortex');
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Contrast Differences: Charge 6 SVCs')
% set(gca, 'YScale', 'log')

%% Plot Charge 6&8 contrasts

figure
hold on
xaxis = bws;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
ylim([5e-13 5e-6]);
% hold on
plot(xaxis,frenchvals,'LineWidth',3,'Color',[0.4660, 0.6740, 0.1880])
% plot(bws,mcmcvals,'LineWidth',3,'Color',[0.4660, 0.6740, 0.1880])
hold on
plot(xaxis,classic8vals,'LineWidth',3,'Color',[0.5 0 0.8])
hold on
plot(xaxis(2:end),stair8vals(2:end),'LineWidth',3,'Color',[0 0.5 0.8])
hold on
plot(xaxis,vortex8vals,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
ax.XAxis.TickValues = [0 0.05 0.1 0.15 0.2];
xticklabels({'0','5%','10%','15%','20%'})

% legend('Wrapped','Sawtooth','Staircase','Classic Vortex');
legend('Galicher','Sawtooth','Staircase','Classic Vortex');
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Contrasts for Charge 8 SVCs')
set(gca, 'YScale', 'log')

% figure
% ax = gca;
% ax.FontSize = 12;
% ax.LineWidth = 2;
% % x.YScale = 'log';
% % hold on
% 
% Tminutes = T./60;
% 
% plot(xaxis,Tminutes,'Color',[0.5 0 0.8],'LineWidth',2)
% % hold on
% % 
% % load staircase6.mat
% % xaxis = res;
% % plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
% % legend('no WFSC','10 iterations of EFC')
% xlabel('Model Resolution [pix per \lambda/D]');
% ylabel('Computation Time (minutes)');
% title('Computing Time no WFSC')
% % set(gca, 'YScale', 'log')




%% ZERNIKE PLOTTING
% 
% vals = vals(:,1:6);
% xaxis = RMSs(1:6);
% vals = vals';
xaxis = RMSs';

% CONVERT RMS to wavs??? (need wavelength)
xaxis = xaxis./550; %for 550nm wavelength

tiptilt = mean(cat(1,vals(1,:),vals(2,:)),1);
defocus = vals(3,:);
astig = mean(cat(1,vals(4,:),vals(5,:)),1);
coma = mean(cat(1,vals(6,:),vals(7,:)),1);


% load chromaticvortexzernikes.mat
% tiptilt2 = mean(cat(1,vals(1,:),vals(2,:)),1)
% focus2 = vals(3,:)
% astig2 = mean(cat(1,vals(4,:),vals(5,:)),1)
% coma2 = mean(cat(1,vals(6,:),vals(7,:)),1)

figure
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
x.YScale = 'log';

hold on
plot(xaxis,tiptilt,'Color',[0.5 0 0.8],'LineWidth',2)
plot(xaxis,defocus,'Color',[0 0 0.8],'LineWidth',2)
plot(xaxis,astig,'Color',[0.5 0 0],'LineWidth',2)
plot(xaxis,coma,'Color',[0 0.7 0],'LineWidth',2)
% plot(xaxis,tiptilt2,'--','Color',[0.5 0 0.8],'LineWidth',2)
% plot(xaxis,focus2,'--','Color',[0 0 0.8],'LineWidth',2)
% plot(xaxis,astig2,'--','Color',[0.5 0 0],'LineWidth',2)
% plot(xaxis,coma2,'--','Color',[0 0.7 0],'LineWidth',2)

% ylim([5e-12 5e-9]);
% xlim([1e-5 1]);
legend({'tip-tilt','defocus','astig','coma'},'Location','northwest')
xlabel('RMS (wavs)');
ylabel('|dE|^2 at 550nm with for 2-4 l/D');
title('Zernike Sensitivities for Staircase Charge 6 SVC')
set(gca, 'YScale', 'log', 'XScale','log')
% ax.XAxis.Exponent = -9;



