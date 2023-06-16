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
% Created on 2022-04-06 by Niyati Desai.
% for Zernike plotting.
% Added contrast curve plotting 2022-06-24
% Added dimple contrast plotting 2022-10-16

%Data saved in niyatid/falco-matlab/2022SimulationData

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
%%
figure
load sawtoothdualzonecontrasts.mat
dualzonevals = vals;
bws(1) = 0;
xaxis = bws;

semilogy(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
set(gca, 'YScale', 'log')
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for Charge 6 Sawtooth Dualzone SVC')


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
set(gcf,'color','w');
ax.FontSize = 15;
ax.LineWidth = 1;
ylim([5e-13 5e-6]);
grid on;
box on;
% hold on
% plot(xaxis,frenchvals,'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880])
% plot(bws,mcmcvals,'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880])
hold on
plot(xaxis,classicvals,'LineWidth',2,'Color',[0.5 0 0.8])
hold on
plot(xaxis(2:end),stairvals(2:end),'LineWidth',2,'LineStyle',':','Color',[0 0.5 0.8])
hold on
plot(xaxis,vortexvals,'LineWidth',2,'LineStyle','--','Color',[0.9290, 0.6940, 0.1250])
hold on
plot(xaxis,dualzonevals,'LineWidth',2,'LineStyle','--','Color',[0.6940, 0.1250, 0.9290])
ax.XAxis.TickValues = [0 0.05 0.1 0.15 0.2];
xticklabels({'0','5%','10%','15%','20%'})

% legend('Wrapped','Sawtooth','Staircase','Classic Vortex');
% legend('Galicher','Sawtooth','Staircase','Classic Vortex');
legend('Sawtooth','Staircase','Classic Vortex','Dual Zone');
legend('Location','southeast')
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Contrasts for Charge 6 SVCs')
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



