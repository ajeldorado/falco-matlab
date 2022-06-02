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

% load frenchzernikes.mat %zernikes across varying bandwidths

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
plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
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

figure
load frenchcontrasts.mat
bws(1) = 0;
frenchvals = vals;
xaxis = bws;
plot(xaxis,vals,'Color',[0 0.5 0.8],'LineWidth',2)
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Raw Contrasts for French Wrapped SVC')
set(gca, 'YScale', 'log')


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


figure
hold on
xaxis = bws;

plot(bws,mcmcvals,'LineWidth',2)
hold on
plot(xaxis,classicvals,'LineWidth',2)
hold on
plot(xaxis,stairvals,'LineWidth',2)
hold on
plot(xaxis,vortexvals,'LineWidth',2)
hold on
plot(xaxis,frenchvals,'LineWidth',2)
legend('MCMC','Classic','Staircase','Vortex','French');
xlabel('Bandwidth');
ylabel('Raw Contrast');
title('Contrasts for various SVCs')
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




%%

vals = vals(:,1:6);
% xaxis = RMSs(1:6);

tiptilt = mean(cat(1,vals(1,:),vals(2,:)),1)
focus = vals(3,:)
astig = mean(cat(1,vals(4,:),vals(5,:)),1)
coma = mean(cat(1,vals(6,:),vals(7,:)),1)


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
plot(xaxis,focus,'Color',[0 0 0.8],'LineWidth',2)
plot(xaxis,astig,'Color',[0.5 0 0],'LineWidth',2)
plot(xaxis,coma,'Color',[0 0.7 0],'LineWidth',2)
% plot(xaxis,tiptilt2,'--','Color',[0.5 0 0.8],'LineWidth',2)
% plot(xaxis,focus2,'--','Color',[0 0 0.8],'LineWidth',2)
% plot(xaxis,astig2,'--','Color',[0.5 0 0],'LineWidth',2)
% plot(xaxis,coma2,'--','Color',[0 0.7 0],'LineWidth',2)
legend({'tiptilt','focus','astig','coma'},'Location','northwest')
xlabel('RMS (nanometers)');
ylabel('|dE|^2 at 550nm with for 2-4 l/D');
title('Zernike Sensitivities for French Wrapped Charge 8 SVC')
set(gca, 'YScale', 'log', 'XScale','log')
% ax.XAxis.Exponent = -9;



