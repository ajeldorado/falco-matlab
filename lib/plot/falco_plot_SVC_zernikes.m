% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to plot Zernike mode analysis.


% REVISION HISTORY:
% --------------
% Modified on 2022-04-06 by Niyati Desai.

% load frenchzernikes.mat %zernikes across varying bandwidths
% xaxis = bws;

load frenchRMSzernikes.mat
vals = vals(:,1:6);
xaxis = RMSs(1:6);

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



