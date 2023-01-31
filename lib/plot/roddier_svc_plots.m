% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script generate figures for Roddier SVC paper.
%--Used to generate chromaticity V curves from saved mat files.


% REVISION HISTORY:
% --------------
% Created on 2023-01-30 by Niyati Desai.


load vortex2jan26.mat
vortcontrasts = rawcontrasts;
temp = sum(rawcontrasts) - 0.5*(rawcontrasts(1)+rawcontrasts(end));
vortavg = temp/(length(rawcontrasts)-1)
facs1 = phasescalefacs;
load roddiersaw2jan26.mat
roddiercontrasts = rawcontrasts;
temp = sum(rawcontrasts) - 0.5*(rawcontrasts(1)+rawcontrasts(end));
roddieravg = temp/(length(rawcontrasts)-1)
facs2 = phasescalefacs;
load dzpmsawtoothjan26.mat
dzpmcontrasts = rawcontrasts;
temp = sum(rawcontrasts) - 0.5*(rawcontrasts(1)+rawcontrasts(end));
dzpmavg = temp/(length(rawcontrasts)-1)
facs3 = phasescalefacs;
load sawtoothjan26.mat
sawtoothcontrasts = rawcontrasts;
temp = sum(rawcontrasts) - 0.5*(rawcontrasts(1)+rawcontrasts(end));
sawtoothavg = temp/(length(rawcontrasts)-1)
facs4 = phasescalefacs;
%     load galicherjan26.mat
%     galichercontrasts = rawcontrasts;
%     facs5 = phasescalefacs;
load wrapped6jan26.mat
wrapped6contrasts = rawcontrasts;
temp = sum(rawcontrasts) - 0.5*(rawcontrasts(1)+rawcontrasts(end));
wrapped6avg = temp/(length(rawcontrasts)-1)
facs6 = phasescalefacs;



figure (11)
xaxis = phasescalefacs;
ax = gca;
ax.FontSize = 15;
ax.LineWidth = 1;
ylim([5e-12 1e-3]);
xlim([0.9 1.1]);
hold on
plot(facs1,vortcontrasts,'LineWidth',2,'Color',[0 0.5 0.8])
hold on
plot (facs4,sawtoothcontrasts,'LineWidth',2,'Color',[0.4 0.7 0.1])
hold on
plot (facs2,roddiercontrasts,'LineWidth',2,'Color',[0.5 0 0.8])
hold on
plot (facs3,dzpmcontrasts,'LineWidth',2,'Color',[0.1 0.9 0.8])
hold on
%     plot (facs5,galichercontrasts,'LineWidth',2,'Color',[0.9 0 0])
%     hold on
plot (facs6,wrapped6contrasts,'LineWidth',2,'Color',[0.9 0.9 0])
hold on
legend('Vortex','Sawtooth','Roddier+Sawtooth','DZPM+Sawtooth','Wrapped6')%,'Galicher')
legend('Location','northeast')
xlabel('Phase Scale Factor');
% % ylabel('Max of final focal plane \lambda/D');
ylabel('Average from 3-10 \lambda/D');
title('Chromaticity for Charge 6 FPMs')
set(gca, 'YScale', 'log')
set(gcf,'color','w');
grid on;