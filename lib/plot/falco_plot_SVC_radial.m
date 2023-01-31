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

%% Focal Plane Gif creator for varying radii

load radiusroddier.mat;

[x,y,z] = size(imcube);
for index=1:z
    figure();
    myim = imcube(:,:,index);
    imagesc(myim/2^16);
    axis tight; axis equal;
    title("Final Focal Plane Image for Radii = " +radii(index));
    set(gca,'tickdir','out')
    set(gcf,'Color','w');
    colorbar;
    caxis([5e-11,1e-8]);
    set(gca,'ColorScale','log')
    saveas(gcf,"radius"+ radii(index)+".png")
end



%% Null Depth comparison w/ vs w/o raddier/dual zone

load nocg1121.mat;
nocg = nulldepths;
load roddier1128.mat;
roddier = nulldepths;
load dualzone1128.mat;
dualzone = nulldepths;
load vortex1128.mat;
vortex = nulldepths;
load sawtooth1128.mat;
sawtooth = nulldepths;
load sawtoothroddier1128.mat;
sawtoothroddier = nulldepths;
load sawtoothdualzone1128.mat;
sawtoothdualzone = nulldepths;

figure
hold on
% xaxis = bws;
xaxis = phasescalefacs;
ax = gca;
ax.FontSize = 15;
ax.LineWidth = 1;
ylim([5e-12 1e-3]);
xlim([0.9 1.12]);
plot(xaxis,nocg, 'b-', 'LineWidth', 2);
hold on;
plot(xaxis,roddier, 'g-', 'LineWidth', 2);
hold on;
plot(xaxis,dualzone, 'r-', 'LineWidth', 2);
hold on;
plot(xaxis,vortex, 'k-', 'LineWidth', 2);
hold on;
plot(xaxis,sawtooth, 'm-', 'LineWidth', 2);
hold on;
plot(xaxis,sawtoothroddier, 'y-', 'LineWidth', 2);
hold on;
plot(xaxis,sawtoothdualzone, 'c-', 'LineWidth', 2);
hold on;


% % legend('Vortex','Sawtooth','Roddier Vortex','Dual Zone Vortex');

% legend('Location','southeast')
legend('No CG','Roddier','Dual Zone','Vortex','Sawtooth','Sawtooth+Roddier','Sawtooth+DualZone');
legend('Location','southeast')
hold on;

xlabel('Bandwidth');
ylabel('Average from 3-10 \lambda/D');
title('Chromaticity')
set(gca, 'YScale', 'log')
set(gcf,'color','w');
grid on;

%% COSINE CHROMATICITY FOR LORENZO

% Load Vortex Dimple Contrast Data

load cos2nulldepths.mat
cosvals = nulldepths;%vals;
load vortex2nulldepths.mat
vortvals = nulldepths;%vals;
load sawtooth2nulldepths.mat
sawtoothvals = nulldepths;%vals;
% load nodimplevort2.mat
% nodimplevortvals = nulldepths;%vals;
% bws(1) = 0;

figure
hold on
% xaxis = bws;
xaxis = phasescalefacs;
ax = gca;
ax.FontSize = 15;
ax.LineWidth = 1;
ylim([5e-12 1e-3]);
xlim([0.9 1.12]);


% plot(xaxis,nodimplevortvals,'LineWidth',2,'Color',[0 0.5 0.8])
% hold on
plot(xaxis,sawtoothvals,'LineWidth',2,'LineStyle','--','Color',[0.4660, 0.6740, 0.1880])
hold on
plot(xaxis,cosvals,'LineWidth',2,'LineStyle','-.','Color',[0.5 0 0.8])
hold on
plot(xaxis,vortvals,'LineWidth',2,'LineStyle',':','Color',[0.9290, 0.6940, 0.1250])
% ax.XAxis.TickValues = [0 0.05 0.1 0.15 0.2];
% xticklabels({'0','5%','10%','15%','20%'})

legend('Sawtooth','Cos','Vort')
legend('Location','southeast')
hold on;
xlabel('Phase Scale Factor');
ylabel('Max of final focal plane \lambda/D');
title('Chromaticity')
set(gca, 'YScale', 'log')
set(gcf,'color','w');
grid on;

%% Contrast Curves from Sims
% 
% load nocg1121.mat;
% load roddier1128.mat;
% load dualzone1128.mat;
% load vortex1128.mat;
% load sawtooth1128.mat;
% load sawtoothroddier1128.mat;
load sawtoothdualzone1128.mat;


% figure(203);
% imagesc(im0/2^16)
% axis image; 
% colorbar;
%

%% Comparing Roddier Radii contrast curves

myLineStyleOrders = {'o','+','*','x','s','d','v','p','h','^'};

newcolors = flip([0.494 0.184 0.556
         0.00 0.447 0.741
         0.301 0.745 0.933
         0.466 0.474 0.288
         0.25 0.80 0.54
         0.929 0.694 0.125
         0.85 0.525 0.098
         0.93 0.14 0.14
         0.80 0.25 0.47
         0.635 0.078 0.184],1);
     
load radiusroddier.mat
[x,y,z] = size(imcube);
for index=1:z

    myim = imcube(:,:,index);
    num_r = 100;
    [r_arr,mymeans] = profiler(myim,num_r,mp.Fend.res);
    
    figure(105)
    ax = gca;
    set(gcf,'color','w');
    thisColor = newcolors(index,:);
    thisMarkerStyle = myLineStyleOrders{index};
    plot(r_arr,mymeans,'-', 'Color', thisColor, ...
        'Marker', thisMarkerStyle, ...
        'LineWidth', 1, 'MarkerSize', 6);
    hold on;
end
xlabel('Angular Separation (\lambda/D)','fontsize',12)
ylabel('Raw Contrasts','fontsize',12)
set(gca, 'YScale', 'log','fontsize',12)
legendCell = cellstr(num2str(radii'));
legend(legendCell)

%% Plot the first profile
load vortex1128.mat
[r_arr,mymeans] = profiler(im3,100,mp.Fend.res);
figure(105)
plot(r_arr, mymeans, 'b-', 'LineWidth', 2);
ax = gca;
set(gcf,'color','w');
set(gcf,'position',[10,10,700,400])
grid on
title('Average Radial Profile', 'FontSize', 10);
xlabel('Angular Separation (\lambda/D)','fontsize',12)
ylabel('Raw Contrasts','fontsize',12)
set(gca, 'YScale', 'log','fontsize',12)
hold on
%% Add additional profiles overlayed on same plot
figure(108)
hold on
plot(r_arr, mymeans, 'c-', 'LineWidth', 2);
%% Assign each profile in legend
legend('No CG','Roddier','Dual Zone','Vortex','Sawtooth','Sawtooth+Roddier','Sawtooth+DualZone');
legend('Location','northeast')
hold on;
%%

%% Profiler

function [r_arr, mymeans] = profiler(im0, num_r,lamOverD) %im, 100, mp.Fend.res
%     im0 = im3;

%     lamOverD = mp.Fend.res;% = 10;
    N = length(im0); %128;
    [X,Y] = meshgrid(-N/2:N/2-1); 
    [THETA,RHO] = cart2pol(X,Y);

    num_r = 100;
    IWA = 0; % inner working angle
    OWA = 15; % outer working angle
    r_arr = linspace(IWA,OWA,num_r);
    mymeans = [];

    delta_r = 0.5; % lam/D %thickness of annulus
    for r=r_arr
        mask = zeros(N,N);

        r_in = r-delta_r; 
        r_out = r+delta_r;
        r_in_pix = r_in*lamOverD;
        r_out_pix = r_out*lamOverD;
    %     mask(RHO>r_in_pix && RHO<r_out_pix) = 1;
        mask = RHO > r_in_pix & RHO<r_out_pix;
        polarim = mask.*im0;
        newmean = mean(nonzeros(polarim));
        mymeans = [mymeans newmean];
    end
end