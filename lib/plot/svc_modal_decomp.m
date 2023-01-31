% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to plot phase profiles for various svc designs
%--and perform modal decompositions.


% REVISION HISTORY:
% --------------
% Created on 2022-04-11 by Niyati Desai.

%% Demo
%%Time specifications:
Fs = 100;                       % samples per second
dt = 1/Fs;                     % seconds per sample
StopTime = 1;                  % seconds
t = linspace(-pi,pi,Fs);%(0:dt:StopTime-dt)';
N = size(t,2);
%%Sine wave:
Fc = 12;                       % hertz
x = cos(Fc*t);
%%Fourier Transform:
X = fftshift(fft(x));
%%Frequency specifications:
dF = Fs/N;                      % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz
%%Plot the spectrum:
figure(311);
plot(f,abs(X)/N);
xlabel('Frequency (in hertz)');
title('Magnitude Response');

%% Set Parameters
clear all
L = 1E7;

N = 2^nextpow2(L)
THETA = linspace(-pi,pi,N);
prof = 0.*THETA;
scale = 1;
linVec = {'--',':','-'};

lambda = 1.1;

if lambda < 1
    i = 1;
elseif lambda > 1
    i = 2;
else
    i = 3;
end

%%Frequency specifications:
Fs = N;
dF = Fs/N;                      % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz

%% Practice sinusoid
p=8;

figure(111)
phase = sin(p*THETA);
plot(THETA,phase,'Color',[0 0.5 0.8],'LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
xlim([-pi pi]);
axis on
ax.XAxis.TickValues = [-pi -pi/2 0 pi/2 pi ];
ax.YAxis.TickValues = [-8*pi -6*pi -4*pi -2*pi 0 2*pi 4*pi 6*pi 8*pi];
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticklabels({'-8\pi','-6\pi','-4\pi','-2\pi','0','2\pi','4\pi','6\pi','8\pi'})
xlabel('Theta'); 
ylabel('Phase');
title("Sinusoid")

% % 

%%Fourier Transform:
fftsin = abs(fftshift(fft(phase)))/N;

%%Plot the spectrum:
figure(311);
plot(f,fftsin);
title("FFT of Sinusoid")



%% Topologies

%cos

charge0 = 6;
type0 = "Cos";
z_m = besselzero(0,1);
z_m = z_m(end);
phase0 = z_m*cos(charge0*THETA);

% figure(1);
% plot(THETA,phase0,'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2)
% ax = gca;
% ax.FontSize = 20;
% ax.LineWidth = 3;
% xlim([-pi pi]);
% axis on
% ax.XAxis.TickValues = [-pi -pi/2 0 pi/2 pi ];
% ax.YAxis.TickValues = [-8*pi -6*pi -4*pi -2*pi 0 2*pi 4*pi 6*pi 8*pi];
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% yticklabels({'-8\pi','-6\pi','-4\pi','-2\pi','0','2\pi','4\pi','6\pi','8\pi'})
% xlabel('Theta'); 
% ylabel('Phase');
% % title("Charge "+charge0+" "+type0 +" Phase Profile")
% title(type0 +" Phase Profile")

t0 = exp(1j/lambda*phase0);

%%Fourier Transform:
myfftcos = abs(fftshift(fft(t0)))/N;

% Plot the spectrum:
% figure (310)
% plot(f,myfftcos);
% ylabel('|C_m|^2'); 
% xlabel('Mode');
% title("Modal Decomposition for Cos SVC")
% set(gca, 'YScale', 'log')


% vortex
charge1 = 6;
phase1 = charge1*THETA;
type1 = "Classic Vortex";


% figure(N)
% % subplot(2,2,2)
% figure(1);
% plot(THETA,phase1,'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2)
% ax = gca;
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% ax.FontSize = 28;
% ax.LineWidth = 1;
% xlim([-pi pi]);
% axis on
% set(gca,'TickDir','out');
% set(gcf,'color','w');
% ax.XAxis.TickValues = [-pi -pi/2 0 pi/2 pi ];
% ax.YAxis.TickValues = [-8*pi -6*pi -4*pi -2*pi 0 2*pi 4*pi 6*pi 8*pi];
% xticklabels({'$-\pi$','$-\frac{\pi}{2}$','$0$','$\frac{\pi}{2}$','$\pi$'})
% yticklabels({'$-8\pi$','$-6\pi$','$-4\pi$','$-2\pi$','$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$'})
% xlabel('$\Theta$ (azimuthal angle)'); 
% ylabel('$\phi$ (phase)');
% title("Charge "+charge1+" "+type1 +" Phase Profile")
% title(type1 +" Phase Profile")

t1 = exp(1j/lambda*phase1);

%%Fourier Transform:
myfftvortex = abs(fftshift(fft(t1)))/N;

%Plot the spectrum:
% figure (311)
% plot(f,myfftvortex);
% ylabel('|C_m|^2'); 
% xlabel('Mode');
% title("Modal Decomposition for Charge 8 Vortex SVC")
% set(gca, 'YScale', 'log')


% sawtooth

charge2 = 6;
domain = (THETA >= 0) & (THETA <= pi);
phase2(domain) = charge2*rem(THETA(domain),2*pi./charge2);
domain = (THETA >= -pi) & (THETA < 0);
phase2(domain) = charge2*rem((THETA(domain)+pi),2*pi./charge2);

newphase = [phase2 phase2 phase2 phase2];
newtheta = linspace(-4*pi,4*pi,N);


type2 = "Sawtooth Vortex";
% subplot(2,2,3)
% figure(2);
% plot(THETA,phase2,'Color',[0.5 0 0.8],'LineWidth',2)
% ax = gca;
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% ax.FontSize = 28;
% ax.LineWidth = 1;
% xlim([-pi pi]);
% axis on
% set(gca,'TickDir','out');
% set(gcf,'color','w');
% ax.XAxis.TickValues = [-pi -pi/2 0 pi/2 pi ];
% ax.YAxis.TickValues = [0 pi 2*pi ];
% xticklabels({'$-\pi$','$-\frac{\pi}{2}$','$0$','$\frac{\pi}{2}$','$\pi$'})
% yticklabels({'$0$','$\pi$','$2\pi$'})
% xlabel('$\Theta$ (azimuthal angle)'); 
% ylabel('$\phi$ (phase)');
% title("Charge "+charge2+" "+type2 +" Phase Profile")
% title(type2 +" Phase Profile")



t2 = exp(1j/lambda*phase2);
%%Fourier Transform:
myfftsawtooth = abs(fftshift(fft(t2)))/N;

%%Plot the spectrum:
% figure(311);
% plot(f,myfftclassical);
% ylabel('|C_m|^2'); 
% xlabel('Mode');
% title("Modal Decomposition for Charge 8 Classically Wrapped SVC")
% set(gca, 'YScale', 'log')


% french wrapped

domain = (THETA > 0) & (THETA < 3*pi/8);
phase3(domain) = 8*THETA(domain);
domain = (THETA > 3*pi/8) & (THETA < pi/2);
phase3(domain) = 8*THETA(domain)-2*pi;
domain = (THETA > pi/2) & (THETA < 5*pi/8);
phase3(domain) = 8*THETA(domain)-4*pi;
domain = (THETA > 5*pi/8) & (THETA <= pi);
phase3(domain) = 8*THETA(domain)-6*pi;

domain = (THETA > -pi) & (THETA < -5*pi/8);
phase3(domain) = 8*(THETA(domain)+pi);
domain = (THETA > -5*pi/8) & (THETA < -pi/2);
phase3(domain) = 8*(THETA(domain)+pi)-2*pi;
domain = (THETA > -pi/2) & (THETA < -3*pi/8);
phase3(domain) = 8*(THETA(domain)+pi)-4*pi;
domain = (THETA > -3*pi/8) & (THETA < 0);
phase3(domain) = 8*(THETA(domain)+pi)-6*pi;

type3 = "Galicher Vortex";
% subplot(2,2,4)
% figure(3);
% plot(THETA,phase3,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2)
% ax = gca;
% ax.FontSize = 20;
% ax.LineWidth = 3;
% xlim([-pi pi]);
% axis on
% ax.XAxis.TickValues = [-pi -5*pi/8 -pi/2 -3*pi/8 0 3*pi/8 pi/2 5*pi/8 pi ];
% ax.YAxis.TickValues = [-pi 0 pi 2*pi 3*pi];
% yticklabels({'-\pi','0','\pi','2\pi','3\pi'})
% xticklabels({'-\pi','^{-5\pi}/_{8}','^{-\pi}/_{2}' ,'^{-3\pi}/_{8}','0','^{3\pi}/_{8}','^{\pi}/_{2}','^{5\pi}/_{8}','\pi'})
% xlabel('Theta'); 
% ylabel('Phase');
% % title("Charge 8 "+type3 +" Phase Profile")
% title(type3 +" Phase Profile")

t3 = exp(1j/lambda*phase3);

%%Fourier Transform:
myfftfrench = abs(fftshift(fft(t3)))/N;

%%Plot the spectrum:
% figure(311)
% plot(f,myfftfrench);
% set(gca, 'YScale', 'log')
% ylabel('|C_m|^2'); 
% xlabel('Mode');
% title("Modal Decomposition for Charge 8 French Wrapped SVC")

% mcmc6 wrapped

domain = (THETA > 0) & (THETA < 1.84799568);
phase4(domain) = 6*THETA(domain);
domain = (THETA > 1.84799568) & (THETA < 2.52559409);
phase4(domain) = 6*THETA(domain) - 2*pi;
domain = (THETA > 2.52559409) & (THETA <= pi);
phase4(domain) = 6*THETA(domain) - 4*pi;
domain = (THETA > -pi+2.52559409) & (THETA < 0);
phase4(domain) = 6*(THETA(domain)+pi)-4*pi;
domain = (THETA > -pi+1.84799568) & (THETA < -pi+2.52559409);
phase4(domain) = 6*(THETA(domain)+pi)-2*pi;
domain = (THETA >= -pi) & (THETA < -pi+1.84799568);
phase4(domain) = 6*(THETA(domain)+pi);


type4 = "Wrapped Vortex";
% subplot(2,2,4)
% figure(4);
% plot(THETA,phase4,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2)
% ax = gca;
% ax.FontSize = 20;
% ax.LineWidth = 3;
% xlim([-pi pi]);
% axis on
% ax.XAxis.TickValues = [-pi -2*pi/5 -1*pi/5 0 3*pi/5 4*pi/5 pi ];
% ax.YAxis.TickValues = [-pi 0 pi 2*pi 3*pi];
% yticklabels({'-\pi','0','\pi','2\pi','3\pi'})
% xticklabels({'-\pi','^{-2\pi}/_{5}','^{-\pi}/_{5}','0','^{3\pi}/_{5}','^{4\pi}/_{5}','\pi'})
% xlabel('Theta'); 
% ylabel('Phase');
% % title("Charge 6 "+type4 +" Phase Profile")
% title(type4 +" Phase Profile")

t4 = exp(1j/lambda*phase4);

%%Fourier Transform:
myfftmcmc = abs(fftshift(fft(t4)))/N;

%%Plot the spectrum:
% figure(311)
% plot(f,myfftmcmc);
% set(gca, 'YScale', 'log')
% ylabel('|C_m|^2'); 
% xlabel('Mode');
% title("Modal Decomposition for Charge 6 MCMC Wrapped SVC")


% staircase
Nsteps = 6;
charge5 = 6;
phase5 = floor(mod((THETA+pi)/(2*pi)*charge5, 1)*Nsteps)/Nsteps*2*pi;


type5 = "Staircase Vortex";
% subplot(2,2,2)
% figure(5);
% plot(THETA,phase5,'Color',[0 0.5 0.8],'LineWidth',2)
% ax = gca;
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% ax.FontSize = 28;
% ax.LineWidth = 1;
% xlim([-pi pi]);
% ylim([0 2*pi*1.03]);
% axis on
% set(gca,'TickDir','out');
% set(gcf,'color','w');
% ax.XAxis.TickValues = linspace(-pi,pi,7);
% ax.YAxis.TickValues = [0 pi 2*pi ];
% xticklabels({'$-\pi$','$-\frac{2\pi}{3}$','$-\frac{\pi}{3}$','$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$'})
% yticklabels({'$0$','$\pi$','$2\pi$'})
% xlabel('$\Theta$ (azimuthal angle)'); 
% ylabel('$\phi$ (phase)');
% title("Charge 6 "+type5 +" Phase Profile")
% title(type5 +" Phase Profile")
 

t5 = exp(1j/lambda*phase5);

%%Fourier Transform:
myfftstaircase = abs(fftshift(fft(t5)))/N;

%%Plot the spectrum:
% figure(311)
% plot(f,myfftstaircase);
% set(gca, 'YScale', 'log')
% ylabel('|C_m|^2'); 
% xlabel('Mode');
% title("Modal Decomposition for Charge 6 Staircase SVC")

%% Plot with stem markers 
% close all

figure(673)
epsilon = 0.05;
hold on
xdata = (0:1:N-1);
ax = gca;
set(gcf,'color','w')
set(gca,'YScale', 'log');
ax.FontSize = 20;
ax.LineWidth = 2;
stem(f,myfftvortex,'d','MarkerSize',10,'LineStyle','-','LineWidth',2,'Color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
stem(f+epsilon,myfftsawtooth,'s','MarkerSize',20,'LineStyle','--','LineWidth',2,'Color',[0.5 0 0.8]);
% stem(f,myfftfrench,'MarkerSize',10,'LineStyle','-','LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]);
% stem(f,myfftmcmc,'MarkerSize',10,'LineStyle','-','LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]);
stem(f-epsilon,myfftstaircase,'*','MarkerSize',15,'LineStyle','-.','LineWidth',2,'Color',[0 0.5 0.8]);
% stem(f,myfftcos,'.','MarkerSize',10,'LineStyle','-','LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]);
if lambda == 1
%     stem(f,myfftvortex,'.','LineStyle',linVec{i},'LineWidth',2,'Color','k');
end

% xlim([-10 23]);
xlim([-1.2 13])
ylim([5E-3 10])
xticks([0:2:12])

%type1-vortex,type2-sawtooth,type3-frenchwrapped,type4-mcmc,type5-staircase,
%type0-cos
legend(type1,type2,type5)
legend('Location','northwest')


% line_type = ['-',"--", ":"];
% line_type = ['-', ":"];

% line_names = ["design wavelength", "0.95 lambda factor", "1.05 lambda factor"];
% line_names = ["design wavelength", "1.05 lambda factor"];
% for j =1:length(line_type)
%     plot([NaN NaN], [NaN NaN],line_type(j), 'Color', 'k', 'DisplayName', line_names(j))
% end

hold off 
% ylabel('|C_m|^2'); 
ylabel('Power');
xlabel('Mode');
% title("Modal Decomposition for Charge 2 SVCs")
