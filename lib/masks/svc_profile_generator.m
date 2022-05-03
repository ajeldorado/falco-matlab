% Copyright 2022, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Script to plot phase maps and profiles for various svc designs.


% REVISION HISTORY:
% --------------
% Modified on 2022-04-11 by Niyati Desai.

%% Demo
%%Time specifications:
Fs = 100;                      % samples per second
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
L = 10000;
N = 2000;
N = 2^nextpow2(L)
THETA = linspace(-pi,pi,N);
prof = 0.*THETA;
scale = 1;
linVec = {'--',':','-'};
lambda = 1.05;

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



%% vortex
charge1 = 6;
phase1 = charge1*THETA;
type1 = "Vortex";


figure(N)
subplot(2,2,2)
plot(THETA,phase1,'Color',[0 0.5 0.8],'LineWidth',2)
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
title("Charge "+charge1+" "+type1 +" Phase Profile")

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



% classical wrapped

charge2 = 6;
domain = (THETA >= 0) & (THETA <= pi);
phase2(domain) = charge2*rem(THETA(domain),2*pi./charge2);
domain = (THETA >= -pi) & (THETA < 0);
phase2(domain) = charge2*rem((THETA(domain)+pi),2*pi./charge2);

newphase = [phase2 phase2 phase2 phase2];
newtheta = linspace(-4*pi,4*pi,N);


type2 = "Classically Wrapped Vortex";
subplot(2,2,3)
plot(THETA,phase2,'Color',[0.5 0 0.8],'LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
xlim([-pi pi]);
axis on
ax.XAxis.TickValues = [-pi -pi/2 0 pi/2 pi ];
ax.YAxis.TickValues = [-2*pi 0 2*pi ];
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticklabels({'-2\pi','0','2\pi'})
xlabel('Theta'); 
ylabel('Phase');
title("Charge "+charge2+" "+type2 +" Phase Profile")



t2 = exp(1j/lambda*phase2);
%%Fourier Transform:
myfftclassical = abs(fftshift(fft(t2)))/N;

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

type3 = "French Wrapped Vortex";
subplot(2,2,4)
plot(THETA,phase3,'Color',[0.5 0.8 0],'LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
xlim([-pi pi]);
axis on
ax.XAxis.TickValues = [-pi -5*pi/8 -pi/2 -3*pi/8 0 3*pi/8 pi/2 5*pi/8 pi ];
ax.YAxis.TickValues = [-pi 0 pi 2*pi 3*pi];
yticklabels({'-\pi','0','\pi','2\pi','3\pi'})
xticklabels({'-\pi','-5\pi/8', '-\pi/2','-3\pi/8','0','3\pi/8','\pi/2','5\pi/8','\pi'})
xlabel('Theta'); 
ylabel('Phase');
title("Charge 8 "+type3 +" Phase Profile")

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


type4 = "MCMC Vortex";
subplot(2,2,4)
plot(THETA,phase4,'Color',[0.5 0.8 0],'LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
xlim([-pi pi]);
axis on
ax.XAxis.TickValues = [-pi -2*pi/5 -1*pi/5 0 3*pi/5 4*pi/5 pi ];
ax.YAxis.TickValues = [-pi 0 pi 2*pi 3*pi];
yticklabels({'-\pi','0','\pi','2\pi','3\pi'})
xticklabels({'-\pi','-2\pi/5','-1\pi/5','0','3\pi/5','4\pi/5','\pi'})
xlabel('Theta'); 
ylabel('Phase');
title("Charge 6 "+type4 +" Phase Profile")

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
subplot(2,2,2)
plot(THETA,phase5,'Color',[0 0.5 0.8],'LineWidth',2)
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
xlim([-pi pi]);
axis on
ax.XAxis.TickValues = linspace(-pi,pi,7);
ax.YAxis.TickValues = [-2*pi 0 2*pi ];
xticklabels({'-\pi','-2\pi/3','-\pi/3','0','\pi/3','2\pi/3','\pi'})
yticklabels({'-2\pi','0','2\pi'})
xlabel('Theta'); 
ylabel('Phase');
title("Charge 6 "+type5 +" Phase Profile")

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

% plot
% close all
% figure(12)
subplot(2,2,1)
hold on
xdata = (0:1:N-1);
% plot(f,myfftvortex,'LineStyle',linVec{i},'LineWidth',1,'Color',[0 0.5 0.8]);
plot(f,myfftstaircase,'LineStyle',linVec{i},'LineWidth',1,'Color',[0 0.5 0.8]);
plot(f,myfftclassical,'LineStyle',linVec{i},'LineWidth',1,'Color',[0.5 0 0.8]);
% plot(f,myfftfrench,'LineStyle',linVec{i},'LineWidth',1,'Color',[0.5 0.8 0]);
plot(f,myfftmcmc,'LineStyle',linVec{i},'LineWidth',1,'Color',[0.5 0.8 0]);
if lambda == 1
    plot(f,myfftvortex,'LineStyle',linVec{i},'LineWidth',1,'Color','k');
end
ax = gca;
xlim([-10 20]);
ylim([1E-6 10])
ax.FontSize = 12;
ax.LineWidth = 2;
legend(type5,type2,type4)
legend('Location','northwest')


line_type = ['-',"--", ":"];

line_names = ["design wavelength", "0.95 lambda factor", "1.05 lambda factor"];
for j =1:length(line_type)
    plot([NaN NaN], [NaN NaN],line_type(j), 'Color', 'k', 'DisplayName', line_names(j))
end

hold off
set(gca, 'YScale', 'log')
ylabel('|C_m|^2'); 
xlabel('Mode');
title("Modal Decomposition for Charge 6 SVCs")

