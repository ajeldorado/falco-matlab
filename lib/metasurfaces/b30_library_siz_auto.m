% CREATED BY: LORENZO KONIG 
% LAST MODIFIED: JAN 2023

% This script reads in a desired phase map and a specific wavelength index
% It finds the corresponding metasurface block size for that phase
% And returns the effective phase at the specific wavelength for
% that metasurface blocksize

function [m,n,m1]=b30_library_siz_auto(phi,n_lam)
    % have to give phase AND wavelength index !
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_diam_LVODCA_sw200_crc1' % to be able to use the data from there
    %load('rev1cd_phase_function_LVODCA_crc1diam_L1_sw200.mat');
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_Si_LVODCA_sw0_crc1' % to be able to use the data from there
    %load('rev1cd_phase_function_LVODCA_crc1Si_L1_sw0.mat');
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_diam_LVODCA_sw200_crc1_COS' % to be able to use the data from there
    %load('rev1cd_phase_function_LVODCA_crc1diam_L1_sw200_COS.mat');
    % ALL USELESS FOR NOW since Lb=1.00 for diam too small !
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_diam_LVODCA_sw200_crc1_Lb1.42_V_no-sw-touch RELEVANT' % to be able to use the data from there
    %load('rev1cd_phase_function_LVODCA_crc1diam_L1_sw200_V.mat'); % using correct period, but very large posts
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_diam_LVODCA_sw200_crc1_Lb1.42_V' % to be able to use the data from there
    %load('rev1cd_phase_function_LVODCA_crc1diam_L1_sw200_V.mat'); % using correct period, and smaller posts such that sw do not touch
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_diam_LVODCA_sw200_crc1_Lb1.42_COS RELEVANT' % to be able to use the data from there
    %load('rev1cd_phase_function_LVODCA_crc1diam_L1_sw200_COS.mat'); % using correct period, and smaller posts such that sw do not touch
    %addpath '/Users/lorenzo/Documents/RETICOLO_V9/rev1cd_diam_LVODCA_sw200_crc1_Lb1.42_COS_no-sw-touch RELEVANT' % to be able to use the data from there
    %load('r20_optimize_rev1cd_phase_function_LVODCA_crc1diam_L1_sw200_COS.mat'); % using correct period, and smaller posts such that sw do not touch
%     load('r20_optimize.mat'); % using correct period, and smaller posts such that sw do not touch
    load('r20_optimize_new.mat'); % using correct period, and smaller posts such that sw do not touch
    load('design.mat','Nlb');
    %phi_c=permute(phi_array(min_ih,min_isn:min_ism+1,:),[2,3,1]);%[    2.8697    2.9464    3.0439   -3.1189+2*pi   -2.9735+2*pi   -2.8010+2*pi   -2.5988+2*pi   -2.3629+2*pi   -2.0881+2*pi   -1.7702+2*pi   -1.4076+2*pi   -1.0038+2*pi   -0.5662+2*pi   -0.1005+2*pi    0.3954+2*pi    0.9276+2*pi    1.4943+2*pi    2.0790+2*pi    2.6708+2*pi   -2.9977+4*pi   -2.3900+4*pi   -1.7801+4*pi   -1.1669+4*pi   -0.6204+4*pi   -0.1226+4*pi    0.3313+4*pi    0.7253+4*pi    1.0519+4*pi    1.3382+4*pi];
    phi_c=permute(phase_function(:,:,:),[2,3,1]);%[    2.8697    2.9464    3.0439   -3.1189+2*pi   -2.9735+2*pi   -2.8010+2*pi   -2.5988+2*pi   -2.3629+2*pi   -2.0881+2*pi   -1.7702+2*pi   -1.4076+2*pi   -1.0038+2*pi   -0.5662+2*pi   -0.1005+2*pi    0.3954+2*pi    0.9276+2*pi    1.4943+2*pi    2.0790+2*pi    2.6708+2*pi   -2.9977+4*pi   -2.3900+4*pi   -1.7801+4*pi   -1.1669+4*pi   -0.6204+4*pi   -0.1226+4*pi    0.3313+4*pi    0.7253+4*pi    1.0519+4*pi    1.3382+4*pi];
    phi_b=phi_c;%+1000;
    for j=1:length(phi_b(1,:))
      add=0;
      for i=1:length(phi_b)-1
        phi_b(i,j)=phi_b(i,j)+add;
        if (phi_c(i+1,j)-phi_c(i,j))<0
            add=add+2*pi; % if there is a phase jump 0-2pi, add 2pi
        end
      end
      phi_b(end,j)=phi_b(end,j)+add;
    end
    phi_a=(phi_b-phi_b(1,:)-pi); % in rad
    T_a=permute(T_function(:,:,:),[2,3,1]);
    m1=interp1(phi_a(:,(Nlb+1)/2),s_vals',phi); % this gives a post size % use wvl No 2 as perfect reference... :')
    m=interp1(s_vals',phi_a(:,n_lam),m1); % THIS GIVES A PHASE UNLIKE THE OTHER SIZ.M FUNCTIONS!!!
    n=interp1(s_vals',T_a(:,n_lam),m1); % this gives a transmission value 
end