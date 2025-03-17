% -------------------------------------------------------------------------
%
% Script to optimize thicknesses for a 2 material vortex achromat


% REVISION HISTORY:
% --------------
% Created on 2025-1-3 by Niyati Desai.

% optimize thicknesses for two materials

% lam input should be in microns
% Material options are (included in /lib_external/materials folder):
% aSi_k
% aSi (alpha silicon)
% cSi (crystalline silicon)
% diamond
% GaP1
% HfO2
% Si3N4
% SiO2
% TiO2
% ZnSe
% fusedsilica

function [thicknesses, phasescalefacs] = optimizeachromat(mat1,mat2,lams)

% mat1 = "diamond";
% mat2 = "TiO2";
% lams = linspace(0.6,0.72,11);
n1 = getn(lams,mat1);
n2 = getn(lams,mat2);
%n2 = ones([1,length(lams)]);

fun = @(x)sum((((x(1)*(n1-1) + x(2)*(n2-1))./lams)-1).^2);
x0 = [0,0]; %[3.5,-2.5];
x = fminsearch(fun,x0)
%fmincon fminbnd 
thicknesses = x;


subcharge1=x(1)*(n1-1)./lams
subcharge2=x(2)*(n2-1)./lams

mysum=subcharge1+subcharge2
phasescalefacs = mysum;

figure()
plot(lams,mysum)
title(mat1+" and "+mat2)
xlabel("wavelength (microns)")
ylabel("phase scale factor")
