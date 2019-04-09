% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate the LUVOIR Design A, Aperture 5, telescope pupil in
%   Matlab using PROPER
% Coordinates and dimensions of the primary, secondary, and hex segments
%   are from Matthew Bolcar (NASA GSFC).
% Coordinates and dimenstions of the secondary mirror support struts were a
%   best-fit match by A.J. Riggs by matching PROPER-made rectangles to the 
%   pupil file from Matthew Bolcar (NASA GSFC).
%
% Modified on 2018-02-25 by A.J. Riggs to be for LUVOIR A aperture 5. 
% Written on  2017-09-07 by A.J. Riggs to generate the first proposed LUVOIR pupil. 
%   Values for the geometry were provided by Matthew Bolcar at NASA GSFC.
%
%--Coordinates of hex segments to skip:
% 1 13 114 115 126 127
% 1 12 113 114 125 126

clear all

%%--Add to the MATLAB Path
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

inputs.Nbeam = 1000;
inputs.magfacD = 1;
inputs.wStrut = 1/100;


pupil = falco_gen_pupil_LUVOIR_A_final(inputs);

figure(2); imagesc(pupil); axis xy equal tight; title('Input Pupil','Fontsize',20); colorbar;

figure(3); imagesc(pupil(2:end,2:end)-fliplr(pupil(2:end,2:end))); axis xy equal tight; title('Symmetry Check: Differencing','Fontsize',20); colorbar
