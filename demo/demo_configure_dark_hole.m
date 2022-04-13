% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to demo ways to define different dark holes for FALCO.
%
% This shows how to define the necessary (and some optional) values in the
% mp structure that falco_configure_dark_hole_region() needs to define the
% dark hole masks stored in mp.Fend.corr.maskBool and mp.Fend.score.maskBool
%
% Note: All possible options for mp.Fend.sides, the variable for which
% side(s) of the star to put the dark hole, are:
% 'left', 'right', 'top', 'up', 'bottom', 'down', 'lr', 'rl', 'leftright', 'rightleft', 'tb', 'bt', 'ud', 'du', 'topbottom', 'bottomtop', 'updown', 'downup'

%% Annular Dark Hole
clear mp

%--Correction region definition 
mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

%--Scoring region definition
mp.Fend.score = mp.Fend.corr;

mp.centering = 'pixel';
mp.Fend.sides = 'leftright';
mp.Fend.shape = 'circle';
mp.Fend.res = 10;
% mp.Fend.xiOffset = 6;

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(11); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;

area = sum(sum(mp.Fend.corr.maskBool));
areaExpected = pi*(mp.Fend.corr.Rout^2 - mp.Fend.corr.Rin^2)*(2*mp.Fend.corr.ang/360)*(mp.Fend.res^2);

% import matlab.unittest.constraints.IsEqualTo
% import matlab.unittest.constraints.RelativeTolerance
% testCase.verifyThat(area, IsEqualTo(areaExpected,'Within', RelativeTolerance(0.001)))

%% Half Annulus Dark Hole
clear mp

%--Correction region definition 
mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
mp.Fend.clockAngDeg = 10; % Clocking angle of the dark hole region

%--Scoring region definition
mp.Fend.score = mp.Fend.corr;

mp.centering = 'pixel';
mp.Fend.sides = 'bottom';
mp.Fend.shape = 'circle';
mp.Fend.res = 10;
% mp.Fend.xiOffset = 6;

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(12); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;


%% D-shape Dark Hole
clear mp

%--Correction region definition 
mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
mp.Fend.clockAngDeg = 10; % Clocking angle of the dark hole region

%--Scoring region definition
mp.Fend.score = mp.Fend.corr;

mp.centering = 'pixel';
mp.Fend.sides = 'top';
mp.Fend.shape = 'D';
mp.Fend.res = 10;
% mp.Fend.xiOffset = 6;

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(13); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;


%% Rectangular Dark Hole
clear mp

%--Correction region definition 
mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
mp.Fend.clockAngDeg = 10; % Clocking angle of the dark hole region

%--Scoring region definition
mp.Fend.score = mp.Fend.corr;

mp.Fend.shape = 'rect';
mp.Fend.sides = 'updown';
mp.centering = 'pixel';
mp.Fend.res = 10;
mp.Fend.Nxi = 250;
mp.Fend.Neta = 400;
% mp.Fend.xiOffset = 6;

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(14); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;


%% Square Dark Hole (with inner circle cut out)
clear mp

%--Correction region definition 
mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
mp.Fend.clockAngDeg = 10; % Clocking angle of the dark hole region

%--Scoring region definition
mp.Fend.score = mp.Fend.corr;

mp.Fend.shape = 'square';
mp.Fend.sides = 'updown';
mp.centering = 'pixel';
mp.Fend.res = 10;
mp.Fend.FOV = 15;
% mp.Fend.xiOffset = 6;

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(15); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;


%% Offset Square Dark Hole
clear mp

%--Correction region definition 
mp.Fend.corr.Rin = 0;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = 4.5;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
mp.Fend.clockAngDeg = 0; % Clocking angle of the dark hole region

%--Scoring region definition
mp.Fend.score = mp.Fend.corr;

mp.Fend.shape = 'square';
mp.Fend.sides = 'updown';
mp.centering = 'pixel';
mp.Fend.res = 10;
mp.Fend.Nxi = 250;
mp.Fend.Neta = 150;
mp.Fend.xiOffset = 6;
mp.Fend.etaOffset = -2;

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(16); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;


%% Multi-Zone, Multi-Shape Dark Hole
% Use this to combine multiple shapes, which can overlap.
clear mp

%--Correction and scoring region definition
mp.Fend.corr.Rin = [2, 2];   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout  = [5, 5];  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = [150, 180];  % angular opening of dark hole correction region [degrees]
mp.Fend.score = mp.Fend.corr;
mp.Fend.FOV = 30;

mp.centering = 'pixel';
mp.Fend.sides = {'lr', 'lr'};
mp.Fend.shape = {'circle', 'square'};
mp.Fend.res = 10;
mp.Fend.xiOffset = [0, 20];

mp.Fend.eval.res = 20;
mp.flagFiber = false;
mp.thput_eval_x = 7;
mp.thput_eval_y = 0;

mp = falco_configure_dark_hole_region(mp);

figure(17); imagesc(mp.Fend.xisDL, mp.Fend.etasDL, mp.Fend.corr.maskBool); axis xy equal tight; colormap gray; drawnow;

area = sum(sum(mp.Fend.corr.maskBool));
areaExpected = pi*(mp.Fend.corr.Rout(1)^2 - mp.Fend.corr.Rin(1)^2)*(2*mp.Fend.corr.ang(1)/360)*(mp.Fend.res^2) + ...
    (4*mp.Fend.corr.Rout(2)^2 - pi*mp.Fend.corr.Rin(2)^2)*(mp.Fend.res^2);
