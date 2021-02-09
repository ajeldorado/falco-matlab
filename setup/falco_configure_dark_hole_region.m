% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office of Technology
% Transfer at the California Institute of Technology.
% -----------------------------------------------------------------------
%
%

function mp = falco_configure_dark_hole_region(mp)

%% Correction Region

CORR.pixresFP = mp.Fend.res;
CORR.centering = mp.centering;
if(isfield(mp.Fend,'FOV'));  CORR.FOV = mp.Fend.FOV;  end
if(isfield(mp.Fend,'xiFOV'));  CORR.xiFOV = mp.Fend.xiFOV;  end
if(isfield(mp.Fend,'etaFOV'));  CORR.etaFOV = mp.Fend.etaFOV;  end
if(isfield(mp.Fend,'Nxi'));  CORR.Nxi = mp.Fend.Nxi;  end
if(isfield(mp.Fend,'Neta'));  CORR.Neta = mp.Fend.Neta;  end

Nzones = length(mp.Fend.corr.Rin);

if Nzones == 1
    CORR.rhoInner = mp.Fend.corr.Rin; %--lambda0/D
    CORR.rhoOuter = mp.Fend.corr.Rout ; %--lambda0/D
    CORR.angDeg = mp.Fend.corr.ang; %--degrees
    CORR.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
    if(isfield(mp.Fend,'shape'));  CORR.shape = mp.Fend.shape;  end
    if(isfield(mp.Fend,'clockAngDeg'));  CORR.clockAngDeg = mp.Fend.clockAngDeg;  end
    if(isfield(mp.Fend,'xiOffset'));  CORR.xiOffset = mp.Fend.xiOffset;  end
    if(isfield(mp.Fend,'etaOffset'));  CORR.etaOffset = mp.Fend.etaOffset;  end
    
    [mp.Fend.corr.mask, mp.Fend.xisDL, mp.Fend.etasDL] = falco_gen_SW_mask(CORR); 

elseif Nzones > 1
    maskCorr = zeros(1, 1);
    for iZone = 1:Nzones
        CORR.rhoInner = mp.Fend.corr.Rin(iZone); %--lambda0/D
        CORR.rhoOuter = mp.Fend.corr.Rout(iZone); %--lambda0/D
        CORR.angDeg = mp.Fend.corr.ang(iZone); %--degrees
        CORR.whichSide = mp.Fend.sides{iZone}; %--which (sides) of the dark hole have open
        if(isfield(mp.Fend,'shape'));  CORR.shape = mp.Fend.shape{iZone};  end
        if(isfield(mp.Fend,'clockAngDeg'));  CORR.clockAngDeg = mp.Fend.clockAngDeg(iZone);  end
        if(isfield(mp.Fend,'xiOffset'));  CORR.xiOffset = mp.Fend.xiOffset(iZone);  end
        if(isfield(mp.Fend,'etaOffset'));  CORR.etaOffset = mp.Fend.etaOffset(iZone);  end
        
        % Combine multiple zones. Use the largest array size
        [maskTemp, ~, ~] = falco_gen_SW_mask(CORR);
        Nrow = max([size(maskTemp, 1), size(maskCorr, 1)]);
        Ncol = max([size(maskTemp, 2), size(maskCorr, 2)]);
        maskCorr = pad_crop(maskCorr, [Nrow, Ncol]) + pad_crop(maskTemp, [Nrow, Ncol]);
    end
    mp.Fend.corr.mask = maskCorr;

    CORR.Nxi = size(maskCorr, 2);
    CORR.Neta = size(maskCorr, 1);
    [~, mp.Fend.xisDL, mp.Fend.etasDL] = falco_gen_SW_mask(CORR); % generate coordinates
end

% Size of the output image 
mp.Fend.Nxi  = size(mp.Fend.corr.mask, 2);
mp.Fend.Neta = size(mp.Fend.corr.mask, 1);

[XIS, ETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
mp.Fend.RHOS = sqrt(XIS.^2 + ETAS.^2);

if mp.flagFiber
    mp = falco_configure_fiber_dark_hole(mp);
end

%% Evaluation Model for Computing Throughput (just need size and coordinates, not map)

mp.Fend.eval.dummy = 1; %--Initialize the structure if it doesn't exist.
CORR.pixresFP = mp.Fend.eval.res; %--Assign the resolution
CORR.Nxi = ceil_even(mp.Fend.eval.res/mp.Fend.res*mp.Fend.Nxi);
CORR.Neta = ceil_even(mp.Fend.eval.res/mp.Fend.res*mp.Fend.Neta);
[~, mp.Fend.eval.xisDL, mp.Fend.eval.etasDL] = falco_gen_SW_mask(CORR);  %--Generate the mask
mp.Fend.eval.Nxi  = CORR.Nxi;
mp.Fend.eval.Neta = CORR.Neta;

% (x,y) location [lambda_c/D] in dark hole at which to evaluate throughput
[XIS,ETAS] = meshgrid(mp.Fend.eval.xisDL - mp.thput_eval_x, mp.Fend.eval.etasDL - mp.thput_eval_y);
mp.Fend.eval.RHOS = sqrt(XIS.^2 + ETAS.^2);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Scoring Region

SCORE.Nxi = mp.Fend.Nxi; % Set same dimensions as for correction region 
SCORE.Neta = mp.Fend.Neta;
SCORE.pixresFP = mp.Fend.res;
SCORE.centering = mp.centering;
% if(isfield(mp.Fend,'FOV'));  SCORE.FOV = mp.Fend.FOV;  end
% if(isfield(mp.Fend,'xiFOV'));  SCORE.xiFOV = mp.Fend.xiFOV;  end
% if(isfield(mp.Fend,'etaFOV'));  SCORE.etaFOV = mp.Fend.etaFOV;  end

if Nzones == 1
    SCORE.rhoInner = mp.Fend.score.Rin; %--lambda0/D
    SCORE.rhoOuter = mp.Fend.score.Rout ; %--lambda0/D
    SCORE.angDeg = mp.Fend.score.ang; %--degrees
    SCORE.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
    if(isfield(mp.Fend,'shape'));  SCORE.shape = mp.Fend.shape;  end
    if(isfield(mp.Fend,'clockAngDeg'));  SCORE.clockAngDeg = mp.Fend.clockAngDeg;  end
    if(isfield(mp.Fend,'xiOffset'));  SCORE.xiOffset = mp.Fend.xiOffset;  end
    if(isfield(mp.Fend,'etaOffset'));  SCORE.etaOffset = mp.Fend.etaOffset;  end

    [mp.Fend.score.mask, ~, ~] = falco_gen_SW_mask(SCORE); 
    
elseif Nzones > 1
    
    maskScore = 0;
    for iZone = 1:Nzones
        SCORE.rhoInner = mp.Fend.score.Rin(iZone); %--lambda0/D
        SCORE.rhoOuter = mp.Fend.score.Rout(iZone); %--lambda0/D
        SCORE.angDeg = mp.Fend.score.ang(iZone); %--degrees
        SCORE.whichSide = mp.Fend.sides{iZone}; %--which (sides) of the dark hole have open
        if(isfield(mp.Fend,'shape'));  SCORE.shape = mp.Fend.shape{iZone};  end
        if(isfield(mp.Fend,'clockAngDeg'));  SCORE.clockAngDeg = mp.Fend.clockAngDeg(iZone);  end
        if(isfield(mp.Fend,'xiOffset'));  SCORE.xiOffset = mp.Fend.xiOffset(iZone);  end
        if(isfield(mp.Fend,'etaOffset'));  SCORE.etaOffset = mp.Fend.etaOffset(iZone);  end
        
        [maskTemp, ~, ~] = falco_gen_SW_mask(SCORE);
        maskScore = maskScore + maskTemp;
    end
    mp.Fend.score.mask = maskScore;
    
end

%--Number of pixels used in the dark hole
mp.Fend.corr.Npix = sum(mp.Fend.corr.mask(:));
mp.Fend.score.Npix = sum(mp.Fend.score.mask(:));

%--Indices of dark hole pixels and logical masks
if mp.flagFiber && mp.flagLenslet
    mp.Fend.corr.inds = find(sum(mp.Fend.lenslet.mask,3) ~= 0);
else
    mp.Fend.corr.inds = find(mp.Fend.corr.mask~=0);
end
mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);

mp.Fend.score.inds = find(mp.Fend.score.mask~=0);
mp.Fend.score.maskBool = logical(mp.Fend.score.mask);


end