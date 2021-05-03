% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office of Technology
% Transfer at the California Institute of Technology.
% -----------------------------------------------------------------------

function mp = falco_configure_dark_hole_region(mp)


CORR.pixresFP = mp.Fend.res;
CORR.centering = mp.centering;
if(isfield(mp.Fend,'FOV'));  CORR.FOV = mp.Fend.FOV;  end
if(isfield(mp.Fend,'xiFOV'));  CORR.xiFOV = mp.Fend.xiFOV;  end
if(isfield(mp.Fend,'etaFOV'));  CORR.etaFOV = mp.Fend.etaFOV;  end
if(isfield(mp.Fend,'Nxi'));  CORR.Nxi = mp.Fend.Nxi;  end
if(isfield(mp.Fend,'Neta'));  CORR.Neta = mp.Fend.Neta;  end


%% Correction Region

Nzones = length(mp.Fend.corr.Rin);

% Will treat all inputs as arrays. If only one dark hole zone, then convert
% strings into cell arrays to be iterable.
if Nzones == 1
    sides{1} = mp.Fend.sides;
    if isfield(mp.Fend, 'shape'); shapes{1} = mp.Fend.shape; else; shapes{1} = 'circle'; end

else
    sides = mp.Fend.sides;
    shapes = mp.Fend.shape;
    if isfield(mp.Fend, 'shape')
        shapes = mp.Fend.shape;
    else
        for ii = 1:Nzones
            shapes{ii} = 'circle';
        end
    end
end

maskCorr = zeros(1, 1);
for iZone = 1:Nzones
    CORR.rhoInner = mp.Fend.corr.Rin(iZone); %--lambda0/D
    CORR.rhoOuter = mp.Fend.corr.Rout(iZone); %--lambda0/D
    CORR.angDeg = mp.Fend.corr.ang(iZone); %--degrees
    CORR.whichSide = sides{iZone}; %--which (sides) of the dark hole have open
    if(isfield(mp.Fend,'shape'));  CORR.shape = shapes{iZone};  end
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

mp.Fend.corr.mask = logical(mp.Fend.corr.mask);
mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);

% Size of the output image 
mp.Fend.Nxi  = size(mp.Fend.corr.mask, 2);
mp.Fend.Neta = size(mp.Fend.corr.mask, 1);

[XIS, ETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
mp.Fend.RHOS = sqrt(XIS.^2 + ETAS.^2);

if mp.flagFiber
    mp = falco_configure_fiber_dark_hole(mp);
end

%% Evaluation Model for Computing Throughput (just need size and coordinates, not mask)

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


%% Scoring Region

Nzones = length(mp.Fend.score.Rin);

SCORE.Nxi = mp.Fend.Nxi; % Set same dimensions as for correction region 
SCORE.Neta = mp.Fend.Neta;
SCORE.pixresFP = mp.Fend.res;
SCORE.centering = mp.centering;

maskScore = 0;
for iZone = 1:Nzones
    SCORE.rhoInner = mp.Fend.score.Rin(iZone); %--lambda0/D
    SCORE.rhoOuter = mp.Fend.score.Rout(iZone); %--lambda0/D
    SCORE.angDeg = mp.Fend.score.ang(iZone); %--degrees
    SCORE.whichSide = sides{iZone}; %--which (sides) of the dark hole have open
    if(isfield(mp.Fend,'shape'));  SCORE.shape = shapes{iZone};  end
    if(isfield(mp.Fend,'clockAngDeg'));  SCORE.clockAngDeg = mp.Fend.clockAngDeg(iZone);  end
    if(isfield(mp.Fend,'xiOffset'));  SCORE.xiOffset = mp.Fend.xiOffset(iZone);  end
    if(isfield(mp.Fend,'etaOffset'));  SCORE.etaOffset = mp.Fend.etaOffset(iZone);  end

    [maskTemp, ~, ~] = falco_gen_SW_mask(SCORE);
    maskScore = maskScore + maskTemp;
end
mp.Fend.score.mask = maskScore;

mp.Fend.score.mask = logical(mp.Fend.score.mask);
mp.Fend.score.maskBool = logical(mp.Fend.score.mask);

%--Number of pixels used in the dark hole
mp.Fend.corr.Npix = sum(mp.Fend.corr.mask(:));
mp.Fend.score.Npix = sum(mp.Fend.score.mask(:));

%--Indices of dark hole pixels and logical masks
if mp.flagFiber && mp.flagLenslet
    mp.Fend.corr.inds = find(sum(mp.Fend.lenslet.mask, 3) ~= 0);
end

%%
mp.Fend.scoreInCorr = mp.Fend.score.mask(mp.Fend.corr.mask); % vector indicating which pixels in vectorized correction region are also in the scoring region

end