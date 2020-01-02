% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%


function mp = falco_configure_dark_hole_region(mp)

%--Set Inputs
maskCorr.pixresFP = mp.Fend.res;
maskCorr.rhoInner = mp.Fend.corr.Rin; %--lambda0/D
maskCorr.rhoOuter = mp.Fend.corr.Rout ; %--lambda0/D
maskCorr.angDeg = mp.Fend.corr.ang; %--degrees
maskCorr.centering = mp.centering;
maskCorr.FOV = mp.Fend.FOV;
maskCorr.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
if(isfield(mp.Fend,'shape'));  maskCorr.shape = mp.Fend.shape;  end
if(isfield(mp.Fend,'clockAngDeg'));  maskCorr.clockAngDeg = mp.Fend.clockAngDeg;  end

%--Compact Model: Generate Software Mask for Correction 
[mp.Fend.corr.mask, mp.Fend.xisDL, mp.Fend.etasDL] = falco_gen_SW_mask(maskCorr); 
mp.Fend.corr.settings = maskCorr; %--Store values for future reference
%--Size of the output image 
%--Need the sizes to be the same for the correction and scoring masks
mp.Fend.Nxi  = size(mp.Fend.corr.mask,2);
mp.Fend.Neta = size(mp.Fend.corr.mask,1);

[XIS,ETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
mp.Fend.RHOS = sqrt(XIS.^2 + ETAS.^2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(mp.flagFiber)
    mp = falco_configure_fiber_dark_hole(mp);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%--Evaluation Model for Computing Throughput (same as Compact Model but
% with different Fend.resolution)
mp.Fend.eval.dummy = 1; %--Initialize the structure if it doesn't exist.
maskCorr.pixresFP = mp.Fend.eval.res; %--Assign the resolution
[mp.Fend.eval.mask, mp.Fend.eval.xisDL, mp.Fend.eval.etasDL] = falco_gen_SW_mask(maskCorr);  %--Generate the mask
mp.Fend.eval.Nxi  = size(mp.Fend.eval.mask,2);
mp.Fend.eval.Neta = size(mp.Fend.eval.mask,1);


% (x,y) location [lambda_c/D] in dark hole at which to evaluate throughput
[XIS,ETAS] = meshgrid(mp.Fend.eval.xisDL - mp.thput_eval_x, mp.Fend.eval.etasDL - mp.thput_eval_y);
mp.Fend.eval.RHOS = sqrt(XIS.^2 + ETAS.^2);

%--Storage array for throughput at each iteration
mp.thput_vec = zeros(mp.Nitr+1,1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%--Software Mask for Scoring Contrast 
%--Set Inputs
maskScore.rhoInner = mp.Fend.score.Rin; %--lambda0/D
maskScore.rhoOuter = mp.Fend.score.Rout ; %--lambda0/D
maskScore.angDeg = mp.Fend.score.ang; %--degrees
maskScore.centering = mp.centering;
maskScore.FOV = mp.Fend.FOV; %--Determines max dimension length
maskScore.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
if(isfield(mp.Fend,'shape'));  maskScore.shape = mp.Fend.shape;  end
if(isfield(mp.Fend,'clockAngDeg'));  maskScore.clockAngDeg = mp.Fend.clockAngDeg;  end

%--Compact Model: Generate Software Mask for Scoring Contrast 
maskScore.Nxi = mp.Fend.Nxi; %--Set min dimension length to be same as for corr 
maskScore.pixresFP = mp.Fend.res;
[mp.Fend.score.mask,~,~] = falco_gen_SW_mask(maskScore); 
mp.Fend.score.settings = maskScore; %--Store values for future reference

%--Number of pixels used in the dark hole
mp.Fend.corr.Npix = sum(sum(mp.Fend.corr.mask));
mp.Fend.score.Npix = sum(sum(mp.Fend.score.mask));

%--Indices of dark hole pixels and logical masks
if(mp.flagFiber && mp.flagLenslet)
    mp.Fend.corr.inds = find(sum(mp.Fend.lenslet.mask,3)~=0);
else
    mp.Fend.corr.inds = find(mp.Fend.corr.mask~=0);
end
mp.Fend.corr.maskBool = logical(mp.Fend.corr.mask);

mp.Fend.score.inds = find(mp.Fend.score.mask~=0);
mp.Fend.score.maskBool = logical(mp.Fend.score.mask);


end