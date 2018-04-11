% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to generate a multi-ring shaped pupil mask in Matlab using PROPER.
%
% Modified by A.J. Riggs on 2017-10-24 to be a function.
% Written by A.J. Riggs on 2017-09-07. 
%
% INPUTS: 
%  rEdgesLeft:  vector of leading (i.e., rising) radial coordinate values of each transmissive ring
%  rEdgesRight: vector of trailing (i.e., falling) radial coordinate values of each transmissive ring
%  dx:      spatial resolution for a pixel (meters)
%  Dsp:     diameter of the beam at the mask (meters)
%  centering:   centering of beam in array. Either 'pixel' or 'interpixel'
%  (optional, to be added later if desired) inputs.flagBinary: flag to specify if pupil should be binary
%
% OUTPUTS:
%  CRMcrop:     2-D square array for the concentric ring mask. Cropped down to the smallest even-sized array with no extra zero padding. 

function CRMcrop = falco_gen_multi_ring_SP(rEdgesLeft,rEdgesRight,dx,Dbeam,centering)


% %--DEBUGGING ONLY: HARD-CODED INPUTS
% clear;
% addpath ~/Repos/FALCO/proper_v3.0.1_matlab_22aug17/
% cd ~/Repos/FALCO/masks/ringSP
% rEdgesLeft  = load('out_RSPLC1D_A_maxTrPH_3115Dpup9860_28WA210_41LS79_10FPres5_BW10N11_c90_Nring10_left.dat');
% rEdgesRight = load('out_RSPLC1D_A_maxTrPH_3115Dpup9860_28WA210_41LS79_10FPres5_BW10N11_c90_Nring10_right.dat');
% dx = .1e-3;%70e-3/1000;
% Dsp = 48e-3;
% centering = 'pixel'; % 'pixel'/'odd' or 'interpixel'/'even'


diam = Dbeam;% width of the beam (m)
NapAcross = Dbeam/dx; % minimum even number of points across to fully contain the actual aperture (if interpixel centered)
if(strcmpi(centering,'pixel'))
    Narray = 2*ceil(1/2*(Dbeam/dx+1/2)); %--number of points across output array. Sometimes requires two more pixels when pixel centered.
else
    Narray = 2*ceil(1/2*Dbeam/dx); %--number of points across output array. Same size as width when interpixel centered.
end

Darray = Narray*dx; %--width of the output array (meters)
bdf = 1; %--beam diameter factor in output array
wl_dummy   = 1e-6;     % wavelength (m); Dummy value--no propagation here, so not used.

switch centering % 0 shift for pixel-centered pupil, or -diam/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -Darray/Narray; 
    case {'pixel'}
        cshift = 0;
end

Nrings = length(rEdgesLeft);


%--Whether to threshold the amplitudes and make the mask binary
% if(isfield(inputs,'flagBinary'))
%     flagBinary = inputs.flagBinary;
% else
    flagBinary = false;
% end

switch centering % 0 shift for pixel-centered pupil, or -diam/2/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -diam/2/Narray; 
    case {'pixel'}
        cshift = 0;
end

% if(Narray>=1000) %--Binary cutoff (dependent on resolution)
%    thresh = 0.95;
% elseif Narray >= 500
%    thresh = 0.75;
% elseif Narray<500
%    thresh = 0.50;
% end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




%--CREATE THE CONCENTRIC RING SHAPED PUPIL ONE RING AT A TIME
ring_cube = zeros(Narray,Narray,Nrings); %--Storage array

for ni=Nrings:-1:1;%length(od_vec)

%--(RE-)INITIALIZE PROPER
bm = prop_begin(Darray, wl_dummy, Narray,'beam_diam_fraction',bdf);
% figure(1); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

%--Outer diameter of ring
ra_OD = (rEdgesRight(ni)*diam); 
cx_OD = 0 + cshift;
cy_OD = 0 + cshift;

bm = prop_circular_aperture(bm, ra_OD,'XC',cx_OD,'YC',cy_OD);%, cx, cy, norm);
% figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;
% figure(3); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;


%--Inner diameter of ring
ra_ID = (rEdgesLeft(ni)*diam);
cx_ID = 0 + cshift;
cy_ID = 0 + cshift;
bm = prop_circular_obscuration(bm, ra_ID,'XC',cx_ID,'YC',cy_ID);%, cx, cy, norm)
% figure(4); imagesc(ifftshift(abs(bm.wf))); axis xy equal tight; colorbar;

ap = ifftshift(abs(bm.wf));
apCrop = padOrCropEven(ap,Narray);
ring_cube(:,:,ni) = apCrop;

% figure(10); imagesc(ap); axis xy equal tight; colorbar; pause(1/20)
    
end

CRM = sum(ring_cube,3);
% figure(21); imagesc(CRM); axis xy equal tight; colorbar;

%--Crop down the mask to get rid of extra zero padding. Speeds up MFTs.
CRMsum = sum(CRM(:));
CRMdiff = 0; counter = 2;
while( abs(CRMdiff) <= 1e-7)
    NcrmCrop = length(CRM)-counter; %--Number of points across the cropped-down Lyot stop
    CRMdiff = CRMsum - sum(sum( padOrCropEven(CRM, NcrmCrop-2) )); %--Subtract an extra 2 to negate the extra step that overshoots.
    counter = counter + 2;
end
CRMcrop = padOrCropEven(CRM,NcrmCrop); %--The cropped-down Lyot stop for the compact model       

% figure(22); imagesc(CRMcrop); axis xy equal tight; colorbar;


end %--END OF FUNCTION
