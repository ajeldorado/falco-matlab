% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to Zernike modes.
% 
%
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2_array = falco_get_Zernike_sensitivities(mp)

%--New fields in structure "mp.eval" that are used:
indsZnoll = mp.eval.indsZnoll;
Rsens = mp.eval.Rsens;

%% Zernikes
Nzern = length(indsZnoll);

rmsZvec = ones(size(indsZnoll))*1e-9; %--RMS values for each Zernike specified in vector indsZnoll [meters] 

ZmapCube = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.

%--Make sure ZmapCube is the right array size
if(size(ZmapCube,1)~=mp.P1.compact.Narr)
    ZmapCubeTemp = zeros(mp.P1.compact.Narr,mp.P1.compact.Narr);
    for zi=1:size(ZmapCube,3)
        ZmapCubeTemp(:,:,zi) = padOrCropEven(ZmapCube(:,:,zi),mp.P1.compact.Narr);
    end
    ZmapCube = ZmapCubeTemp; 
    clear ZmapCubeTemp
end

%% Make scoring masks

Nannuli = size(Rsens,1);
masks = zeros(size(mp.Fend.corr.mask,1), size(mp.Fend.corr.mask,2), Nannuli);

for ni = 1:Nannuli
    %--Make scoring masks for the annular regions
    maskCorr.pixresFP = mp.Fend.res;
    maskCorr.rhoInner = Rsens(ni,1); %--lambda0/D
    maskCorr.rhoOuter = Rsens(ni,2) ; %--lambda0/D
    maskCorr.angDeg = mp.Fend.corr.ang; %--degrees
    maskCorr.centering = mp.centering;
    maskCorr.FOV = mp.Fend.FOV;
    maskCorr.whichSide = 'both'; %--which (sides) of the dark hole have open
    %--Compact Model: Generate Software Mask for Correction 
    [masks(:,:,ni), ~, ~] = falco_gen_SW_mask(maskCorr);
end

%% Get nominal, unaberrated final E-field at each wavelength.

mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp); %--Input E-fields

E0cube = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.Nsbp); %--initialize

for si = 1:mp.Nsbp   
    modvar.sbpIndex = si;
    modvar.whichSource = 'star';     
    E0cube(:,:,si) = model_compact(mp, modvar);
end

%% Get aberrated, final focal-plane E-fields for each Zernike mode

EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,Nzern,mp.Nsbp); %--initialize
dEZarray = zeros(size(EZarray));

for iz = 1:Nzern
    for si = 1:mp.Nsbp      
        mp.P1.compact.E(:,:,si) = exp(1i*2*pi/mp.sbp_centers(si)*rmsZvec(iz)*ZmapCube(:,:,iz));%.*mp.P1.compact.E(:,:,si);       
        modvar.sbpIndex = si;
        modvar.whichSource = 'star';     
        EZarray(:,:,iz,si) = model_compact(mp, modvar);
        dEZarray(:,:,iz,si) = EZarray(:,:,iz,si)-E0cube(:,:,si);
    end
end

dE2cube = mean(abs(dEZarray).^2,4); % |dE|^2 averaged over wavelength

%% Compute Zernike sensitivity values averaged across the whole dark hole

dE2_array = zeros(Nzern,Nannuli);
for iz = 1:Nzern
   dEtemp = dE2cube(:,:,iz);
   
   for ni = 1:Nannuli
       dE2_array(iz,ni) = mean(dEtemp( logical(masks(:,:,ni))) );
   end
   
end

%--Print Zernike sensitivity results to command line
for iz = 1:Nzern
    fprintf('|dE|^2 at %dnm with %dnm RMS of    Z%d =',round(mp.lambda0*1e9),  round(1e9*rmsZvec(iz)), indsZnoll(iz) )
    
    for ni = 1:Nannuli
       mean(dEtemp( logical(masks(:,:,ni))) );
       fprintf('\t%.2e (%.1f-%.1f l/D)',dE2_array(iz,ni), Rsens(ni,1), Rsens(ni,2) )
    end
    fprintf('\n')
end

end %--END OF FUNCTION