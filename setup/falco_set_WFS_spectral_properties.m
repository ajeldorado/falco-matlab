% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Spectral weighting of images. 
%
%
% REVISION HISTORY:
% --------------
% Created by G. Ruane on 2020-04-21 based on falco_set_spectral_properties
% ---------------

function mp = falco_set_WFS_spectral_properties(mp)

if(~isfield(mp.wfs,'lambda0')); mp.wfs.lambda0 = mp.lambda0; end
if(~isfield(mp.wfs,'Nsbp')); mp.wfs.Nsbp = 1; end
if(~isfield(mp.wfs,'Nwpsbp')); mp.wfs.Nwpsbp = 1; end
if(~isfield(mp.wfs,'fracBW')); mp.wfs.fracBW = 0.01; end

lambda0 = mp.wfs.lambda0;
Nsbp = mp.wfs.Nsbp;
Nwpsbp = mp.wfs.Nwpsbp;
fracBW = mp.wfs.fracBW;

%% Compact Model

%--Center(-ish) wavelength indices (ref = reference). (Only the center if
%  an odd number of wavelengths is used.)
mp.wfs.si_ref = ceil(Nsbp/2);

%--Wavelengths used for Compact Model (and Jacobian Model)
mp.wfs.sbp_weights = ones(Nsbp,1);
if(mp.wfs.flagSim) %--Set ctrl wavelengths evenly between endpoints (inclusive) of the total bandpass.
    if(Nsbp==1)
        mp.wfs.sbp_centers = lambda0;
    else
        mp.wfs.sbp_centers = lambda0*linspace(1-fracBW/2,1+fracBW/2,Nsbp);
        mp.wfs.sbp_weights(1) = 1/2; %--Give end sub-bands half weighting
        mp.wfs.sbp_weights(end) = 1/2; %--Give end sub-bands half weighting
    end
else %--For cases with multiple sub-bands: Choose wavelengths to be at subbandpass centers since the wavelength samples will span to the full extent of the sub-bands.
    mp.wfs.fracBWcent2cent = fracBW*(1-1/Nsbp); %--Bandwidth between centers of endpoint subbandpasses.
    mp.wfs.sbp_centers = lambda0*linspace(1-mp.wfs.fracBWcent2cent/2,1+mp.wfs.fracBWcent2cent/2,Nsbp); %--Space evenly at the centers of the subbandpasses.
end
mp.wfs.sbp_weights = mp.wfs.sbp_weights/sum(mp.wfs.sbp_weights); %--Normalize the sum of the weights

fprintf(' WFS using %d discrete wavelength(s) in each of %d sub-bandpasses over a %.1f%% total bandpass \n', Nwpsbp, Nsbp, 100*fracBW);
fprintf('WFS sub-bandpasses are centered at wavelengths [nm]:\t '); fprintf('%.2f  ',1e9*mp.wfs.sbp_centers); fprintf('\n\n');

%% Bandwidth and Wavelength Specs: Full Model

%--Center(-ish) wavelength indices (ref = reference). (Only the center if an odd number of wavelengths is used.)
mp.wfs.wi_ref = ceil(Nwpsbp/2);

%--Wavelength factors/weights within each sub-bandpass. For full model only
mp.wfs.lambda_weights = ones(Nwpsbp,1); %--Initialize as all ones. Weights within a single sub-bandpass
if(Nwpsbp==1) %--Give equal weighting to all wavelengths
    mp.wfs.dlam = 0; %--Delta lambda between every wavelength in the sub-band in the full model
else %--Give half weighting to the end wavelengths
    %--Spectral weighting in image
    mp.wfs.lambda_weights(1) = 1/2; %--Give end wavelengths half weighting
    mp.wfs.lambda_weights(end) = 1/2; %--Give end wavelengths half weighting
    mp.wfs.fracBWsbp = fracBW/Nsbp; %--Bandwidth per sub-bandpass
    %--Indexing of wavelengths in each sub-bandpass
    sbp_facs = linspace(1-mp.wfs.fracBWsbp/2,1+mp.wfs.fracBWsbp/2, Nwpsbp); %--Factor applied to lambda0 only
    mp.wfs.dlam = (sbp_facs(2) - sbp_facs(1))*lambda0; %--Delta lambda between every wavelength in the full model 
end
mp.wfs.lambda_weights = mp.wfs.lambda_weights/sum(mp.wfs.lambda_weights); %--Normalize sum of the weights (within the sub-bandpass)

%--Make vector of all wavelengths and weights used in the full model
lambdas = zeros(Nsbp*Nwpsbp,1);
lambda_weights_all = zeros(Nsbp*Nwpsbp,1);
mp.wfs.lambdasMat = zeros(Nsbp,Nwpsbp);
mp.wfs.indsLambdaMat = zeros(Nsbp*Nwpsbp,2);
counter = 0;
for si=1:Nsbp
    mp.wfs.lambdasMat(si,:) = (-(Nwpsbp-1)/2:(Nwpsbp-1)/2)*mp.wfs.dlam + mp.wfs.sbp_centers(si); 
    for wi=1:Nwpsbp
        counter = counter+1;
        lambdas(counter) = mp.wfs.lambdasMat(si,wi);
        lambda_weights_all(counter) = mp.wfs.sbp_weights(si)*mp.wfs.lambda_weights(wi);
        mp.wfs.indsLambdaMat(counter,:) = [si,wi];
    end
end

%--Get rid of redundant wavelengths in the complete list, and sum weights for repeated wavelengths
% indices of unique wavelengths
[~, inds_unique] = unique(round(1e12*lambdas)); %--Check equality at the picometer level for wavelength
mp.wfs.indsLambdaUnique = inds_unique;
% indices of duplicate wavelengths
duplicate_inds = setdiff( 1:length(lambdas) , inds_unique);
% duplicate weight values
duplicate_values = lambda_weights_all(duplicate_inds);

%--Shorten the vectors to contain only unique values. Combine weights for repeated wavelengths.
mp.wfs.lambdas = lambdas(inds_unique);
mp.wfs.lambda_weights_all = lambda_weights_all(inds_unique);
for idup=1:length(duplicate_inds)
    lambda = lambdas(duplicate_inds(idup));
    weight = lambda_weights_all(duplicate_inds(idup));
    ind = find(abs(mp.full.lambdas-lambda)<=1e-11);
    mp.wfs.lambda_weights_all(ind) = mp.wfs.lambda_weights_all(ind) + weight;
end
mp.wfs.NlamUnique = length(inds_unique);

end