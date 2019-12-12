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
% Created by A.J. Riggs on 2019-12-10 by extracting material from falco_init_ws.m.
% ---------------

function mp = falco_set_spectral_weights(mp)

%% Compact Model

if(isfield(mp,'Nwpsbp')==false);  mp.Nwpsbp = 1;  end

%--Center(-ish) wavelength indices (ref = reference). (Only the center if
%  an odd number of wavelengths is used.)
mp.si_ref = ceil(mp.Nsbp/2);

%--Wavelengths used for Compact Model (and Jacobian Model)
mp.sbp_weights = ones(mp.Nsbp,1);
if(mp.Nwpsbp==1 && mp.flagSim) %--Set ctrl wavelengths evenly between endpoints (inclusive) of the total bandpass.
    if(mp.Nsbp==1)
        mp.sbp_centers = mp.lambda0;
    else
        mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
        mp.sbp_weights(1) = 1/2; %--Give end sub-bands half weighting
        mp.sbp_weights(end) = 1/2; %--Give end sub-bands half weighting
    end
else %--For cases with multiple sub-bands: Choose wavelengths to be at subbandpass centers since the wavelength samples will span to the full extent of the sub-bands.
    mp.fracBWcent2cent = mp.fracBW*(1-1/mp.Nsbp); %--Bandwidth between centers of endpoint subbandpasses.
    mp.sbp_centers = mp.lambda0*linspace(1-mp.fracBWcent2cent/2,1+mp.fracBWcent2cent/2,mp.Nsbp); %--Space evenly at the centers of the subbandpasses.
end
mp.sbp_weights = mp.sbp_weights/sum(mp.sbp_weights); %--Normalize the sum of the weights

fprintf(' Using %d discrete wavelength(s) in each of %d sub-bandpasses over a %.1f%% total bandpass \n', mp.Nwpsbp, mp.Nsbp,100*mp.fracBW);
fprintf('Sub-bandpasses are centered at wavelengths [nm]:\t '); fprintf('%.2f  ',1e9*mp.sbp_centers); fprintf('\n\n');

%% Bandwidth and Wavelength Specs: Full Model

%--Center(-ish) wavelength indices (ref = reference). (Only the center if an odd number of wavelengths is used.)
mp.wi_ref = ceil(mp.Nwpsbp/2);

%--Wavelength factors/weights within each sub-bandpass. For full model only
mp.full.lambda_weights = ones(mp.Nwpsbp,1); %--Initialize as all ones. Weights within a single sub-bandpass
if(mp.Nwpsbp==1) %--Give equal weighting to all wavelengths
    mp.full.dlam = 0; %--Delta lambda between every wavelength in the sub-band in the full model
else %--Give half weighting to the end wavelengths
    %--Spectral weighting in image
    mp.full.lambda_weights(1) = 1/2; %--Give end wavelengths half weighting
    mp.full.lambda_weights(end) = 1/2; %--Give end wavelengths half weighting
    mp.fracBWsbp = mp.fracBW/mp.Nsbp; %--Bandwidth per sub-bandpass
    %--Indexing of wavelengths in each sub-bandpass
    sbp_facs = linspace(1-mp.fracBWsbp/2,1+mp.fracBWsbp/2,mp.Nwpsbp); %--Factor applied to lambda0 only
    mp.full.dlam = (sbp_facs(2) - sbp_facs(1))*mp.lambda0; %--Delta lambda between every wavelength in the full model 
end
mp.full.lambda_weights = mp.full.lambda_weights/sum(mp.full.lambda_weights); %--Normalize sum of the weights (within the sub-bandpass)

%--Make vector of all wavelengths and weights used in the full model
lambdas = zeros(mp.Nsbp*mp.Nwpsbp,1);
lambda_weights_all = zeros(mp.Nsbp*mp.Nwpsbp,1);
mp.full.lambdasMat = zeros(mp.Nsbp,mp.Nwpsbp);
mp.full.indsLambdaMat = zeros(mp.Nsbp*mp.Nwpsbp,2);
counter = 0;
for si=1:mp.Nsbp
    mp.full.lambdasMat(si,:) = (-(mp.Nwpsbp-1)/2:(mp.Nwpsbp-1)/2)*mp.full.dlam + mp.sbp_centers(si); 
    for wi=1:mp.Nwpsbp
        counter = counter+1;
        lambdas(counter) = mp.full.lambdasMat(si,wi);
        lambda_weights_all(counter) = mp.sbp_weights(si)*mp.full.lambda_weights(wi);
        mp.full.indsLambdaMat(counter,:) = [si,wi];
    end
end

%--Get rid of redundant wavelengths in the complete list, and sum weights for repeated wavelengths
% indices of unique wavelengths
[~, inds_unique] = unique(round(1e12*lambdas)); %--Check equality at the picometer level for wavelength
mp.full.indsLambdaUnique = inds_unique;
% indices of duplicate wavelengths
duplicate_inds = setdiff( 1:length(lambdas) , inds_unique);
% duplicate weight values
duplicate_values = lambda_weights_all(duplicate_inds);

%--Shorten the vectors to contain only unique values. Combine weights for repeated wavelengths.
mp.full.lambdas = lambdas(inds_unique);
mp.full.lambda_weights_all = lambda_weights_all(inds_unique);
for idup=1:length(duplicate_inds)
    lambda = lambdas(duplicate_inds(idup));
    weight = lambda_weights_all(duplicate_inds(idup));
    ind = find(abs(mp.full.lambdas-lambda)<=1e-11);
    mp.full.lambda_weights_all(ind) = mp.full.lambda_weights_all(ind) + weight;
end
mp.full.NlamUnique = length(inds_unique);

%%
% %--Make vector of all wavelengths and weights used in the full model
% mp.full.lambdas = zeros(mp.Nsbp*mp.Nwpsbp,1);
% mp.full.lambda_weights_all = zeros(mp.Nsbp*mp.Nwpsbp,1);
% mp.full.lambdaInds = zeros(mp.Nsbp*mp.Nwpsbp,2); %--Indices as a guide
% counter = 0;
% for si=1:mp.Nsbp
%     for wi=1:mp.Nwpsbp
%         counter = counter+1;
%         mp.full.lambdas(counter) = mp.sbp_centers(si)*mp.full.sbp_facs(wi);
%         mp.full.lambda_weights_all = mp.sbp_weights(si)*mp.full.lambda_weights(wi);
%         mp.full.lambdaInds(counter,:) = [si,wi]; %--Indices
%     end
% end
% mp.full.Nlambdas = mp.Nsbp*mp.Nwpsbp;


end