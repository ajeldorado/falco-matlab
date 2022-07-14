function [ev, jacStruct] = get_ekf_jacobian(mp,ev,jacStruct)


% Find values to convert images back to counts rather than normalized
% intensity
ev.peak_psf_counts = zeros(1,mp.Nsbp);
ev.e_scaling = zeros(1,mp.Nsbp);

%     sbp_texp  = tb.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
%     numCoadds = tb.info.sbp_numCoadds(si);% Number of coadds to use for each sub-bandpass
%     PSFpeak   = tb.info.PSFpeaks(si);% counts per second 
%     PSFpeak_counts = PSFpeak*sbp_texp*numCoadds;

for iSubband = 1:mp.Nsbp

    % potentially set mp.detector.peakFluxVec(si) * mp.detector.tExpUnprobedVec(si) set to mp.tb.info.sbp_texp(si)*mp.tb.info.PSFpeaks(si);
    % to have cleaner setup
    ev.peak_psf_counts(iSubband) = mp.tb.info.sbp_texp(iSubband)*mp.tb.info.PSFpeaks(iSubband);
    ev.e_scaling(iSubband) = sqrt(mp.tb.info.PSFpeaks(iSubband));

end

% Get e_scaling
% TODO: need e_scaling

% Rearrange jacobians
jacStruct = rearrange_jacobians(jacStruct,mp);

end

function jacStruct = rearrange_jacobians(jacStruct,mp)

% active_dms = ismember([1:9],mp.dm_ind);

% assume we are max using 2 DMs for estimation
% jacStruct.G_tot = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele*active_dms(1) + mp.dm2.Nele*active_dms(2),Nsbp);

G1 = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele,mp.Nsbp);
G2 = zeros(2*size(jacStruct.G1,1),mp.dm2.Nele,mp.Nsbp);

% Set up jacobian so real and imag components alternate and jacobian from
% each DM is stacked
for iSubband = 1:mp.Nsbp
    
    if any(mp.dm_ind == 1)
        G1_comp = jacStruct.G1(:,:,iSubband);
        G1_split = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele);
        G1_split(1:2:end,:) = real(G1_comp);
        G1_split(2:2:end,:) = imag(G1_comp);

        G1(:,:,iSubband) = G1_split;

    else
        G1 = [];
    end

    if any(mp.dm_ind == 2)
        G2_comp = jacStruct.G2(:,:,iSubband);
        G2_split = zeros(2*size(jacStruct.G2,1),mp.dm2.Nele);
        G2_split(1:2:end,:) = real(G2_comp);
        G2_split(2:2:end,:) = imag(G2_comp);

        G2(:,:,iSubband) = G2_split;

    else
        G2 = [];

    end
end

jacStruct.G_tot = [G1, G2];


end