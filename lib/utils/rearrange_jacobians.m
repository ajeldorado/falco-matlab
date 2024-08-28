

function G_tot = rearrange_jacobians(mp,jacStruct,dm_inds)

if strcmpi(mp.estimator, 'pairwise-bb') 
    mp.Nsbp = mp.Nsbp_bb; 
end


G1 = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele,mp.Nsbp);
G2 = zeros(2*size(jacStruct.G2,1),mp.dm2.Nele,mp.Nsbp);



% Set up jacobian so real and imag components alternate and jacobian from
% each DM is stacked

for iSubband = 1:mp.Nsbp
    
    if any(dm_inds == 1) 
        G1_comp = jacStruct.G1(:,:,iSubband);
        G1_split = zeros(2*size(jacStruct.G1,1),mp.dm1.Nele);
        G1_split(1:2:end,:) = real(G1_comp);
        G1_split(2:2:end,:) = imag(G1_comp);

        G1(:,:,iSubband) = G1_split;

    else
        G1 = [];
    end

    if any(dm_inds == 2)
        G2_comp = jacStruct.G2(:,:,iSubband);
        G2_split = zeros(2*size(jacStruct.G2,1),mp.dm2.Nele);
        G2_split(1:2:end,:) = real(G2_comp);
        G2_split(2:2:end,:) = imag(G2_comp);

        G2(:,:,iSubband) = G2_split;

    else
        G2 = [];

    end
end

G_tot = [G1, G2];

if strcmpi(mp.estimator, 'pairwise-bb') 
    mp.Nsbp = 1; 
end

end