% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Apply spatial weighting to the Jacobian.
%
% INPUTS
% ------
% mp : structure of model parameters
% jacStruct : structure containing the Jacobians
%
% OUTPUTS
% -------
% jacStruct : structure containing the Jacobians

function jacStruct = falco_apply_spatial_weighting_to_Jacobian(mp, jacStruct)
    
    for iStar = 1:mp.compact.star.count
        
        if any(mp.dm_ind == 1)
            jacStruct.G1(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G1(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm1.Nele, mp.jac.NmodePerStar]);
        end
        
        if any(mp.dm_ind == 2)
            jacStruct.G2(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G2(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm2.Nele, mp.jac.NmodePerStar]);
        end
        
        if any(mp.dm_ind == 8)
            jacStruct.G8(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G8(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm8.Nele, mp.jac.NmodePerStar]);
        end
        
        if any(mp.dm_ind == 9)
            jacStruct.G9(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G9(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm9.Nele, mp.jac.NmodePerStar]);
        end
        
    end
    
end