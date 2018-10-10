% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to produce the complex transmission of a focal plane mask made
% of metal and dielectric layered on fused silica.

function [FPMcube,DM8surf,DM9surf] = falco_gen_HLC_FPM_complex_trans_cube(mp,modelType)



%--Different variables for compact and full models
switch modelType
    case 'compact'
        complexTransCube = mp.complexTransCompact; % size = (Ndiel,Nmetal,mp.Nsbp)
        lamVec = mp.sbp_centers;
    case 'full'
        complexTransCube = mp.complexTransFull; % size = (Ndiel,Nmetal,mp.Nsbp*mp.Nwpsbp)
        lamVec = mp.full.lambdas;
    otherwise 
        error('falco_gen_HLC_FPM_surf_from_cube: Error! Must specify full or compact model.');
end
Nlam = length(lamVec);

t_diel_bias = mp.t_diel_bias_nm*1e-9; % bias thickness of dielectric [m]


%--Generate thickness profiles of each layer
DM8surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm8,modelType); %--Metal layer profile [m]
DM9surf = t_diel_bias + falco_gen_HLC_FPM_surf_from_cube(mp.dm9,modelType); %--Dielectric layer profile [m]
Nxi = size(DM8surf,2);
Neta = size(DM8surf,1);

if(numel(DM8surf)~=numel(DM9surf))
    error('falco_gen_HLC_FPM_complex_trans_cube: FPM surface arrays must be same size!')
end    


%--Obtain the indices of the nearest thickness values in the complex transmission datacube.
DM8transInd = falco_discretize_FPM_surf(DM8surf, mp.t_metal_nm_vec, mp.dt_metal_nm);
DM9transInd = falco_discretize_FPM_surf(DM9surf, mp.t_diel_nm_vec,  mp.dt_diel_nm);

%--Look up values (WAY faster than interpolating)
FPMcube = zeros(Neta, Nxi, Nlam); %--Initialize output array of FPM's complex transmission    
for ilam=1:Nlam 
    for ix = 1:Nxi
        for iy = 1:Neta
            ind_metal = DM8transInd(iy,ix);
            ind_diel = DM9transInd(iy,ix);
            FPMcube(iy,ix,ilam) = complexTransCube(ind_diel,ind_metal,ilam);
        end
    end
end

%     %--Interpolate values
%     fprintf('Beginning computation of complex transmission for FPM...\t'); tic
%     [X,Y] = meshgrid(mp.t_metal_nm_vec,mp.t_diel_nm_vec); %--Coordinates for metal thickness and dielectric thickness
%     FPMcube = zeros(Neta, Nxi, Nlam); %--Initialize output array of FPM's complex transmission
%     for il=1:Nlam        
%         for ix = 1:Nxi
%             for iy = 1:Neta
%                 FPMcube(iy,ix,il) = interp2(X, Y, squeeze(complexTransCube(:,:,il)), DM8surf(iy,ix), DM9surf(iy,ix));
%             end
%         end
%     end
%     fprintf('done.\tTime = %.2fs\n',toc); 


end %--END OF FUNCTION
