% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to produce the complex transmission of a focal plane mask made
% of metal and dielectric layered on fused silica.

function [FPMcube,DM8surf,DM9surf] = falco_gen_EHLC_FPM_complex_trans_cube(mp,modelType)

%--Different variables for compact and full models
switch modelType
    case 'compact'
        complexTransCube = mp.complexTransCompact; % size = (Ndiel,Nmetal,mp.Nsbp)
        lamVec = mp.sbp_centers;
        Narray = mp.F3.compact.Nxi;
        iris = mp.F3.compact.iris;
        t_diel_nom_m = 1e-9*mp.F3.compact.t_diel_nom_nm;
    case 'full'
        complexTransCube = mp.complexTransFull; % size = (Ndiel,Nmetal,mp.Nsbp*mp.Nwpsbp)
        lamVec = mp.full.lambdas;
        Narray = mp.F3.full.Nxi;
        iris = mp.F3.full.iris;
        t_diel_nom_m = 1e-9*mp.F3.full.t_diel_nom_nm;
    otherwise 
        error('falco_gen_EHLC_FPM_surf_from_cube: Error! Must specify full or compact model.');
end
Nlam = length(lamVec);

%--Generate thickness profiles of each layer
DM8surf = padOrCropEven( falco_gen_EHLC_FPM_surf_from_cube(mp.dm8,modelType), Narray); %--Metal layer profile [m]
DM9surf = t_diel_nom_m + padOrCropEven( falco_gen_EHLC_FPM_surf_from_cube(mp.dm9,modelType), Narray); %--Dielectric layer profile [m]

%--Obtain the indices of the nearest thickness values in the complex transmission datacube.
DM8transInd = falco_discretize_FPM_surf(DM8surf, mp.t_metal_nm_vec, mp.dt_metal_nm);
DM9transInd = falco_discretize_FPM_surf(DM9surf, mp.t_diel_nm_vec,  mp.dt_diel_nm);
 
%--Look up values
FPMcube = zeros(Narray, Narray, Nlam); %--Initialize output array of FPM's complex transmission    
for ilam=1:Nlam 
    for ix = 1:Narray
        for iy = 1:Narray
            ind_metal = DM8transInd(iy,ix);
            ind_diel = DM9transInd(iy,ix);
            FPMcube(iy,ix,ilam) = complexTransCube(ind_diel,ind_metal,ilam);
        end
    end
end

%--Apply the outer opaque iris
for ilam=1:Nlam     
    FPMcube(:,:,ilam) = iris.*FPMcube(:,:,ilam);    
end

end %--END OF FUNCTION