% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to produce the complex transmission of a focal plane mask made
% of metal and dielectric layered on fused silica.

function FPMmat = falco_gen_HLC_FPM_complex_trans_mat(mp,si,wi,modelType)

%--Different variables for compact and full models
switch modelType
    case 'compact'
        complexTransCube = mp.complexTransCompact; % size = (Ndiel,Nmetal,mp.Nsbp)
        ilam = si; %--Index of the wavelength in mp.sbp_centers
    case 'full'
        complexTransCube = mp.complexTransFull; % size = (Ndiel,Nmetal,mp.Nsbp*mp.Nwpsbp)
        ilam = (si-1)*mp.Nwpsbp + wi; %--Index of the wavelength in mp.lam_array
    otherwise 
        error('falco_gen_HLC_FPM_surf_from_cube: Error! Must specify full or compact model.');
end

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

%--Look up values
FPMmat = zeros(Neta, Nxi); %--Initialize output array of FPM's complex transmission    
for ix = 1:Nxi
    for iy = 1:Neta
        ind_metal = DM8transInd(iy,ix);
        ind_diel = DM9transInd(iy,ix);
        FPMmat(iy,ix) = complexTransCube(ind_diel,ind_metal,ilam);
    end
end

%     %--Coordinates for metal thickness and dielectric thickness
%     [X,Y] = meshgrid(mp.t_metal_nm_vec,mp.t_diel_nm_vec);
%     
%     %--Interpolate values
%     FPMmat = zeros(Neta, Nxi); %--Initialize output array of FPM's complex transmission    
% %     fprintf('Beginning computation of complex transmission for FPM...\t'); tic
%     for ix = 1:Nxi
%             fprintf('Beginning computation of complex transmission for FPM row %d (%d points)...\t',ix,Nxi); tic
%         for iy = 1:Neta
%             FPMmat(iy,ix) = interp2(X, Y, squeeze(complexTransCube(:,:,ilam)), DM8surf(iy,ix), DM9surf(iy,ix),'nearest');
%         end
%         tend = toc;
%         fprintf('done.\tTime = %.2fs\t(%.2fs per pixel)\n',tend,tend/Nxi);
%     end
% %     fprintf('done.\tTime = %.2fs\n',toc);
    

end %--END OF FUNCTION
