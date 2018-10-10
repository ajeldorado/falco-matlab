% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to produce the complex transmission of a focal plane mask made
% of metal and dielectric layered on fused silica.
%
% Modified on 2018-08-13 by A.J. Riggs for extended PMGI and nickel layers
% (of independent array sizes) and to include an outer opaque iris.
% Created on 2018-05-27 by A.J. Riggs for the HLC occulter.

function FPMmat = falco_gen_EHLC_FPM_complex_trans_mat(mp,si,wi,modelType)

%--Different variables for compact and full models
switch modelType
    case 'compact'
        complexTransCube = mp.complexTransCompact; % size = (Ndiel,Nmetal,mp.Nsbp)
        ilam = si; %--Index of the wavelength in mp.sbp_centers
        Narray = mp.F3.compact.Nxi;
        iris = mp.F3.compact.iris;
        t_diel_nom_m = 1e-9*mp.F3.compact.t_diel_nom_nm;
    case 'full'
        complexTransCube = mp.complexTransFull; % size = (Ndiel,Nmetal,mp.Nsbp*mp.Nwpsbp)
        ilam = (si-1)*mp.Nwpsbp + wi; %--Index of the wavelength in mp.lam_array
        Narray = mp.F3.full.Nxi;
        iris = mp.F3.full.iris;
        t_diel_nom_m = 1e-9*mp.F3.full.t_diel_nom_nm;
    otherwise 
        error('falco_gen_EHLC_FPM_surf_from_cube: Error! Must specify full or compact model.');
end


%--Generate thickness profiles of each layer
DM8surf = padOrCropEven( falco_gen_EHLC_FPM_surf_from_cube(mp.dm8,modelType), Narray); %--Metal layer profile [m]
DM9surf = t_diel_nom_m + padOrCropEven( falco_gen_EHLC_FPM_surf_from_cube(mp.dm9,modelType), Narray); %--Dielectric layer profile [m]

%--Obtain the indices of the nearest thickness values in the complex transmission datacube.
DM8transInd = falco_discretize_FPM_surf(DM8surf, mp.t_metal_nm_vec, mp.dt_metal_nm);
DM9transInd = falco_discretize_FPM_surf(DM9surf, mp.t_diel_nm_vec,  mp.dt_diel_nm);

%--Look up values
FPMmat = zeros(Narray, Narray); %--Initialize output array of FPM's complex transmission    
for ix = 1:Narray
    for iy = 1:Narray
        ind_metal = DM8transInd(iy,ix);
        ind_diel = DM9transInd(iy,ix);
        FPMmat(iy,ix) = complexTransCube(ind_diel,ind_metal,ilam);
    end
end

%--Apply the outer opaque iris
FPMmat = FPMmat.*iris;

% Nxi8 = size(DM8surf,2);
% Neta8 = size(DM8surf,1);
% Nxi9 = size(DM9surf,2);
% Neta9 = size(DM9surf,1);
% if( (Nxi8~=Neta8) || (Nxi9~=Neta9)  )
%     error('falco_gen_EHLC_FPM_complex_trans_cube: FPM arrays must be square.')
% end
% 
% %--Pad up to the larger size (of either DM8surf or DM9surf)
% Nsubarr = max([Nxi8,Nxi9]); 
% DM8surf = padOrCropEven(DM8surf,Nsubarr);
% DM9surf = padOrCropEven(DM9surf,Nsubarr);
% 
% 
% %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
% DM8transInd = falco_discretize_FPM_surf(DM8surf, mp.t_metal_nm_vec, mp.dt_metal_nm);
% DM9transInd = falco_discretize_FPM_surf(DM9surf, mp.t_diel_nm_vec,  mp.dt_diel_nm);
%   
% %--Look up values
% FPMsubarray = zeros(Nsubarr, Nsubarr); %--Initialize output array of FPM's complex transmission    
% for ix = 1:Nsubarr
%     for iy = 1:Nsubarr
%         ind_metal = DM8transInd(iy,ix);
%         ind_diel = DM9transInd(iy,ix);
%         FPMsubarray(iy,ix) = complexTransCube(ind_diel,ind_metal,ilam);
%     end
% end
% 
% 
% %%--OUTER FPM
% 
% %--Complex transmission of the points outside the inner, modulated FPM (just fused silica a bias layer of dielectric).
% ind_metal = falco_discretize_FPM_surf(0, mp.t_metal_nm_vec, mp.dt_metal_nm); %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
% ind_diel = falco_discretize_FPM_surf(t_diel_bias, mp.t_diel_nm_vec,  mp.dt_diel_nm); %--Obtain the indices of the nearest thickness values in the complex transmission datacube.
% transOuterFPM = mp.complexTransCompact(ind_diel,ind_metal,ilam); %--Complex transmission of the points outside the FPM (just fused silica with neither dielectric nor metal).            
% 
% 
% %--Pad to include the outer part of the FPM
% FPMmat = padOrCropEven(FPMsubarray,Narray,'extrapval',transOuterFPM);
% 
% %--Apply the outer opaque iris
% FPMmat = FPMmat.*iris;
% 
% 
% 
%     
% 
% 
% %     %--Coordinates for metal thickness and dielectric thickness
% %     [X,Y] = meshgrid(mp.t_metal_nm_vec,mp.t_diel_nm_vec);
% %     
% %     %--Interpolate values
% %     FPMmat = zeros(Neta, Nxi); %--Initialize output array of FPM's complex transmission    
% % %     fprintf('Beginning computation of complex transmission for FPM...\t'); tic
% %     for ix = 1:Nxi
% %             fprintf('Beginning computation of complex transmission for FPM row %d (%d points)...\t',ix,Nxi); tic
% %         for iy = 1:Neta
% %             FPMmat(iy,ix) = interp2(X, Y, squeeze(complexTransCube(:,:,ilam)), DM8surf(iy,ix), DM9surf(iy,ix),'nearest');
% %         end
% %         tend = toc;
% %         fprintf('done.\tTime = %.2fs\t(%.2fs per pixel)\n',tend,tend/Nxi);
% %     end
% % %     fprintf('done.\tTime = %.2fs\n',toc);
    

end %--END OF FUNCTION
