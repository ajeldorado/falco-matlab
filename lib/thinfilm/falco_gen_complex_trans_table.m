% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Calculate thin-film complex transmission data cube. The three dimensions
% are for metal thickness, dielectric thickness, and wavelength.
%
% function [complexTrans] = falco_gen_complex_trans_table(mp,varargin)%
%
% REQUIRED INPUT:
% -mp = structure of model parameters
%
% OPTIONAL INPUT:
% -'reflection','refl','lowfs','lowfsc' = keyword to return the complex 
% reflection coefficient instead of the complex transmission coefficient.
%
% OUTPUTS:
% -complexTransCompact = complex field reflection coefficients sampled within the
%   possible thicknesses of metal and dielectric and at the chosen
%   wavelengths. For compact model.
% -complexTransFull = same as above, but for full model.
%
% REVISION HISTORY:
% -------------------------------------------------------------------------
% Modified on 2019-01-28 by A.J. Riggs to use pol=2 for the mean of the
% complex transmission for the two polarizations.
% Modified on 2018-05-01 by A.J. Riggs.
% Created on 2018-03 by Erkin Sidick.
% -------------------------------------------------------------------------

function [complexTransCompact, complexTransFull] = falco_gen_complex_trans_table(mp,varargin)

% Set default values of input parameters
flagRefl = false; % flag to take the value in reflection instead of transmission
substrate = 'FS'; % material name of the mask substrate [FS or N-BK7]
%--Enable different arguments values by using varargin
icav = 0; % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'reflection','refl','lowfs','lowfsc'}
        flagRefl = true; % flag to take the value in reflection instead of transmission
      case {'substrate','sub','glass'}
        icav = icav + 1;
        substrate = varargin{icav}; % material name of the mask substrate [FS or N-BK7]
      otherwise
        error('falco_gen_complex_trans_table: Unknown keyword: %s\n', varargin{icav});
    end
end

%%

mp.F3.metal = 'Ni';
mp.F3.diel = 'PMGI';

fn_cube_compact = sprintf('%s/data/material/ct_cube_%s_Ti%.1fnm_%s_%.1fto%.1fby%.2f_%s_%.1fto%.1fby%.2f_wvl%dnm_BW%.1fN%d_%.1fdeg_compact.mat',...
    mp.path.falco,substrate,mp.t_Ti_nm,mp.F3.metal,min(mp.t_metal_nm_vec), max(mp.t_metal_nm_vec), mp.dt_metal_nm, ...
    mp.F3.diel, min(mp.t_diel_nm_vec),  max(mp.t_diel_nm_vec),  mp.dt_diel_nm, ...
    (1e9*mp.lambda0),100*mp.fracBW,mp.Nsbp,mp.aoi);

fn_cube_full = sprintf('%s/data/material/ct_cube_%s_Ti%.1fnm_%s_%.1fto%.1fby%.2f_%s_%.1fto%.1fby%.2f_wvl%dnm_BW%.1f_%dN%d_%.1fdeg_full.mat',...
    mp.path.falco,substrate,mp.t_Ti_nm,mp.F3.metal,min(mp.t_metal_nm_vec), max(mp.t_metal_nm_vec), mp.dt_metal_nm, ...
    mp.F3.diel, min(mp.t_diel_nm_vec),  max(mp.t_diel_nm_vec),  mp.dt_diel_nm, ...
    (1e9*mp.lambda0),100*mp.fracBW,mp.Nsbp,mp.Nwpsbp,mp.aoi);

if(flagRefl)
    fn_cube_compact = [fn_cube_compact(1:end-4),'_refl.mat'];
    fn_cube_full = [fn_cube_full(1:end-4),'_refl.mat'];
end

t_Ti_m = 1e-9*mp.t_Ti_nm; %--Static base layer of titanium beneath any nickel.
aoi = mp.aoi;
d0fac = mp.FPM.d0fac;
Nsbp = mp.Nsbp;
t_diel_m_vec = 1e-9*mp.t_diel_nm_vec; %--PMGI thickness range
t_metal_m_vec = 1e-9*mp.t_metal_nm_vec; %--nickel thickness range

Nmetal = length(mp.t_metal_nm_vec);
Ndiel  = length(mp.t_diel_nm_vec);

%--Compact Model: Load the data if it has been generated before; otherwise generate it.
if(exist(fn_cube_compact,'file')==2)
    load(fn_cube_compact,'complexTransCompact');
    fprintf('Loaded complex transmission datacube: %s\n', fn_cube_compact)    
else

    fprintf('Computing thin film equations for compact model:\n')
    complexTransCompact = zeros(Ndiel,Nmetal,mp.Nsbp);
    sbp_centers = mp.sbp_centers;
    
    %--Parallel/distributed computing
    if(mp.flagParfor) 
        parfor si = 1:Nsbp
            lam = sbp_centers(si);
            d0 = lam * d0fac; % Max thickness of PMGI + Ni
            for imetal = 1:Nmetal        
                [tCoef,rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_m, t_metal_m_vec(imetal), t_diel_m_vec, d0, 2);    
                if(flagRefl)
                    complexTransCompact(:,imetal,si) = rCoef;  
                else
                    complexTransCompact(:,imetal,si) = tCoef;
                end
            end
            fprintf('\tDone computing wavelength %d of %d.\n',si,Nsbp);
        end
    %--Regular (serial) computing
    else 
        for si = 1:Nsbp
            lam = sbp_centers(si);
            d0 = lam * mp.FPM.d0fac; % Max thickness of PMGI + Ni
            
            [tCoef,rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_m, t_metal_m_vec, t_diel_m_vec, d0, 2); 
            if(flagRefl)
                complexTransCompact(:,:,si) = rCoef;       
            else
                complexTransCompact(:,:,si) = tCoef;
            end
            fprintf('\tDone computing wavelength %d of %d.\n',si,Nsbp);
        end        
    end
    
    %--Save out for future use
    save(fn_cube_compact,'complexTransCompact')
    fprintf('Saved complex transmission datacube: %s\n', fn_cube_compact)

end

%--Full Model: Load the data if it has been generated before; otherwise generate it.
if(exist(fn_cube_full,'file')==2)
    load(fn_cube_full,'complexTransFull');
    fprintf('Loaded complex transmission datacube: %s\n', fn_cube_full)    
else
    
    fprintf('Computing thin film equations for full model:\n')
    if(mp.Nwpsbp==1)
        complexTransFull = complexTransCompact;
    else
        complexTransFull = zeros(Ndiel,Nmetal,mp.Nsbp*mp.Nwpsbp);
        lambdas = mp.full.lambdas;
        
        %--Parallel/distributed computing
        if(mp.flagParfor) 
            parfor li = 1:length(lambdas)
                lam = lambdas(li);
                d0 = lam * d0fac; % Max thickness of PMGI + Ni
                
                for imetal = 1:Nmetal
                    [tCoef,rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_m, t_metal_m_vec(imetal), t_diel_m_vec, d0, 2);            
                    if(flagRefl)
                        complexTransFull(:,imetal,li) = tCoef;
                    else
                        complexTransFull(:,imetal,li) = rCoef;
                    end
                end
                fprintf('\tDone computing wavelength %d of %d.\n',li,length(lambdas));
            end
        else %--Regular (serial) computing
            for li = 1:length(lambdas)
                lam = lambdas(li);
                d0 = lam * mp.FPM.d0fac; % Max thickness of PMGI + Ni
                
                [tCoef,rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_m, t_metal_m_vec, t_diel_m_vec, d0, 2); 
                if(flagRefl)
                    complexTransFull(:,:,li) = rCoef;
                else
                    complexTransFull(:,:,li) = tCoef;
                end
                fprintf('\tDone computing wavelength %d of %d.\n',li,length(lambdas));
            end
        end
    end
    
    %--Save out for future use
    save(fn_cube_full,'complexTransFull')
    fprintf('Saved complex transmission datacube: %s\n', fn_cube_full)

end % END OF FUNCTION