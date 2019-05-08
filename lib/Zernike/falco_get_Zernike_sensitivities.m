% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to Zernike modes introduced at the entrance pupil.
% 
%
% Modified on 2019-05-02 by A.J. Riggs to use the full model, including the
% option of a PROPER full model.
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2_array = falco_get_Zernike_sensitivities(mp)


rmsVal = 1e-9; %--RMS values for each Zernike specified in vector indsZnoll [meters] 

indsZnoll = mp.eval.indsZnoll;
Rsens = mp.eval.Rsens; %--Radii ranges for the zernike sensitivity calcuations. they are allowed to overlap
Nannuli = size(Rsens,1);
Nzern = length(indsZnoll);

%%--Make scoring masks
masks = zeros(size(mp.Fend.corr.mask,1), size(mp.Fend.corr.mask,2), Nannuli);
for ni = 1:Nannuli
    %--Make scoring masks for the annular regions
    maskStruct.pixresFP = mp.Fend.res;
    maskStruct.rhoInner = Rsens(ni,1); %--lambda0/D
    maskStruct.rhoOuter = Rsens(ni,2) ; %--lambda0/D
    maskStruct.angDeg = mp.Fend.corr.ang; %--degrees
    maskStruct.centering = mp.centering;
    maskStruct.FOV = mp.Fend.FOV;
    maskStruct.whichSide = 'both'; %--which (sides) of the dark hole have open
    %--Generate Software Mask for Correction 
    [masks(:,:,ni), ~, ~] = falco_gen_SW_mask(maskStruct);
end


if(mp.full.flagPROPER)  %--PROPER full models
    
    %--Note: Generate the Zernike modes directly in the PROPER model
    
    %--Get nominal, unaberrated final E-field at each wavelength.
    E0cube = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.Nsbp); %--initialize
    for si = 1:mp.Nsbp   
        modvar.sbpIndex = si;
        modvar.wpsbpIndex = mp.wi_ref; %--for now, just do at center wavelength
        modvar.whichSource = 'star';     
        E0cube(:,:,si) = model_full(mp, modvar);
    end
    
    %--Get aberrated, final focal-plane E-fields for each Zernike mode
    EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,Nzern,mp.Nsbp); %--initialize
    dEZarray = zeros(size(EZarray));
    
    if(isfield(mp.full,'zindex')==false);  mp.full.zindex = [];  mp.full.zval_m = [];  end %--Initialize the Zernike modes to incude as empty if it doesn't exist
    zindex0 = mp.full.zindex; %--Save the original
    zval_m0 = mp.full.zval_m; %--Save the original
    for iz = 1:Nzern
        for si = 1:mp.Nsbp  
            mp.full.zindex = zindex0; %--Reset
             mp.full.zval_m = zval_m0; %--Reset
                       
            %--Set the Zernike index as used in the PROPER full model
            if(any(zindex0==indsZnoll(iz))) %--Add the delta to an existing entry
                zind = find(zindex0==indsZnoll(iz));
                mp.full.zval_m(zind) = mp.full.zval_m(zind) + 1e-9;
                
            else %--Concatenate the Zenike modes to the vector if it isn't included already
                mp.full.zindex = [mp.full.zindex(:),indsZnoll(iz)];
                mp.full.zval_m = [zval_m0(:), 1e-9]; % [meters]
            end
            
            %--Call the full model
            modvar.sbpIndex = si;
            modvar.wpsbpIndex = mp.wi_ref; %--for now, just do at center wavelength
            modvar.whichSource = 'star';
            EZarray(:,:,iz,si) = model_full(mp, modvar);
            dEZarray(:,:,iz,si) = EZarray(:,:,iz,si)-E0cube(:,:,si); %--Compute the delta E-field
        end
    end
    
    
else %--FALCO models

    %%--Generate Zernike map datacube
    ZmapCube = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    %--Make sure ZmapCube is the right array size
    if(size(ZmapCube,1)~=mp.P1.compact.Narr)
        ZmapCubeTemp = zeros(mp.P1.compact.Narr,mp.P1.compact.Narr);
        for zi=1:size(ZmapCube,3)
            ZmapCubeTemp(:,:,zi) = padOrCropEven(ZmapCube(:,:,zi),mp.P1.compact.Narr);
        end
        ZmapCube = ZmapCubeTemp; 
        clear ZmapCubeTemp
    end

    %%--Get nominal, unaberrated final E-field at each wavelength.
    mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp); %--Input E-fields
    E0cube = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.Nsbp); %--initialize
    for si = 1:mp.Nsbp   
        modvar.sbpIndex = si;
        modvar.whichSource = 'star';     
        E0cube(:,:,si) = model_compact(mp, modvar);
    end

    %%--Get aberrated, final focal-plane E-fields for each Zernike mode
    EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,Nzern,mp.Nsbp); %--initialize
    dEZarray = zeros(size(EZarray));
    for iz = 1:Nzern
        for si = 1:mp.Nsbp      
            mp.P1.compact.E(:,:,si) = exp(1i*2*pi/mp.sbp_centers(si)*rmsVal*ZmapCube(:,:,iz));%.*mp.P1.compact.E(:,:,si);       
            modvar.sbpIndex = si;
            modvar.whichSource = 'star';     
            EZarray(:,:,iz,si) = model_compact(mp, modvar);
            dEZarray(:,:,iz,si) = EZarray(:,:,iz,si)-E0cube(:,:,si);
        end
    end
    
    
end


%% Compute Zernike sensitivity values averaged across each annulus (or annular sector) in the dark hole

dE2cube = mean(abs(dEZarray).^2,4); % |dE|^2 averaged over wavelength
dE2_array = zeros(Nzern,Nannuli);
for iz = 1:Nzern
   dEtemp = dE2cube(:,:,iz);

   for ni = 1:Nannuli
       dE2_array(iz,ni) = mean(dEtemp( logical(masks(:,:,ni))) );
   end

end

%--Print Zernike sensitivity results to command line
for iz = 1:Nzern
    fprintf('|dE|^2 at %dnm with %dnm RMS of    Z%d =',round(mp.lambda0*1e9),  round(1e9*rmsVal), indsZnoll(iz) )

    for ni = 1:Nannuli
       mean(dEtemp( logical(masks(:,:,ni))) );
       fprintf('\t%.2e (%.1f-%.1f l/D)',dE2_array(iz,ni), Rsens(ni,1), Rsens(ni,2) )
    end
    fprintf('\n')
end

end %--END OF FUNCTION