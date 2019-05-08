% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to Zernike modes introduced at the entrance pupil.
% 
% Modified on 2019-05-08 by A.J. Riggs to use all wavelengths and to
% parallelize all the E-field calculations when using a PROPER full model.
% Modified on 2019-05-02 by A.J. Riggs to use the full model, including the
% option of a PROPER full model.
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2mat = falco_get_Zernike_sensitivities(mp)

indsZnoll = mp.eval.indsZnoll;
Rsens = mp.eval.Rsens; %--Radii ranges for the zernike sensitivity calcuations. they are allowed to overlap
Nannuli = size(Rsens,1);
Nzern = length(indsZnoll);

%%--Make scoring masks
maskCube = zeros(mp.Fend.Neta,mp.Fend.Nxi, Nannuli); %zeros(size(mp.Fend.corr.mask,1), size(mp.Fend.corr.mask,2), Nannuli);
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
    [maskCube(:,:,ni), ~, ~] = falco_gen_SW_mask(maskStruct);
end


if(mp.full.flagPROPER)  %--Use self-contained full models written in PROPER
    
    %--Zernike modes are made directly in the PROPER model. Not needed here.
    
    %--Number of polarization states used
    mp.full.dummy = 1; %--Initialize if this doesn't exist
    if(isfield(mp.full,'pol_conds'))  
        Npol = length(mp.full.pol_conds);  
    else
        Npol = 1;
    end
    
    %% Get unaberrated E-fields
    %--Loop over all wavelengths and polarizations        
    inds_list = allcomb(1:mp.full.NlamUnique,1:Npol).'; %--dimensions: [2 x mp.full.NlamUnique*Npol ]
    Nvals = size(inds_list,2);
    
    %--Get nominal, unaberrated final E-field at each wavelength and polarization
    E0array = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.full.NlamUnique,Npol); %--initialize
        
    tic; fprintf('Computing unaberrated E-fields...\t');
    %--Obtain all the images in parallel
    if(mp.flagParfor)
        parfor ni=1:Nvals;  Estruct{ni} = falco_get_single_sim_Efield_LamPol(ni,inds_list,mp);  end
    else
        for ni=Nvals:-1:1;  Estruct{ni} = falco_get_single_sim_Efield_LamPol(ni,inds_list,mp);  end
    end
    fprintf('done. Time = %.2f s\n',toc);

    %--Reorganize the output
    for ni=1:Nvals  
        ilam = inds_list(1,ni);
        ipol = inds_list(2,ni);
        E0array(:,:,ilam,ipol) = Estruct{ni};
    end 
    clear Estruct
    
    %% Get E-fields with Zernike aberrations
    %--Loop over all wavelengths, polarizations, and Zernike modes   
    inds_list_zern = allcomb(1:mp.full.NlamUnique,1:Npol,1:Nzern).'; %--dimensions: [3 x mp.full.NlamUnique*Npol*Nzern ]
    NvalsZern = size(inds_list_zern,2);

    %--Get nominal, unaberrated final E-field at each wavelength and polarization
    %EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.full.NlamUnique,Npol,Nzern); %--initialize
    dEZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.full.NlamUnique,Npol,Nzern); %--initialize 
    
    %--Obtain all the images in parallel
    tic; fprintf('Computing aberrated E-fields for Zernike sensitivities...\t');
    if(mp.flagParfor)
        parfor ni=1:NvalsZern;  Estruct{ni} = falco_get_single_sim_Efield_LamPolZern(ni,inds_list_zern,mp);  end
    else
        for ni=NvalsZern:-1:1;  Estruct{ni} = falco_get_single_sim_Efield_LamPolZern(ni,inds_list_zern,mp);  end
    end
    fprintf('done. Time = %.2f s\n',toc);
    
    %--Reorganize the output
    for ni=1:NvalsZern
        ilam = inds_list_zern(1,ni);
        ipol = inds_list_zern(2,ni);
        izern = inds_list_zern(3,ni);
        dEZarray(:,:,ilam,ipol,izern) = Estruct{ni}  - E0array(:,:,ilam,ipol); %--Compute the delta E-field
%         EZarray(:,:,ilam,ipol,izern) = Estruct{ni}; 
%         dEZarray(:,:,ilam,ipol,izern) = EZarray(:,:,ilam,ipol,izern) - E0array(:,:,ilam,ipol); %--Compute the delta E-field
    end
    clear Estruct

% %     %--Get aberrated, final focal-plane E-fields for each Zernike mode
% %     EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,Nzern,Nvals); %--initialize
% %     dEZarray = zeros(size(EZarray));
%     
%     if(isfield(mp.full,'zindex')==false);  mp.full.zindex = [];  mp.full.zval_m = [];  end %--Initialize the Zernike modes to include as empty if the variable doesn't exist already
%     zindex0 = mp.full.zindex; %--Save the original
%     zval_m0 = mp.full.zval_m; %--Save the original
%     for iz = 1:Nzern
%         for ni = 1:mp.Nsbp  
%             mp.full.zindex = zindex0; %--Reset
%             mp.full.zval_m = zval_m0; %--Reset
%                        
%             %--Put the Zernike index and coefficent in the vectors used by the PROPER full model
%             if(any(zindex0==indsZnoll(iz))) %--Add the delta to an existing entry
%                 zind = find(zindex0==indsZnoll(iz));
%                 mp.full.zval_m(zind) = mp.full.zval_m(zind) + 1e-9;
%             else %--Concatenate the Zenike modes to the vector if it isn't included already
%                 mp.full.zindex = [mp.full.zindex(:),indsZnoll(iz)];
%                 mp.full.zval_m = [zval_m0(:), 1e-9]; % [meters]
%             end
%             
%             %--Call the full model
%             ilam = inds_list(1,ni);
%             ipol = inds_list(2,ni);
%             mp.full.polaxis = mp.full.pol_conds(ipol);
%             si = mp.full.indsLambdaMat(ilam,1);
%             wi = mp.full.indsLambdaMat(ilam,2);
%             modvar.sbpIndex = si;
%             modvar.wpsbpIndex = wi;
%             modvar.whichSource = 'star';
%             EZarray(:,:,iz,ni) = model_full(mp, modvar);
%             dEZarray(:,:,iz,ni) = EZarray(:,:,iz,ni)-E0cube(:,:,ni); %--Compute the delta E-field
%         end
%     end
    
    
else %--Use a full model included in FALCO
%--NOT UPDATED YET
%     %%--Generate Zernike map datacube
%     ZmapCube = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
%     %--Make sure ZmapCube is the right array size
%     if(size(ZmapCube,1)~=mp.P1.compact.Narr)
%         ZmapCubeTemp = zeros(mp.P1.compact.Narr,mp.P1.compact.Narr);
%         for zi=1:size(ZmapCube,3)
%             ZmapCubeTemp(:,:,zi) = padOrCropEven(ZmapCube(:,:,zi),mp.P1.compact.Narr);
%         end
%         ZmapCube = ZmapCubeTemp; 
%         clear ZmapCubeTemp
%     end
% 
%     %%--Get nominal, unaberrated final E-field at each wavelength.
%     mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp); %--Input E-fields
%     E0cube = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.Nsbp); %--initialize
%     for si = 1:mp.Nsbp   
%         modvar.sbpIndex = si;
%         modvar.whichSource = 'star';     
%         E0cube(:,:,si) = model_compact(mp, modvar);
%     end
% 
%     %%--Get aberrated, final focal-plane E-fields for each Zernike mode
%     EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,Nzern,mp.Nsbp); %--initialize
%     dEZarray = zeros(size(EZarray));
%     for iz = 1:Nzern
%         for si = 1:mp.Nsbp      
%             mp.P1.compact.E(:,:,si) = exp(1i*2*pi/mp.sbp_centers(si)*mp.full.ZrmsVal*ZmapCube(:,:,iz));%.*mp.P1.compact.E(:,:,si);       
%             modvar.sbpIndex = si;
%             modvar.whichSource = 'star';     
%             EZarray(:,:,iz,si) = model_compact(mp, modvar);
%             dEZarray(:,:,iz,si) = EZarray(:,:,iz,si)-E0cube(:,:,si);
%         end
%     end
    
    
end


%% Compute Zernike sensitivity values averaged across each annulus (or annular sector) in the dark hole

dE2cube = squeeze(mean(mean(abs(dEZarray).^2,4),3)); % |dE|^2 averaged over wavelength and polarization state
dE2mat = zeros(Nzern,Nannuli);
for iz = 1:Nzern
   dEtemp = dE2cube(:,:,iz);
   for ia = 1:Nannuli;  dE2mat(iz,ia) = mean(dEtemp( logical(maskCube(:,:,ia))) );  end
end

%--Print Zernike sensitivity results to command line
for iz = 1:Nzern
    fprintf('|dE|^2 at %dnm with %dnm RMS of    Z%d =',round(mp.lambda0*1e9),  round(1e9*mp.full.ZrmsVal), indsZnoll(iz) )
    for ia = 1:Nannuli
       mean(dEtemp( logical(maskCube(:,:,ia))) );
       fprintf('\t%.2e (%.1f-%.1f l/D)',dE2mat(iz,ia), Rsens(ia,1), Rsens(ia,2) )
    end
    fprintf('\n')
end

end %--END OF FUNCTION


%% Get the stellar E-field for the specified wavelength, polarization, and Zernike aberration
function Estar = falco_get_single_sim_Efield_LamPolZern(ni,inds_list_zern,mp)

ilam  = inds_list_zern(1,ni);
ipol  = inds_list_zern(2,ni);
izern = inds_list_zern(3,ni);

indsZnoll = mp.eval.indsZnoll;
% Rsens = mp.eval.Rsens; %--Radii ranges for the zernike sensitivity calcuations. they are allowed to overlap
% Nannuli = size(Rsens,1);
% Nzern = length(indsZnoll);

if(mp.full.flagPROPER)
    %--Initialize the Zernike modes to include as empty if the variable doesn't exist already
    if(isfield(mp.full,'zindex')==false)  
        mp.full.zindex = [];  
        mp.full.zval_m = [];  
    end 
    zindex0 = mp.full.zindex; %--Save the original
    zval_m0 = mp.full.zval_m; %--Save the original

    %--Put the Zernike index and coefficent in the vectors used by the PROPER full model
    if(any(zindex0==indsZnoll(izern))) %--Add the delta to an existing entry
        zind = find(zindex0==indsZnoll(izern));
        mp.full.zval_m(zind) = mp.full.zval_m(zind) + mp.full.ZrmsVal;
    else %--Concatenate the Zenike modes to the vector if it isn't included already
        mp.full.zindex = [mp.full.zindex(:),indsZnoll(izern)];
        mp.full.zval_m = [zval_m0(:), mp.full.ZrmsVal]; % [meters]
    end
end %--End of mp.full.flagPROPER if statement

%--Get the stellar E-field
modvar.sbpIndex   = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),1);
modvar.wpsbpIndex = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),2);
mp.full.polaxis = mp.full.pol_conds(ipol);
modvar.whichSource = 'star';

Estar = model_full(mp, modvar);
% Iout = (abs(Estar).^2); %--Apply spectral weighting outside this function
    
end %--END OF FUNCTION


%% Get the stellar E-field for the specified wavelength and polarization
function Estar = falco_get_single_sim_Efield_LamPol(ni,inds_list,mp)

ilam = inds_list(1,ni);
ipol = inds_list(2,ni);

%--Get the stellar E-field
modvar.sbpIndex   = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),1);
modvar.wpsbpIndex = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),2);
mp.full.polaxis = mp.full.pol_conds(ipol);
modvar.whichSource = 'star';

Estar = model_full(mp, modvar);
% Iout = (abs(Estar).^2); %--Apply spectral weighting outside this function
    
end %--END OF FUNCTION