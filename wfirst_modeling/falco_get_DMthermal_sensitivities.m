% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to 1mK DM thermal drift.
% 
% Modified on 2019-07-23 by A.J. Riggs to compute the DM gain error from
% DM thermal drift.
% Modified on 2019-05-08 by A.J. Riggs to use all wavelengths and to
% parallelize all the E-field calculations when using a PROPER full model.
% Modified on 2019-05-02 by A.J. Riggs to use the full model, including the
% option of a PROPER full model.
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2vec = falco_get_DMthermal_sensitivities(mp)

Rsens = mp.eval.Rsens; %--Radii ranges for the sensitivity calcuations. They are allowed to overlap
Nannuli = size(Rsens,1);

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
    
%% Get E-fields with DM thermal drift applied

%--Add in the bias voltage (50V assuming 5nm/V, so 250V assuming 1nm/V)
bias = 250;
mp.dm1.V = bias + mp.dm1.V;
mp.dm2.V = bias + mp.dm2.V;

%--Add the dV from 1mK temperature drift assuming a rate of 2.6 percent/Kelvin.
mp.dm1.V = (1+2.6/100/1000)*mp.dm1.V - bias; %--Subtract bias to avoid numerical issues
mp.dm2.V = (1+2.6/100/1000)*mp.dm2.V - bias; %--Subtract bias to avoid numerical issues

%--Loop over all wavelengths and polarizations        
inds_list = allcomb(1:mp.full.NlamUnique,1:Npol).'; %--dimensions: [2 x mp.full.NlamUnique*Npol ]
Nvals = size(inds_list,2);

%--Get nominal, unaberrated final E-field at each wavelength and polarization
E1array = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.full.NlamUnique,Npol); %--initialize

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
    E1array(:,:,ilam,ipol) = Estruct{ni};
end 
clear Estruct

dEZarray = E1array - E0array;

%% Compute DM thermal sensitivity values averaged across each annulus (or annular sector) in the dark hole

dE2mat = squeeze(mean(mean(abs(dEZarray).^2,4),3)); % |dE|^2 averaged over wavelength and polarization state
dE2vec = zeros(Nannuli,1);

for ia = 1:Nannuli;  dE2vec(ia) = mean(dE2mat( logical(maskCube(:,:,ia))) );  end

%--Print sensitivity results to command line
fprintf('|dE|^2 from 1 mK thermal drift =')
for ia = 1:Nannuli
   fprintf('\t%.2e (%.1f-%.1f l/D)',dE2vec(ia), Rsens(ia,1), Rsens(ia,2) )
end
fprintf('\n')

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