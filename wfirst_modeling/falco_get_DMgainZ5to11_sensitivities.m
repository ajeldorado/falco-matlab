% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to DM gain error when fitting Zernike modes 5 to 11.
% 
% Modified on 2019-07-23 by A.J. Riggs to compute the DM gain error from
% Zernikes rather than the Zernikes themselves.
% Modified on 2019-05-08 by A.J. Riggs to use all wavelengths and to
% parallelize all the E-field calculations when using a PROPER full model.
% Modified on 2019-05-02 by A.J. Riggs to use the full model, including the
% option of a PROPER full model.
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2mat = falco_get_DMgainZ5to11_sensitivities(mp)

Rsens = mp.eval.Rsens; %--Radii ranges for the zernike sensitivity calcuations. They are allowed to overlap
Nannuli = size(Rsens,1);
indsZnoll = 5:11;
Nzern = length(indsZnoll);

% dm1V0 = mp.dm1.V; %--Store unperturbed DM1 command for resetting.
scaleFac = 1e-9/mean(mp.dm1.VtoH(:)); %--Change to 1nm/V gain. DM maps computed for (10% of) 1V RMS, giving 1nm RMS.

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

%%--Generate Zernike map datacube at DM actuator resolution
dVzernCube = falco_gen_norm_zernike_maps(46.3,'interpixel',indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
%--Make sure ZmapCube is padded or cropped to the right array size
if(size(dVzernCube,1)~=mp.dm1.Nact)
    ZmapCubeTemp = zeros(mp.dm1.Nact,mp.dm1.Nact);
    for zi=1:size(dVzernCube,3)
        ZmapCubeTemp(:,:,zi) = padOrCropEven(dVzernCube(:,:,zi),mp.dm1.Nact);
    end
    dVzernCube = ZmapCubeTemp; 
    clear ZmapCubeTemp
end    
dVzernCube = scaleFac*dVzernCube; %--Adjust voltage to get 1nm RMS of each Zernike.   

%%--Gain error maps: normal distribution, sigma=0.06, mean = 0;
Nrand = 10; %--Number of random 2-D maps
sigmaVal = 0.06;
gainErrorCube = sigmaVal*randn(mp.dm1.Nact,mp.dm1.Nact,Nrand);
gainErrorCube(gainErrorCube<-sigmaVal) = -sigmaVal; %--Clip the distribution
gainErrorCube(gainErrorCube>sigmaVal) = sigmaVal; %--Clip the distribution

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
    
%% Get E-fields with DM gain errors from fitting Zernike modes
%--Loop over all wavelengths, polarizations, and Zernike modes   
inds_list_zern = allcomb(1:mp.full.NlamUnique,1:Npol,1:Nzern,1:Nrand).'; %--dimensions: [4 x mp.full.NlamUnique*Npol*Nzern*Nrand ]
NvalsAb = size(inds_list_zern,2);

%--Get nominal, unaberrated final E-field at each wavelength and polarization
dEZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.full.NlamUnique,Npol,Nzern,Nrand); %--initialize 

%--Obtain all the images in parallel
tic; fprintf('Computing aberrated E-fields for DM gain sensitivities...\t');
if(mp.flagParfor)
    parfor ni=1:NvalsAb;  Estruct{ni} = falco_get_single_sim_Efield_LamPolZernGain(ni,inds_list_zern,dVzernCube,gainErrorCube,mp);  end
else
    for ni=NvalsAb:-1:1;  Estruct{ni} = falco_get_single_sim_Efield_LamPolZernGain(ni,inds_list_zern,dVzernCube,gainErrorCube,mp);  end
end
fprintf('done. Time = %.2f s\n',toc);
    
    %--Reorganize the output
    for ni=1:NvalsAb
        ilam = inds_list_zern(1,ni);
        ipol = inds_list_zern(2,ni);
        izern = inds_list_zern(3,ni);
        igain = inds_list_zern(4,ni);
        dEZarray(:,:,ilam,ipol,izern,igain) = Estruct{ni}  - E0array(:,:,ilam,ipol); %--Compute the delta E-field
    end
    clear Estruct

%% Compute Zernike sensitivity values averaged across each annulus (or annular sector) in the dark hole

dE2cube = squeeze(mean(mean( mean(abs(dEZarray).^2,6) ,4),3)); % |dE|^2 averaged over gain map, then polarization state, and then wavelength.
dE2mat = zeros(Nzern,Nannuli);
for iz = 1:Nzern
   dEtemp = dE2cube(:,:,iz);
   for ia = 1:Nannuli;  dE2mat(iz,ia) = mean(dEtemp( logical(maskCube(:,:,ia))) );  end
end

%--Print Zernike sensitivity results to command line
for iz = 1:Nzern
    fprintf('|dE|^2 at %dnm from 10%% DM gain error on 1nm RMS of    Z%d =',round(mp.lambda0*1e9), indsZnoll(iz) )
    for ia = 1:Nannuli
       mean(dEtemp( logical(maskCube(:,:,ia))) );
       fprintf('\t%.2e (%.1f-%.1f l/D)',dE2mat(iz,ia), Rsens(ia,1), Rsens(ia,2) )
    end
    fprintf('\n')
end

end %--END OF FUNCTION


%% Get the stellar E-field for the specified wavelength, polarization, and DM gain aberration for a given Zernike mode
function Estar = falco_get_single_sim_Efield_LamPolZernGain(ni,inds_list_zern,dVzernCube,gainErrorCube,mp)

ilam  = inds_list_zern(1,ni);
ipol  = inds_list_zern(2,ni);
izern = inds_list_zern(3,ni);
igain = inds_list_zern(4,ni);

%--Error in DM gain, tested as a delta voltage applied to DM1:
dDM1V = dVzernCube(:,:,izern).*gainErrorCube(:,:,igain);
mp.dm1.V = mp.dm1.V + dDM1V;

%--Get the stellar E-field
si = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),1);
wi = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),2);
modvar.sbpIndex   = si;
modvar.wpsbpIndex = wi;
mp.full.polaxis = mp.full.pol_conds(ipol);
modvar.whichSource = 'star';

Estar = model_full(mp, modvar);
    
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