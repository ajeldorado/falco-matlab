% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to lateral shear at the entrance pupil.
% 
% Modified on 2019-07-09 by A.J. Riggs from
% falco_get_Zernike_sensitivities to falco_get_pupil_shear_sensitivities 
% Modified on 2019-05-08 by A.J. Riggs to use all wavelengths and to
% parallelize all the E-field calculations when using a PROPER full model.
% Modified on 2019-05-02 by A.J. Riggs to use the full model, including the
% option of a PROPER full model.
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2mat = falco_get_pupil_shear_sensitivities(mp)

Rsens = mp.eval.Rsens; %--Radii ranges for the zernike sensitivity calcuations. They are allowed to overlap
Nannuli = size(Rsens,1);
Nshear = 2; %--Number of lateral shear directions (2 for X and Y)

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
    
%% Get E-fields with Pupil Shear 
%--Loop over all wavelengths, polarizations, and pupil shears
inds_list_shear = allcomb(1:mp.full.NlamUnique,1:Npol,1:Nshear).'; %--dimensions: [3 x mp.full.NlamUnique*Npol*Nshear ]
NvalsShear = size(inds_list_shear,2);

%--Get nominal, unaberrated final E-field at each wavelength and polarization
dEZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.full.NlamUnique,Npol,Nshear); %--initialize 

%--Obtain all the images in parallel
tic; fprintf('Computing aberrated E-fields for pupil shear sensitivities...\t');
if(mp.flagParfor)
    parfor ni=1:NvalsShear;  Estruct{ni} = falco_get_single_sim_Efield_LamPolShear(ni,inds_list_shear,mp);  end
else
    for ni=NvalsShear:-1:1;  Estruct{ni} = falco_get_single_sim_Efield_LamPolShear(ni,inds_list_shear,mp);  end
end
fprintf('done. Time = %.2f s\n',toc);
    
    %--Reorganize the output
    for ni=1:NvalsShear
        ilam = inds_list_shear(1,ni);
        ipol = inds_list_shear(2,ni);
        izern = inds_list_shear(3,ni);
        dEZarray(:,:,ilam,ipol,izern) = Estruct{ni}  - E0array(:,:,ilam,ipol); %--Compute the delta E-field
    end
    clear Estruct

%% Compute Zernike sensitivity values averaged across each annulus (or annular sector) in the dark hole

dE2cube = squeeze(mean(mean(abs(dEZarray).^2,4),3)); % |dE|^2 averaged over wavelength and polarization state
dE2mat = zeros(Nshear,Nannuli);
for iz = 1:Nshear
   dEtemp = dE2cube(:,:,iz);
   for ia = 1:Nannuli;  dE2mat(iz,ia) = mean(dEtemp( logical(maskCube(:,:,ia))) );  end
end

%--Print pupil shear sensitivity results to command line
for iz = 1:Nshear
    if(iz==1)
        whichAxis = 'X';
    elseif(iz==2)
        whichAxis = 'Y';
    end
    fprintf('|dE|^2 at %dnm with %.1f microns of    %s pupil shear =',round(mp.lambda0*1e9), mp.full.pupilShearVal*1e6 , whichAxis )
    for ia = 1:Nannuli
       mean(dEtemp( logical(maskCube(:,:,ia))) );
       fprintf('\t%.2e (%.1f-%.1f l/D)',dE2mat(iz,ia), Rsens(ia,1), Rsens(ia,2) )
    end
    fprintf('\n')
end

end %--END OF FUNCTION


%% Get the stellar E-field for the specified wavelength, polarization, and Zernike aberration
function Estar = falco_get_single_sim_Efield_LamPolShear(ni,inds_list_shear,mp)

ilam    = inds_list_shear(1,ni);
ipol    = inds_list_shear(2,ni);
ishear  = inds_list_shear(3,ni);

%--Get the stellar E-field
si = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),1);
wi = mp.full.indsLambdaMat(mp.full.indsLambdaUnique(ilam),2);
modvar.sbpIndex   = si;
modvar.wpsbpIndex = wi;
mp.full.polaxis = mp.full.pol_conds(ipol);
modvar.whichSource = 'star';

if(mp.full.flagPROPER)
    %--Initialize the shear variables if they don't exist
    if(isfield(mp.full,'cgi_x_shift_m')==false);  mp.full.cgi_x_shift_m = 0;  end
    if(isfield(mp.full,'cgi_y_shift_m')==false);  mp.full.cgi_y_shift_m = 0;  end

    %--Apply the extra pupil shear
    if(ishear==1)
        mp.full.cgi_x_shift_m = mp.full.cgi_x_shift_m + mp.full.pupilShearVal; % Apply more X shear
    elseif(ishear==2)
        mp.full.cgi_y_shift_m = mp.full.cgi_y_shift_m + mp.full.pupilShearVal; % Apply more X shear
    end
else
   error('*** Pupil shear option not available outside the PROPER model at the moment. ***'); 
end

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