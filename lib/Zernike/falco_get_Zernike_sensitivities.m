% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run after a FALCO trial to compute the |dE|^2 sensitivities 
% of a coronagraph to Zernike modes.
% 
%
% Modified on 2018-12-11 by A.J. Riggs to be a function.
% Written by A.J. Riggs on 2018-08-10.

function dE2_array = falco_get_Zernike_sensitivities(mp)

%--New fields in structure "mp.eval" that are used:
indsZnoll = mp.eval.indsZnoll;
Rsens = mp.eval.Rsens;

% %% OPTIONAL INPUTS: Keyword flags for "plot" and "parfor"
%   flagPlot  = false;        % flag to show plots or not
%   flagParfor = false;       % flag to use parfor or not
% 
%   icav = 0;             % index in cell array varargin
%   while icav < size(varargin, 2)
%     icav = icav + 1;
%     switch lower(varargin{icav})
%       case {'plot'}
%         flagPlot  = true;   % flag to show plots or not
%       case {'parfor'}
%         flagParfor = 1;     % flag to use parfor or not
%       otherwise
%         error('func_eval_Zernike_sensitivities: Unknown keyword: %s\n', varargin{icav});
%     end
%   end
% 
% %% Updated fonfig file name
% fConfig2 = [fn_config(1:end-4) '_eval.mat'];
% 
% %% Load final FALCO result
% 
% %--Set input parameters
% load(fn_config,'mp');
% 
% %--Load and keep final DM settings from the design trial
% temp = load(fn_snippet,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% if(any(mp.dm_ind==8));  mp.dm8.V = temp.out.DM8V;  end
% if(any(mp.dm_ind==9));  mp.dm9.V = temp.out.DM9V;  end
% 
% clear temp
% 
% %--Flags
% mp.flagPlot = flagPlot; % flag to make plots in falco_init_ws.m
% mp.flagParfor = flagParfor; % flag to use parfor (in eval, only for generating material datacube)
% 
% %% Save new config and initialize the workspace
% %--Save out new config file
% save(fConfig2,'mp');
% 
% %%--Get configuration data from a function file
% [mp,~] = falco_init_ws(fConfig2);
% 

%% Zernikes
Nzern = length(indsZnoll);

rmsZvec = ones(size(indsZnoll))*1e-9;  %--RMS values for each Zernike specified in vector indsZnoll [meters] 

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

% figure(1); imagesc(ZmapCube(:,:,5)); axis xy equal tight; colorbar;

%% Get rid of original tip/tilt offsets (just change wavelength for eval)
% if(isfield(mp,'ttx'))
%     mp = rmfield(mp,'ttx');
% end

%% Make scoring masks

Nannuli = size(Rsens,1);


masks = zeros(size(mp.Fend.corr.mask,1), size(mp.Fend.corr.mask,2), Nannuli);

for ni = 1:Nannuli
    %--Make scoring masks for the annular regions
    maskCorr.pixresFP = mp.Fend.res;
    maskCorr.rhoInner = Rsens(ni,1); %--lambda0/D
    maskCorr.rhoOuter = Rsens(ni,2) ; %--lambda0/D
    maskCorr.angDeg = mp.Fend.corr.ang; %--degrees
    maskCorr.centering = mp.centering;
    maskCorr.FOV = mp.Fend.FOV;
    maskCorr.whichSide = 'both';%mp.Fend.sides; %--which (sides) of the dark hole have open
    %--Compact Model: Generate Software Mask for Correction 
    [masks(:,:,ni), ~, ~] = falco_gen_SW_mask(maskCorr); 
    
%     indsSWmask{ni} = find(masks(:,:,ni)==1);
end

%% Get nominal, unaberrated final E-field at each wavelength.

mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp); %--Input E-fields

E0cube = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.Nsbp); %--initialize


% fprintf('Computing unaberrated E-field at wavelength\t')
for si = 1:mp.Nsbp
%     fprintf('%d/%d  ',si,mp.Nsbp)
    
    modvar.sbpIndex = si;
    modvar.whichSource = 'star';     

    E0cube(:,:,si) = model_compact(mp, modvar);
%     I0c = abs(E0cube(:,:,si)).^2;
%     
%     if(flagPlot)
%         figure(501); imagesc(log10(I0c.*mp.Fend.corr.mask),[-10 -7]); axis xy equal tight; colorbar;
%         title(sprintf('Unaberrated PSF at %d nm',round(mp.sbp_centers(si)*1e9)),'Fontsize',20);
%         set(gca,'FontSize',20); drawnow;
%         pause(0.1)
%     end
end
% fprintf('   done.\n');
        


%% Get aberrated, final focal-plane E-fields for each Zernike mode


EZarray = zeros(mp.Fend.Neta,mp.Fend.Nxi,Nzern,mp.Nsbp); %--initialize
dEZarray = zeros(size(EZarray));

for iz = 1:Nzern
%     fprintf('Computing E-field with Z%d at wavelength\t',indsZnoll(iz))

    %mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nsbp); %--reset to ones

    for si = 1:mp.Nsbp
%         fprintf('%d/%d  ',si,mp.Nsbp)
        
        mp.P1.compact.E(:,:,si) = exp(1i*2*pi/mp.sbp_centers(si)*rmsZvec(iz)*ZmapCube(:,:,iz));%.*mp.P1.compact.E(:,:,si);
        
        modvar.sbpIndex = si;
        modvar.whichSource = 'star';     

        EZarray(:,:,iz,si) = model_compact(mp, modvar);
        dEZarray(:,:,iz,si) = EZarray(:,:,iz,si)-E0cube(:,:,si);
%         IZ = abs(EZarray(:,:,iz,si)).^2;
%         if(flagPlot)
%             figure(490); imagesc(log10(IZ.*mp.Fend.corr.mask),[-10 -7]); axis xy equal tight; colorbar;
%             title(sprintf('PSF with %dnm of Z%d at %d nm',round(1e9*rmsZvec(iz)), indsZnoll(iz),round(mp.sbp_centers(si)*1e9)),'Fontsize',20);
%             set(gca,'FontSize',20); drawnow;
%         end
        
    end
%     fprintf('   done.\n');

end

dE2cube = mean(abs(dEZarray).^2,4); % |dE|^2 averaged over wavelength

%% Compute Zernike sensitivity values averaged across the whole dark hole

dE2_array = zeros(Nzern,Nannuli);
for iz = 1:Nzern
   dEtemp = dE2cube(:,:,iz);
   
   for ni = 1:Nannuli
       dE2_array(iz,ni) = mean(dEtemp( logical(masks(:,:,ni))) );
   end
   
end

%--Print Zernike sensitivity results to command line
% fprintf('\n');
for iz = 1:Nzern
    fprintf('|dE|^2 at %dnm with %dnm RMS of    Z%d =',round(mp.lambda0*1e9),  round(1e9*rmsZvec(iz)), indsZnoll(iz) )
    
    for ni = 1:Nannuli
       mean(dEtemp( logical(masks(:,:,ni))) );
       fprintf('\t%.2e (%.1f-%.1f l/D)',dE2_array(iz,ni), Rsens(ni,1), Rsens(ni,2) )
    end
    fprintf('\n')
end

% %--Plot
% for iz = 1:Nzern
% 
%         if(flagPlot)
%             figure(500+iz); imagesc(log10(dE2cube(:,:,iz))); axis xy equal tight; colorbar;
%             title(sprintf('Band-avg |dE|^2 with %dnm of Z%d at %d nm',round(1e9*rmsZvec(iz)), indsZnoll(iz),round(mp.lambda0*1e9)),'Fontsize',20);
%             drawnow;
%             pause(0.1);
%             
%         end
% end

end %--END OF FUNCTION

