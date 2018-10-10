% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to estimate the final focal plane electric field with 
%  pair-wise probing and batch process estimation.
%
%--INPUTS
% x_in: Scalar value.
%
%--OUTPUTS
%  x_out: even-valued integer value
%
%--VERSION CHANGE HISTORY
% Modified on 2018-04-23 by A.J. Riggs from the Princeton HCIL lab code.
% Written on 2015-02-19 by A.J. Riggs.
%


function [estvar,ip] = falco_est_batch(DM,mp,ep,ip,modvar,estvar)

% Definitions:
if(isfield(ep,'probeDM')==2)
    Nact = mp.dm2.Nact;
else
    Nact = mp.dm1.Nact;
end

DM1Vnom = mp.dm1.V;
DM2Vnom = mp.dm2.V;
wi_ref = mp.wi_ref;
Nsbp = mp.Nsbp;
sbp_centers = mp.sbp_centers;
corr_inds = mp.F4.corr.inds; %--IMPORTANT: Images are from full model but control is on the compact model, so resolutions at F4 must be same for compact and full models.
score_inds = mp.F4.score.inds; %--IMPORTANT: Images are from full model but control is on the compact model, so resolutions at F4 must be same for compact and full models.

% Itr = estvar.Itr;
Npairs = ep.Npairs; % % Number of image PAIRS for DM Diversity or Kalman filter initialization
IcubeOut = zeros(mp.F4.Neta,mp.F4.Nxi,1+2*Npairs);
VcubeOut = zeros(mp.F4.Neta,mp.F4.Nxi,2,1+2*Npairs);

% Nprobes=Npairs*2;   
probePhaseVec = [0 Npairs];
for k = 1:Npairs-1
    probePhaseVec = [probePhaseVec probePhaseVec(end)-(Npairs-1)];% % #ok<AGstroke2ROW>
    probePhaseVec = [probePhaseVec probePhaseVec(end)+(Npairs)]; % #ok<AGstroke2ROW>
end
probePhaseVec = probePhaseVec*pi/(Npairs);
% ProbePhsVec = (offset + ProbePhsVec)*pi/(Npairs);

% ImNum0 = ip.ImNum;

%% Initialize output arrays
Eest_array = zeros(length(corr_inds),Nsbp);
Iinco_array = Eest_array;

fprintf('Estimating electric field with batch process estimation ...\n'); tic;

IunprobedSum = 0;
for si=1:Nsbp
fprintf('Wavelength: %u/%u ... ',si,Nsbp);

% % Valid for all calls to model_compact.m:
modvar.flagGenMat = 0; 
modvar.sbpIndex = si;
modvar.wpsbpIndex = wi_ref;
modvar.whichSource = 'star';

    
%% Measure current contrast level average, and on each side of Image Plane
% Reset DM commands to the unprobed state:
mp.dm1.V = DM1Vnom;
mp.dm2.V = DM2Vnom;

%% Separate out values of images at dark hole pixels and delta DM voltage settings
Iplus  = zeros( [length(corr_inds), Npairs]); % Pixels of plus probes' intensities
Iminus = zeros( [length(corr_inds), Npairs]); % Pixels of minus probes' intensities
DM1Vplus  = zeros( [Nact,Nact, Npairs]);
DM1Vminus = zeros( [Nact,Nact, Npairs]);
DM2Vplus  = zeros( [Nact,Nact, Npairs]);
DM2Vminus = zeros( [Nact,Nact, Npairs]);

%% Compute probe shapes and take probed images:


if( (mp.flagSim==true) ) %--For simulations
    
    %--Take first image (for unprobed DM settings).
    modvar.whichImg = 1;
%     [Iunprobed,ip] = hcil_getImage(DM_config, ip, mp, modvar);
    [Iunprobed,ip] = falco_get_image(mp, modvar, DM,ip);
    IupVec = Iunprobed(corr_inds); % Pixels of unprobed image in dark hole
    % ip.num_im = num_im0;
    IunprobedSum = IunprobedSum+Iunprobed; %--What is this?
    
    IcubeOut(:,:, modvar.whichImg) = Iunprobed;
    VcubeOut(:,:,1, modvar.whichImg) = mp.dm1.V;
    VcubeOut(:,:,2, modvar.whichImg) = mp.dm2.V;

    ep.InormMean(estvar.Itr) = mean(Iunprobed(score_inds));
    ep.InormCorrMean(estvar.Itr) = mean(Iunprobed(corr_inds));
    fprintf('Measured unprobed Inorm (Corr / Score): %.2e \t%.2e \n',ep.InormCorrMean(estvar.Itr),ep.InormMean(estvar.Itr));    

    % Set (approximate) probe intensity based on current measured Inorm
    ep.InormProbeMax = 1e-4;
    InormProbe = min( [sqrt(ep.InormMean(estvar.Itr)*1e-5), ep.InormProbeMax]);
    fprintf('Chosen probe intensity: %.2e \n',InormProbe);    
    
    
    iodd=1; ieven=1; sumIprobed = 0;
    for qq=1:2*Npairs
        
        dx = mp.Ddm/2/mp.Ndm;
        xs = (-mp.Ndm:mp.Ndm-1)'*dx + dx/2;
        XS = repmat(xs.',2*mp.Ndm,1);
        YS = XS.';            

        ProbeH=hcil_makeProbeSurf(ep.ProbeArea,mp.Ddm,...
            sbp_centers(si),probePhaseVec(qq),ep.ProbeoffsetX,ep.ProbeoffsetY,XS,YS,InormProbe); %Function used to compute the Probe Shape. Output is in meters.
        ProbeH = 1*ProbeH; % Scale the probe amplitude empirically if needed
        dDMVprobe = imresize(ProbeH,Nact/length(XS));  % crude way for now, in units of meters
                  
        % Select which DM to use for probing. (DM2 works better)
        if(ep.whichDMprobe == 1)
            dDM1Vprobe = dDMVprobe./mp.dm1.VtoH; % Now in volts
            dDM2Vprobe = 0;
        elseif(ep.whichDMprobe == 2)
            dDM1Vprobe = 0;        
            dDM2Vprobe = dDMVprobe./mp.dm1.VtoH; % Now in volts
        end
%         dDM2Vprobe = dDM2Vprobe.*padarray(ones(32),[5,5],0);
        % cd(mp.folders.m); fn = sprintf('map_probe%d_DM1.fits',qq); fitswrite(dDM1Vprobe,fn);
        mp.dm1.V = DM1Vnom+dDM1Vprobe;
        mp.dm2.V = DM2Vnom+dDM2Vprobe;

        modvar.whichImg = 1+qq;
        %[Im,ip] = hcil_getImage(DM_config, ip, mp.folders, mp, modvar);
        [Im,ip] = falco_get_image(mp, modvar, DM, ip);
        IcubeOut(:,:,modvar.whichImg) = Im;
        VcubeOut(:,:,1,modvar.whichImg) = mp.dm1.V;
        VcubeOut(:,:,2,modvar.whichImg) = mp.dm2.V;
        probeSign = ['-','+'];
        fprintf('Actual Probe %d%s Contrast is: %.2e \n',ceil(qq/2),probeSign(mod(qq,2)+1),mean(Im(corr_inds)));
        sumIprobed = sumIprobed + mean(Im(corr_inds));

        % Store values:
        if mod(qq,2)==1  % Odd; for plus probes
            DM1Vplus(:,:,iodd) = dDM1Vprobe + DM1Vnom;
            DM2Vplus(:,:,iodd) = dDM2Vprobe + DM2Vnom;
            Iplus(:,iodd) = Im(corr_inds);
            iodd=iodd+1;
        elseif mod(qq,2)==0  % Even; for minus probes
            DM1Vminus(:,:,ieven) = dDM1Vprobe + DM1Vnom; 
            DM2Vminus(:,:,ieven) = dDM2Vprobe + DM2Vnom; 
            Iminus(:,ieven) = Im(corr_inds);
            ieven=ieven+1;      
        end
    end
    meanIprobed = sumIprobed/(2*Npairs); % Mean probe contrast

end

% figure(24); imagesc(log10(abs(Im))); axis xy equal tight; colorbar;
% pause(1);

% %% Use model-based or measurement-based probes, or both
% Eest = zeros(length(corr_inds),1);
% if( strcmpi(ep.probeAmpType,'modelJac') )
%         %%G2 = zeros(length(corr_inds),mp.Nact^2,Nsbp)
%         dDM2VprobePlus = DM2Vplus-repmat(DM2Vnom,[1,1,Npairs]);
%         JacDM2 = ep.G2(:,:,si);
%         HmodelAll = zeros(Npairs,2,length(corr_inds)); % Observation matrix
%         for qq=1:Npairs
%             dDM2VprobeVec = reshape(dDM2VprobePlus(:,:,qq),[mp.Nact^2,1]);
%             dElin = JacDM2*dDM2VprobeVec; % dimensions: [Npix x 1]
%             for jj=1:length(corr_inds)
%                 HmodelAll(qq,:,jj) =  [real(dElin(jj)); imag(dElin(jj))];
%             end
%         end
%         if( strcmpi(ep.probeAmpType,'modelJac') )
%             for jj=1:length(corr_inds)
%                 Epix = pinv(HmodelAll(:,:,jj))*z(:,jj);
%                 Eest(jj) = Epix(1) + 1i*Epix(2); 
%             end
%         end
%     
% elseif( strcmpi(ep.probeAmpType,'meas') || strcmpi(ep.probeAmpType,'both') )

%% Propagate each DM setting to the image plane. (model based)
% For unprobed field based on model:
mp.dm1.V = DM1Vnom;
mp.dm2.V = DM2Vnom; % Added July 9, 2014
Eunprobed = hcil_model(mp, DM_config,modvar);
EupVec = Eunprobed(corr_inds);

% For probed fields based on model:
Eplus  = zeros(size(Iplus ));
Eminus = zeros(size(Iminus));
for qq=1:Npairs
    % For plus probes:
    mp.dm1.V = squeeze( DM1Vplus(:,:,qq));
    mp.dm2.V = squeeze( DM2Vplus(:,:,qq));
    Etemp = hcil_model(mp, DM_config, modvar);
    Eplus(:,qq) = Etemp(corr_inds);
    % For minus probes:
    mp.dm1.V = squeeze( DM1Vminus(:,:,qq));
    mp.dm2.V = squeeze( DM2Vminus(:,:,qq));
    Etemp = hcil_model(mp, DM_config, modvar);
    Eminus(:,qq) = Etemp(corr_inds);
end

%% Create delta E-fields for each probe image. Then create Npairs phase angles.
dEplus  = Eplus -repmat(EupVec,[1,Npairs]);
dEminus = Eminus-repmat(EupVec,[1,Npairs]);
dphdm  = zeros( [length(corr_inds), Npairs]); % phase angles
for qq=1:Npairs
    dphdm(:,qq) = atan2( imag(dEplus(:,qq))-imag(dEminus(:,qq)),real(dEplus(:,qq))-real(dEminus(:,qq)) );
end

%% Calculate probe amplitudes and measurement vector.
ampSq = (Iplus+Iminus)/2 - repmat(IupVec,[1,Npairs]);  % square of probe E-field amplitudes
ampSq(ampSq<0) = 0;  % If probe amplitude is zero, amplitude is zero there.
amp = sqrt(ampSq);   % E-field amplitudes, dimensions: [length(corr_inds), Npairs]
isnonzero = all(amp,2);
z = ( (Iplus-Iminus)/4).';  % Measurement vector, dimensions: [Npairs,length(corr_inds)]

for qq=1:Npairs % Display the actual probe intensity
    ampSq2D = zeros(mp.Neta,mp.Nxi); ampSq2D(corr_inds) = ampSq(:,qq); 
    fprintf('*** Avg meas probe %d intensity =\t%.3e \n',qq,mean(ampSq2D(corr_inds)));
end

%% Compute electric field in dark hole pixel by pixel.
zerosCounter = 0;
for ipix=1:length(corr_inds)
    H = zeros(Npairs,2); % Observation matrix
    if(isnonzero(ipix)==1) % Leave Eest for a pixel as zero if any probe amplitude is zero there.
        for qq=1:Npairs
            H(qq,:) = amp(ipix,qq)*[cos(dphdm(ipix,qq)), sin(dphdm(ipix,qq)) ];
        end
    else
        zerosCounter = zerosCounter+1;
        if( strcmpi(ep.probeAmpType,'both') )
            dE = dEplus(ipix,:).';
            H = [real(dE),imag(dE)]; %nonlinear way of calculating
        end
    end
    Epix = pinv(H)*z(:,ipix);
    Eest(ipix) = Epix(1) + 1i*Epix(2);
end
if(estvar.Itr>2)
    Eest(abs(Eest).^2>1e-4) = 0;  % If estimate is too bright, the estimate was bad.
end
% Eest = Eest.*isnonzero; 
fprintf('%d of %d pixels had zero probe amplitude. \n',zerosCounter,length(corr_inds)); 

% end % End of if-elseif statements for ep.probeAmpType
%% Compute the incoherent light
Eest_array(:,si) = Eest;
Iinco_array(:,si) =  IupVec-abs(Eest).^2;

end % End of loop over the wavelengths

% Calculate time spent on exposures
% estvar.timeTotal = Nsbp*(1+2*Npairs)*ip.ExpTime*ip.num_im; % Since unprobed imaged is half that time.

estvar.Eest_array = Eest_array;
estvar.Iinco_array = Iinco_array;
ep.InormMean(estvar.Itr)=ep.InormMean(estvar.Itr);
estvar.meanIprobed = meanIprobed;
% estvar.amp = amp/sqrt(InormProbe); % Normalized probe amplitude
estvar.IupMean= IunprobedSum/Nsbp; %--Bandpass-averaged image

estvar.ampSqMean = mean(mean(ampSq,1),2); % Mean probe intensity
% estvar.HmodelAll = HmodelAll;
% estvar.dDM2VprobePlus = dDM2VprobePlus;
% estvar.dElin = dElin;

% Calculate starlight estimated contrast
Istar2D = zeros(mp.Neta,mp.Nxi);
Istar2D(corr_inds) = abs(Eest).^2;
% Istar2D(Istar2D<0) = 0;
estvar.estContrast = sum(sum(mp.ScoreMask.*Istar2D))/sum(sum(mp.ScoreMask));


if(isfield(ep,'flagSave'))
    if(ep.flagSave) % Save images and DM commands out for later
%     cd(mp.folders.data)
%     %fnI = sprintf('img_run%03d_estItr%03d_im%04dto%04d.fits',ip.RunNum,estvar.Itr,ImNum0,ip.ImNum);
%     %fnV = sprintf('dmV_run%03d_estItr%03d_im%04dto%04d.fits',ip.RunNum,estvar.Itr,ImNum0,ip.ImNum);
%     fnI = sprintf('img_run%03d_estItr%03d.fits',mp.RunNum,estvar.Itr);
%     fnV = sprintf('dmV_run%03d_estItr%03d.fits',mp.RunNum,estvar.Itr);
%     fitswrite(IcubeOut,fnI);
%     fitswrite(VcubeOut,fnV);
%     cd(mp.folders.m)
    end
end

fprintf(' done. Time: %.3f\n',toc);




