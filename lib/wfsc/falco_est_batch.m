% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to estimate the final focal plane electric field via 
%  pair-wise probing and batch process estimation.
%
%--Reference for the algorithm and its usage:
% A. Give?on, B. Kern, and S. Shaklan, "Pair-wise, deformable mirror, 
% image plane-based diversity electric field estimation for high contrast 
% coronagraphy," in Proceedings of SPIE, vol. 8151, p. 815110, 2011.
% 
%
%--INPUTS
% x_in: Scalar value.
%
%--OUTPUTS
%  x_out: even-valued integer value
%
%--REVISION HISTORY
% Modified on 2019-02-06 by A.J. Riggs for the updated FALCO syntax.
% Modified on 2018-04-23 by A.J. Riggs from the Princeton HCIL lab code.
% Created on 2015-02-19 by A.J. Riggs at Princeton University.
%
%
%--New variables
%    -mp.est.probe.Npairs:      Number of pairwise probe pairs to use
%    -mp.est.probe.whichDM:      Which DM # to use for probing. 1 or 2. Default is 1
%    -mp.est.probe.axis:    which axis to have the phase discontinuity along [x or y]

function [ip] = falco_est_batch(mp)

%--Select number of actuators across based on chosen DM for the probing
if(mp.est.probe.whichDM==1)
    Nact = mp.dm1.Nact;
elseif(mp.est.probe.whichDM==2)
    Nact = mp.dm2.Nact;
end

%--Store the initial DM commands
if(any(mp.dm_ind==1));  DM1Vnom = mp.dm1.V;  end
if(any(mp.dm_ind==2));  DM2Vnom = mp.dm2.V;  end

% Definitions:
Npairs = mp.est.probe.Npairs; % % Number of image PAIRS for DM Diversity or Kalman filter initialization
ip.Icube = zeros(mp.F4.Neta,mp.F4.Nxi,1+2*Npairs);
if(any(mp.dm_ind==1));  ip.Vcube.dm1 = zeros(mp.dm1.Nact,mp.dm1.Nact,1+2*Npairs);  end
if(any(mp.dm_ind==2));  ip.Vcube.dm2 = zeros(mp.dm2.Nact,mp.dm2.Nact,1+2*Npairs);  end

%--Generate evenly spaced probes along the complex unit circle
% NOTE: Nprobes=Npairs*2;   
probePhaseVec = [0 Npairs];
for k = 1:Npairs-1
    probePhaseVec = [probePhaseVec probePhaseVec(end)-(Npairs-1)];% % #ok<AGstroke2ROW>
    probePhaseVec = [probePhaseVec probePhaseVec(end)+(Npairs)]; % #ok<AGstroke2ROW>
end
probePhaseVec = probePhaseVec*pi/(Npairs);
% ProbePhsVec = (offset + ProbePhsVec)*pi/(Npairs);

switch lower(mp.est.probe.axis)
    case 'y'
        badAxisVec = repmat('y',[2*Npairs,1]);
    case 'x'
        badAxisVec = repmat('x',[2*Npairs,1]);
    case{'alt','xy','alternate'}
        badAxisVec = repmat('x',[2*Npairs,1]);
        badAxisVec(3:4:end) = 'y';
        badAxisVec(4:4:end) = 'y';
end

%% Initialize output arrays
ip.Eest = zeros(mp.F4.corr.Npix,mp.Nsbp);
ip.IincoEst = zeros(mp.F4.corr.Npix,mp.Nsbp);
ip.I0mean = 0;
ip.IprobedMean = 0;

%% Get images and perform estimates in each sub-bandpass

fprintf('Estimating electric field with batch process estimation ...\n'); tic;

for si=1:mp.Nsbp
    fprintf('Wavelength: %u/%u ... ',si,mp.Nsbp);

    % % Valid for all calls to model_compact.m:
    modvar.sbpIndex = si;
    modvar.whichSource = 'star';

    %% Measure current contrast level average, and on each side of Image Plane
    % Reset DM commands to the unprobed state:
    mp.dm1.V = DM1Vnom;
    mp.dm2.V = DM2Vnom;

    %% Separate out values of images at dark hole pixels and delta DM voltage settings
    Iplus  = zeros( [mp.F4.corr.Npix, Npairs]); % Pixels of plus probes' intensities
    Iminus = zeros( [mp.F4.corr.Npix, Npairs]); % Pixels of minus probes' intensities
    DM1Vplus  = zeros( [Nact,Nact, Npairs]);
    DM1Vminus = zeros( [Nact,Nact, Npairs]);
    DM2Vplus  = zeros( [Nact,Nact, Npairs]);
    DM2Vminus = zeros( [Nact,Nact, Npairs]);

    %% Compute probe shapes and take probed images:


    if( (mp.flagSim==true) ) %--Do if in simulation

        %--Take initial, unprobed image (for unprobed DM settings).
        whichImg = 1;
        I0 = falco_get_sbp_image(mp,si);
        I0vec = I0(mp.F4.corr.maskBool); % Vectorize the correction region pixels
        ip.I0mean = ip.I0mean+I0/mp.Nsbp; %--Getting the sub-bandpass-averaged Inorm

        %--Store values for first image and its DM commands
        ip.Icube(:,:, whichImg) = I0;
        if(any(mp.dm_ind==1));  ip.Vcube.dm1(:,:, whichImg) = mp.dm1.V;  end
        if(any(mp.dm_ind==2));  ip.Vcube.dm2(:,:, whichImg) = mp.dm2.V;  end

        %--Compute the average Inorm
        ip.score.Inorm = mean(I0(mp.F4.score.maskBool));
        ip.corr.Inorm  = mean(I0(mp.F4.corr.maskBool));
        fprintf('Measured unprobed Inorm (Corr / Score): %.2e \t%.2e \n',ip.corr.Inorm,ip.score.Inorm);    

        % Set (approximate) probe intensity based on current measured Inorm
        ip.InormProbeMax = 1e-4;
        InormProbe = min( [sqrt(ip.score.Inorm*1e-5), ip.InormProbeMax]);
        fprintf('Chosen probe intensity: %.2e \n',InormProbe);    

        
        %--Perform the probing
        iOdd=1; iEven=1; %--Initialize index counters
        for iProbe=1:2*Npairs           
            
            %--Generate the command map for the probe
            probeCmd = falco_gen_pairwise_probe(mp,InormProbe,probePhaseVec(iProbe),badAxisVec(iProbe));

    %         dx = mp.Ddm/2/mp.Ndm;
    %         xs = (-mp.Ndm:mp.Ndm-1)'*dx + dx/2;
    %         XS = repmat(xs.',2*mp.Ndm,1);
    %         YS = XS.'; 
    %         ProbeH=hcil_makeProbeSurf(ep.ProbeArea,mp.Ddm,...
    %             mp.sbp_centers(si),probePhaseVec(qq),ep.ProbeoffsetX,ep.ProbeoffsetY,XS,YS,InormProbe); %Function used to compute the Probe Shape. Output is in meters.
    %         ProbeH = 1*ProbeH; % Scale the probe amplitude empirically if needed
    %         probeCmd = imresize(ProbeH,Nact/length(XS));  % crude way for now, in units of meters

            %--Select which DM to use for probing. Allocate probe to that DM
            if(mp.est.probe.whichDM == 1)
                dDM1Vprobe = probeCmd./mp.dm1.VtoH; % Now in volts
                dDM2Vprobe = 0;
            elseif(mp.est.probe.whichDM == 2)
                dDM1Vprobe = 0;        
                dDM2Vprobe = probeCmd./mp.dm1.VtoH; % Now in volts
            end
            if(any(mp.dm_ind==1));  mp.dm1.V = DM1Vnom+dDM1Vprobe;  end
            if(any(mp.dm_ind==2));  mp.dm2.V = DM2Vnom+dDM2Vprobe;  end
            %figure(202); imagesc(dDM1Vprobe); axis xy equal tight; colorbar; set(gca,'Fontsize',20); drawnow;

            %--Take probed image
            Im = falco_get_sbp_image(mp,si);
            whichImg = 1+iProbe; %--Increment image counter
            ip.IprobedMean = ip.IprobedMean + mean(Im(mp.F4.corr.maskBool))/(2*Npairs); %--Inorm averaged over wavelength for the probed images
            if(mp.flagPlot);  figure(203); imagesc(log10(Im),[-8 -3]); axis xy equal tight; colorbar; set(gca,'Fontsize',20); drawnow;  end

            %--Store probed image and its DM settings
            ip.Icube(:,:,whichImg) = Im;
            if(any(mp.dm_ind==1));  ip.Vcube.dm1(:,:,whichImg) = mp.dm1.V;  end
            if(any(mp.dm_ind==2));  ip.Vcube.dm2(:,:,whichImg) = mp.dm2.V;  end
            
            %--Report results
            probeSign = ['-','+'];
            fprintf('Actual Probe %d%s Contrast is: %.2e \n',ceil(iProbe/2),probeSign(mod(iProbe,2)+1),mean(Im(mp.F4.corr.maskBool)));

            %--Assign image to positive or negative probe collection:
            if mod(iProbe,2)==1  % Odd; for plus probes
                if(any(mp.dm_ind==1));  DM1Vplus(:,:,iOdd) = dDM1Vprobe + DM1Vnom;  end
                if(any(mp.dm_ind==2));  DM2Vplus(:,:,iOdd) = dDM2Vprobe + DM2Vnom;  end
                Iplus(:,iOdd) = Im(mp.F4.corr.maskBool);
                iOdd=iOdd+1;
            elseif mod(iProbe,2)==0  % Even; for minus probes
                if(any(mp.dm_ind==1));  DM1Vminus(:,:,iEven) = dDM1Vprobe + DM1Vnom;  end
                if(any(mp.dm_ind==2));  DM2Vminus(:,:,iEven) = dDM2Vprobe + DM2Vnom;  end 
                Iminus(:,iEven) = Im(mp.F4.corr.maskBool);
                iEven=iEven+1;      
            end
        end

    else %--Testbed mode

    end %--End of "if statement" for simulation or testbed


    %% Perform the estimation
    % figure(24); imagesc(log10(abs(Im))); axis xy equal tight; colorbar;
    % pause(1);

    % %% Use model-based or measurement-based probes, or both
    % Eest = zeros(mp.F4.corr.Npix,1);
    % if( strcmpi(ep.probeAmpType,'modelJac') )
    %         %%G2 = zeros(mp.F4.corr.Npix,mp.Nact^2,mp.Nsbp)
    %         dDM2VprobePlus = DM2Vplus-repmat(DM2Vnom,[1,1,Npairs]);
    %         JacDM2 = ep.G2(:,:,si);
    %         HmodelAll = zeros(Npairs,2,mp.F4.corr.Npix); % Observation matrix
    %         for qq=1:Npairs
    %             dDM2VprobeVec = reshape(dDM2VprobePlus(:,:,qq),[mp.Nact^2,1]);
    %             dElin = JacDM2*dDM2VprobeVec; % dimensions: [Npix x 1]
    %             for jj=1:mp.F4.corr.Npix
    %                 HmodelAll(qq,:,jj) =  [real(dElin(jj)); imag(dElin(jj))];
    %             end
    %         end
    %         if( strcmpi(ep.probeAmpType,'modelJac') )
    %             for jj=1:mp.F4.corr.Npix
    %                 Epix = pinv(HmodelAll(:,:,jj))*z(:,jj);
    %                 Eest(jj) = Epix(1) + 1i*Epix(2); 
    %             end
    %         end
    %     
    % elseif( strcmpi(ep.probeAmpType,'meas') || strcmpi(ep.probeAmpType,'both') )

    %% Propagate each DM setting to the image plane. (model based)
    % For unprobed field based on model:
    if(any(mp.dm_ind==1));  mp.dm1.V = DM1Vnom;  end
    if(any(mp.dm_ind==2));  mp.dm2.V = DM2Vnom;  end % Added July 9, 2014
    E0 = model_compact(mp, modvar);
    E0vec = E0(mp.F4.corr.maskBool);

    %--For probed fields based on model:
    Eplus  = zeros(size(Iplus ));
    Eminus = zeros(size(Iminus));
    for iProbe=1:Npairs
        % For plus probes:
        if(any(mp.dm_ind==1));  mp.dm1.V = squeeze( DM1Vplus(:,:,iProbe));  end
        if(any(mp.dm_ind==2));  mp.dm2.V = squeeze( DM2Vplus(:,:,iProbe));  end
        Etemp = model_compact(mp, modvar);
        Eplus(:,iProbe) = Etemp(mp.F4.corr.maskBool);
        % For minus probes:
        if(any(mp.dm_ind==1));  mp.dm1.V = squeeze( DM1Vminus(:,:,iProbe));  end
        if(any(mp.dm_ind==2));  mp.dm2.V = squeeze( DM2Vminus(:,:,iProbe));  end
        Etemp = model_compact(mp, modvar);
        Eminus(:,iProbe) = Etemp(mp.F4.corr.maskBool);
    end

    %% Create delta E-fields for each probe image. Then create Npairs phase angles.
    dEplus  = Eplus  - repmat(E0vec,[1,Npairs]);
    dEminus = Eminus - repmat(E0vec,[1,Npairs]);
    dphdm  = zeros( [mp.F4.corr.Npix, Npairs]); %--phases of the probes
    for iProbe=1:Npairs
        dphdm(:,iProbe) = atan2( imag(dEplus(:,iProbe))-imag(dEminus(:,iProbe)),real(dEplus(:,iProbe))-real(dEminus(:,iProbe)) );
    end

    %% Calculate probe amplitudes and measurement vector. (Refer again to Give'on+ SPIE 2011 to undersand why.)
    ampSq = (Iplus+Iminus)/2 - repmat(I0vec,[1,Npairs]);  % square of probe E-field amplitudes
    ampSq(ampSq<0) = 0;  % If probe amplitude is zero, amplitude is zero there.
    amp = sqrt(ampSq);   % E-field amplitudes, dimensions: [mp.F4.corr.Npix, Npairs]
    isnonzero = all(amp,2);
    z = ( (Iplus-Iminus)/4).';  % Measurement vector, dimensions: [Npairs,mp.F4.corr.Npix]

    for iProbe=1:Npairs % Display the actual probe intensity
        ampSq2D = zeros(mp.F4.Neta,mp.F4.Nxi); ampSq2D(mp.F4.corr.maskBool) = ampSq(:,iProbe); 
        fprintf('*** Mean measured Inorm for probe #%d  =\t%.3e \n',iProbe,mean(ampSq2D(mp.F4.corr.maskBool)));
        if(mp.flagPlot);  figure(201); imagesc(ampSq2D); axis xy equal tight; colorbar; set(gca,'Fontsize',20); drawnow; pause(1);  end
    end
    

    %% Compute electric field in the dark hole pixel-by-pixel.
    Eest = zeros(mp.F4.corr.Npix,1);
    zerosCounter = 0;
    for ipix=1:mp.F4.corr.Npix
        H = zeros(Npairs,2); % Observation matrix
        if(isnonzero(ipix)==1) % Leave Eest for a pixel as zero if any probe amplitude is zero there.
            for iProbe=1:Npairs
                H(iProbe,:) = amp(ipix,iProbe)*[cos(dphdm(ipix,iProbe)), sin(dphdm(ipix,iProbe)) ];
            end
        else
            zerosCounter = zerosCounter+1;
%             if( strcmpi(ep.probeAmpType,'both') )
%                 dE = dEplus(ipix,:).';
%                 H = [real(dE),imag(dE)]; %nonlinear way of calculating
%             end
        end
        Epix = pinv(H)*z(:,ipix);
        Eest(ipix) = Epix(1) + 1i*Epix(2);
    end
    %if(estvar.Itr>2)
    Eest(abs(Eest).^2 > 1e-2) = 0;  % If estimate is too bright, the estimate was probably bad. !!!!!!!!!!!!!!BE VERY CAREFUL WITH THIS HARD-CODED VALUE!!!!!!!!!!!!!!!
    %end
    % Eest = Eest.*isnonzero; 
    fprintf('%d of %d pixels were given zero probe amplitude. \n',zerosCounter,mp.F4.corr.Npix); 

    % end % End of if-elseif statements for ep.probeAmpType
    %% Compute the incoherent light
    ip.Eest(:,si) = Eest;
    ip.IincoEst(:,si) =  I0vec-abs(Eest).^2;

end % End of loop over the wavelengths

% Calculate time spent on exposures
% estvar.timeTotal = mp.Nsbp*(1+2*Npairs)*ip.ExpTime*ip.num_im; % Since unprobed imaged is half that time.

ip.ampSqMean = mean(ampSq(:)); %mean(mean(ampSq,1),2); % Mean probe intensity

% Calculate starlight estimated contrast
ip.Istar2D = zeros(mp.F4.Neta,mp.F4.Nxi);
ip.Istar2D(mp.F4.corr.maskBool) = abs(ip.Eest).^2;
% Istar2D(Istar2D<0) = 0;

%--Calculate the mean normalized intensity over the whole dark hole at all
% wavelengths.
ip.Iest = abs(ip.Eest).^2;
ip.InormEst = mean(ip.Iest(:));

%--Save out all probing images and their DM voltage command maps
% if(isfield(ep,'flagSave'))
%     if(ep.flagSave) % Save images and DM commands out for later
% %     cd(mp.path.data)
% %     %fnI = sprintf('img_run%03d_estItr%03d_im%04dto%04d.fits',ip.RunNum,estvar.Itr,ImNum0,ip.ImNum);
% %     %fnV = sprintf('dmV_run%03d_estItr%03d_im%04dto%04d.fits',ip.RunNum,estvar.Itr,ImNum0,ip.ImNum);
% %     fnI = sprintf('img_run%03d_estItr%03d.fits',mp.RunNum,estvar.Itr);
% %     fnV = sprintf('dmV_run%03d_estItr%03d.fits',mp.RunNum,estvar.Itr);
% %     fitswrite(ip.Icube,fnI);
% %     fitswrite(VcubeOut,fnV);
% %     cd(mp.path.falco)
%     end
% end

fprintf(' done. Time: %.3f\n',toc);



end %--END OF FUNCTION


