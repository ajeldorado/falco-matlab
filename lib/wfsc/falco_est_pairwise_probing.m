% Copyright 2018, 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to estimate the final focal plane electric field via 
%  pair-wise probing and batch process estimation.
%
%--References for the algorithms and their usage:
% A. Give'on, B. Kern, and S. Shaklan, "Pair-wise, deformable mirror, 
% image plane-based diversity electric field estimation for high contrast 
% coronagraphy," in Proceedings of SPIE, vol. 8151, p. 815110, 2011.
%
% T. D. Groff and N. J. Kasdin, "Kalman filtering techniques for focal plane 
% electric field estimation," Journal of the Optical Society of America A, 
% vol. 30, no. 1, pp. 128-139, 2013.
%
%
%--INPUTS
% x_in: Scalar value.
%
%--OUTPUTS
%  x_out: even-valued integer value
%
%--REVISION HISTORY
% Modified on 2019-02-25 by A.J. Riggs to include the batch process
%   estimator and Kalman filter in the same file since their setup is the same.
% Modified on 2019-02-06 by A.J. Riggs for the updated FALCO syntax.
% Modified on 2018-04-23 by A.J. Riggs from the Princeton HCIL lab code.
% Created on 2015-02-19 by A.J. Riggs at Princeton University.
%
%
%--New variables for Kalman filter
%  - mp.est.ItrStartKF:  Which correction iteration to start recursive estimate
%  - mp.est.tExp
%  - mp.est.num_im
%  - mp.readNoiseStd
%  - mp.peakCountsPerPixPerSec
%  - mp.est.Qcoef
%  - mp.est.Rcoef

function [ev] = falco_est_pairwise_probing(mp,ev,varargin)

%--If there is a third input, it is the Jacobian structure
if( size(varargin, 2)==1 )
    jacStruct = varargin{1};
end

%--"ev" is passed in only for the Kalman filter. Clear it for the batch
% process to avoid accidentally using old data.
switch lower(mp.estimator)
    case{'pwp-bp'}
        clear ev
end

%--Select number of actuators across based on chosen DM for the probing
if(mp.est.probe.whichDM==1)
    Nact = mp.dm1.Nact;
elseif(mp.est.probe.whichDM==2)
    Nact = mp.dm2.Nact;
end

%--Store the initial DM commands
if(any(mp.dm_ind==1));  DM1Vnom = mp.dm1.V;  end
if(any(mp.dm_ind==2));  DM2Vnom = mp.dm2.V;  else; DM2Vnom = zeros(size(mp.dm1.V)); end

% Definitions:
Npairs = mp.est.probe.Npairs; % % Number of image PAIRS for DM Diversity or Kalman filter initialization
ev.Icube = zeros(mp.Fend.Neta,mp.Fend.Nxi,1+2*Npairs);
if(any(mp.dm_ind==1));  ev.Vcube.dm1 = zeros(mp.dm1.Nact,mp.dm1.Nact,1+2*Npairs);  end
if(any(mp.dm_ind==2));  ev.Vcube.dm2 = zeros(mp.dm2.Nact,mp.dm2.Nact,1+2*Npairs);  end

%--Generate evenly spaced probes along the complex unit circle
% NOTE: Nprobes=Npairs*2;   
probePhaseVec = [0 Npairs];
for k = 1:Npairs-1
    probePhaseVec = [probePhaseVec probePhaseVec(end)-(Npairs-1)]; % #ok<AGstroke2ROW>
    probePhaseVec = [probePhaseVec probePhaseVec(end)+(Npairs)]; % #ok<AGstroke2ROW>
end
probePhaseVec = probePhaseVec*pi/(Npairs);

switch lower(mp.est.probe.axis)
    case 'y'
        badAxisVec = repmat('y',[2*Npairs,1]);
    case 'x'
        badAxisVec = repmat('x',[2*Npairs,1]);
    case{'alt','xy','alternate'}
        badAxisVec = repmat('x',[2*Npairs,1]);
        badAxisVec(3:4:end) = 'y';
        badAxisVec(4:4:end) = 'y';
    case 'multi'
        badAxisVec = repmat('m',[2*Npairs,1]);
end

%% Initialize output arrays
ev.Eest = zeros(mp.Fend.corr.Npix,mp.Nsbp);
ev.IincoEst = zeros(mp.Fend.corr.Npix,mp.Nsbp);
ev.I0mean = 0;
ev.IprobedMean = 0;

%% Get images and perform estimates in each sub-bandpass

fprintf('Estimating electric field with batch process estimation ...\n'); tic;

for si=1:mp.Nsbp
    fprintf('Wavelength: %u/%u ... ',si,mp.Nsbp);

    % Valid for all calls to model_compact.m:
    modvar.sbpIndex = si;
    modvar.whichSource = 'star';

    %% Measure current contrast level average, and on each side of Image Plane
    % Reset DM commands to the unprobed state:
    mp.dm1.V = DM1Vnom;
    mp.dm2.V = DM2Vnom;
    %% Separate out values of images at dark hole pixels and delta DM voltage settings
    Iplus  = zeros( [mp.Fend.corr.Npix, Npairs]); % Pixels of plus probes' intensities
    Iminus = zeros( [mp.Fend.corr.Npix, Npairs]); % Pixels of minus probes' intensities
    DM1Vplus  = zeros( [Nact,Nact, Npairs]);
    DM1Vminus = zeros( [Nact,Nact, Npairs]);
    DM2Vplus  = zeros( [Nact,Nact, Npairs]);
    DM2Vminus = zeros( [Nact,Nact, Npairs]);

    %% Compute probe shapes and take probed images:

    %--Take initial, unprobed image (for unprobed DM settings).
    whichImg = 1;
    I0 = falco_get_sbp_image(mp,si);
    I0vec = I0(mp.Fend.corr.maskBool); % Vectorize the correction region pixels
    ev.I0mean = ev.I0mean+I0/mp.Nsbp; %--Getting the sub-bandpass-averaged Inorm

    %--Store values for first image and its DM commands
    ev.Icube(:,:, whichImg) = I0;
    if(any(mp.dm_ind==1));  ev.Vcube.dm1(:,:, whichImg) = mp.dm1.V;  end
    if(any(mp.dm_ind==2));  ev.Vcube.dm2(:,:, whichImg) = mp.dm2.V;  end

    %--Compute the average Inorm in the scoring and correction regions
    ev.score.Inorm = mean(I0(mp.Fend.score.maskBool));
    ev.corr.Inorm  = mean(I0(mp.Fend.corr.maskBool));
    fprintf('Measured unprobed Inorm (Corr / Score): %.2e \t%.2e \n',ev.corr.Inorm,ev.score.Inorm);    

    % Set (approximate) probe intensity based on current measured Inorm
    if(mp.flagFiber)
        ev.InormProbeMax = 1e-5;
        InormProbe = min([sqrt(max(I0)*1e-8), ev.InormProbeMax]);
    else
        ev.InormProbeMax = 1e-4;
        InormProbe = min( [sqrt(max(I0vec)*1e-5), ev.InormProbeMax]); %--Change this to a high percentile value (e.g., 90%) instead of the max to avoid being tricked by noise
    end
    fprintf('Chosen probe intensity: %.2e \n',InormProbe);    

    %--Perform the probing
    iOdd=1; iEven=1; %--Initialize index counters
    for iProbe=1:2*Npairs           

        %--Generate the command map for the probe
        probeCmd = falco_gen_pairwise_probe(mp,InormProbe,probePhaseVec(iProbe),badAxisVec(iProbe));

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

        %--Take probed image
        if(mp.flagFiber)
            Im = falco_get_sbp_image_fiber(mp,si);
        else
            Im = falco_get_sbp_image(mp,si);
        end
        ImNonneg = Im; ImNonneg(Im<0) = 0;
        whichImg = 1+iProbe; %--Increment image counter
        ev.IprobedMean = ev.IprobedMean + mean(Im(mp.Fend.corr.maskBool))/(2*Npairs); %--Inorm averaged over all the probed images

        %--Store probed image and its DM settings
        ev.Icube(:,:,whichImg) = Im;
        if(any(mp.dm_ind==1));  ev.Vcube.dm1(:,:,whichImg) = mp.dm1.V;  end
        if(any(mp.dm_ind==2));  ev.Vcube.dm2(:,:,whichImg) = mp.dm2.V;  end

        %--Report results
        probeSign = ['-','+'];
        fprintf('Actual Probe %d%s Contrast is: %.2e \n',ceil(iProbe/2),probeSign(mod(iProbe,2)+1),mean(Im(mp.Fend.corr.maskBool)));

        %--Assign image to positive or negative probe collection:
        if mod(iProbe,2)==1  % Odd; for plus probes
            if(any(mp.dm_ind==1));  DM1Vplus(:,:,iOdd) = dDM1Vprobe + DM1Vnom;  end
            if(any(mp.dm_ind==2));  DM2Vplus(:,:,iOdd) = dDM2Vprobe + DM2Vnom;  end
            Iplus(:,iOdd) = Im(mp.Fend.corr.maskBool);
            iOdd=iOdd+1;
        elseif mod(iProbe,2)==0  % Even; for minus probes
            if(any(mp.dm_ind==1));  DM1Vminus(:,:,iEven) = dDM1Vprobe + DM1Vnom;  end
            if(any(mp.dm_ind==2));  DM2Vminus(:,:,iEven) = dDM2Vprobe + DM2Vnom;  end 
            Iminus(:,iEven) = Im(mp.Fend.corr.maskBool);
            iEven=iEven+1;
        end
    end

    %% Calculate probe amplitudes and measurement vector. (Refer again to Give'on+ SPIE 2011 to undersand why.)
    ampSq = (Iplus+Iminus)/2 - repmat(I0vec,[1,Npairs]);  % square of probe E-field amplitudes
    ampSq(ampSq<0) = 0;  % If probe amplitude is zero, amplitude is zero there.
    amp = sqrt(ampSq);   % E-field amplitudes, dimensions: [mp.Fend.corr.Npix, Npairs]
    isnonzero = all(amp,2);
    zAll = ( (Iplus-Iminus)/4).';  % Measurement vector, dimensions: [Npairs,mp.Fend.corr.Npix]
    ampSq2Dcube = zeros(mp.Fend.Neta,mp.Fend.Nxi,mp.est.probe.Npairs);
    for iProbe=1:Npairs % Display the actual probe intensity
        ampSq2D = zeros(mp.Fend.Neta,mp.Fend.Nxi); ampSq2D(mp.Fend.corr.maskBool) = ampSq(:,iProbe); 
        ampSq2Dcube(:,:,iProbe) = ampSq2D;
        fprintf('*** Mean measured Inorm for probe #%d  =\t%.3e \n',iProbe,mean(ampSq2D(mp.Fend.corr.maskBool)));
    end

    %% Plot relevant data for all the probes
    falco_plot_pairwise_probes(mp,ev,DM1Vplus-repmat(DM1Vnom,[1,1,size(DM1Vplus,3)]),ampSq2Dcube)

    %% Perform the estimation
    
    if(mp.est.flagUseJac) %--Use Jacobian for estimation. This is fully model-based if the Jacobian is purely model-based, or it is better if the Jacobian is adaptive based on empirical data.
        
        dEplus  = zeros(size(Iplus ));
        for iProbe=1:Npairs
            if(mp.est.probe.whichDM == 1)
                dV = DM1Vplus(:,:,iProbe)-DM1Vnom;
                dEplus(:,iProbe) = squeeze(jacStruct.G1(:,:,si))*dV(mp.dm1.act_ele);
            elseif(mp.est.probe.whichDM == 2)
                dV = DM2Vplus(:,:,iProbe)-DM2Vnom;
                dEplus(:,iProbe) = squeeze(jacStruct.G2(:,:,si))*dV(mp.dm2.act_ele);
            end
        end
        
    else %--Get the probe phase from the model and the probe amplitude from the measurements

        % For unprobed field based on model:
        if(any(mp.dm_ind==1));  mp.dm1.V = DM1Vnom;  end
        if(any(mp.dm_ind==2));  mp.dm2.V = DM2Vnom;  end % Added July 9, 2014
        if(mp.flagFiber)
            [~, E0] = model_compact(mp, modvar);
        else
            E0 = model_compact(mp, modvar);
        end
        E0vec = E0(mp.Fend.corr.maskBool);

        %--For probed fields based on model:
        Eplus  = zeros(size(Iplus ));
        Eminus = zeros(size(Iminus));
        for iProbe=1:Npairs
            % For plus probes:
            if(any(mp.dm_ind==1));  mp.dm1.V = squeeze( DM1Vplus(:,:,iProbe));  end
            if(any(mp.dm_ind==2));  mp.dm2.V = squeeze( DM2Vplus(:,:,iProbe));  end
            if(mp.flagFiber)
                [~, Etemp] = model_compact(mp, modvar);
            else
                Etemp = model_compact(mp, modvar);
            end
            Eplus(:,iProbe) = Etemp(mp.Fend.corr.maskBool);
            % For minus probes:
            if(any(mp.dm_ind==1));  mp.dm1.V = squeeze( DM1Vminus(:,:,iProbe));  end
            if(any(mp.dm_ind==2));  mp.dm2.V = squeeze( DM2Vminus(:,:,iProbe));  end
            if(mp.flagFiber)
                [~, Etemp] = model_compact(mp, modvar);
            else
                Etemp = model_compact(mp, modvar);
            end
            Eminus(:,iProbe) = Etemp(mp.Fend.corr.maskBool);
        end

        %%--Create delta E-fields for each probe image. Then create Npairs phase angles.
        dEplus  = Eplus  - repmat(E0vec,[1,Npairs]);
        dEminus = Eminus - repmat(E0vec,[1,Npairs]);
        dphdm  = zeros( [mp.Fend.corr.Npix, Npairs]); %--phases of the probes
        for iProbe=1:Npairs
            dphdm(:,iProbe) = atan2( imag(dEplus(:,iProbe))-imag(dEminus(:,iProbe)),real(dEplus(:,iProbe))-real(dEminus(:,iProbe)) );
        end
        
    end 
    
%% Batch process the measurements to estimate the electric field in the dark hole. Done pixel by pixel.

if( strcmpi(mp.estimator,'pwp-bp') || (strcmpi(mp.estimator,'pwp-kf') && ev.Itr<mp.est.ItrStartKF ) )
    Eest = zeros(mp.Fend.corr.Npix,1);
    zerosCounter = 0;
    for ipix=1:mp.Fend.corr.Npix
        
        if(mp.est.flagUseJac) 
            dE = dEplus(ipix,:).';
            H = [real(dE),imag(dE)];
        else
            H = zeros(Npairs,2); % Observation matrix
            if(isnonzero(ipix)==1) % Leave Eest for a pixel as zero if any probe amplitude is zero there.
                for iProbe=1:Npairs
                    H(iProbe,:) = amp(ipix,iProbe)*[cos(dphdm(ipix,iProbe)), sin(dphdm(ipix,iProbe)) ];
                end
            else
                zerosCounter = zerosCounter+1;
            end
        end
        Epix = pinv(H)*zAll(:,ipix); %--Batch process estimation
        Eest(ipix) = Epix(1) + 1i*Epix(2);
    end
    Eest(abs(Eest).^2 > 1e-2) = 0;  % If estimate is too bright, the estimate was probably bad. !!!!!!!!!!!!!!BE VERY CAREFUL WITH THIS HARD-CODED VALUE!!!!!!!!!!!!!!!
    fprintf('%d of %d pixels were given zero probe amplitude. \n',zerosCounter,mp.Fend.corr.Npix); 

    %--Initialize the state and state covariance estimates for Kalman
    %filter. The state is the real and imag parts of the E-field.
    if(strcmpi(mp.estimator,'pwp-kf'))
        %--Re-organize the batch-processed E-field estimate into the first state estimate for the Kalman filter
        xOld = zeros(2*mp.Fend.corr.Npix,1);
        for ii = 1:mp.Fend.corr.Npix
            xOld(2*(ii-1)+1:2*(ii-1)+2) = [real(Eest(ii)); imag(Eest(ii))];
        end
        ev.xOld = xOld; %--Save out for returning later

        %--Initialize the state covariance matrix (2x2 for each dark hole pixel)
        ev.Pold_KF_array = repmat(mp.est.Pcoef0*eye(2), [mp.Fend.corr.Npix,1,mp.Nsbp]);
    end

end   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--Kalman Filter Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
if(strcmpi(mp.estimator,'pwp-kf') && (ev.Itr>=mp.est.ItrStartKF) )
    
    xOld = ev.xOld;
    Pold = ev.Pold_KF_array(:,:,si);
    
    %--Construct the observation matrix, H, for all pixels
    Hall = zeros(Npairs,2,mp.Fend.corr.Npix);
    for ipix=1:mp.Fend.corr.Npix
        if(mp.est.flagUseJac)
            dE = dEplus(ipix,:).'; % dE based on the Jacobian only
            H = [real(dE),imag(dE)]; % Observation matrix
        else
            H = zeros(Npairs,2); % Observation matrix
            if(isnonzero(ipix)==1) % Leave Eest for a pixel as zero if any probe amplitude is zero there.
                for qq=1:Npairs
                    H(qq,:) = amp(ipix,qq)*[cos(dphdm(ipix,qq)), sin(dphdm(ipix,qq)) ];
                end
            end
        end
        Hall(:,:,ipix) = H;    
    end
    
    %--Compute the change in E-field since last correction iteration
    if(mp.est.flagUseJac) %--Use Jacobian to compute delta E-field from previous correction step to now.
        if(mp.est.probe.whichDM == 1)
            dE = squeeze(jacStruct.G1(:,:,si))*mp.dm1.dV(mp.dm1.act_ele);
        elseif(mp.est.probe.whichDM == 2)
            dE = squeeze(jacStruct.G2(:,:,si))*mp.dm2.dV(mp.dm2.act_ele);
        end
    else
        % For Xminus, use nonlinear dynamics instead of Gamma. This means
        % difference the output of model_compact rather than using the
        % Jacobian.
        %--Previous unprobed field based on model:
        if(any(mp.dm_ind==1));  mp.dm1.V = DM1Vnom-mp.dm1.dV;  end
        if(any(mp.dm_ind==2));  mp.dm2.V = DM2Vnom-mp.dm2.dV;  end
        if(mp.flagFiber)
            [~, Eprev] = model_compact(mp, modvar);
        else
            Eprev = model_compact(mp, modvar);
        end
        EprevVec = Eprev(mp.Fend.corr.maskBool);
        dE = E0vec-EprevVec; % Change in unprobed E-field between correction iterations
    end

    %--Construct dX, the change in state, from dE
    dX = zeros(size(xOld));
    for ii=1:mp.Fend.corr.Npix
       dX(2*(ii-1)+1:2*(ii-1)+2) = [real(dE(ii)); imag(dE(ii))];
    end

    %--Compute Sensor Noise, R. You can calculate this from the properties
    % of your optical system.
	% ncounts_readout = mp.readNoiseStd;  % Read noise in counts (ADU). Might need to be in photo-electrons
    % ncounts_shot = sqrt(ev.IprobedMean*mp.peakCountsPerPixPerSec);
    % Dark current not included here (yet).
    ncounts_std = sqrt( (sqrt(2)*ev.IprobedMean*mp.peakCountsPerPixPerSec*mp.est.tExp + mp.readNoiseStd^2)/mp.est.num_im);
    Rvar = (ncounts_std/(mp.peakCountsPerPixPerSec*mp.est.tExp))^2; % Don't forget to square it since R = E<n*n.'>. This is a variable scalar
    Rmat = mp.est.Rcoef*Rvar*eye(Npairs); % A.J.'s way, used in v2
    fprintf('Sensor noise coefficient: %.3e\n',Rmat(1,1));

    %--Compute Process Noise, Q. This you just have to tune (via the scalar mp.est.Qcoef).
    Q00 = ev.corr.Inorm; %--Set coefficient to the current measured contrast. Can change to the last estimated contrast.
    Q = Q00*mp.est.Qcoef*repmat( eye(2), [mp.Fend.corr.Npix, 1]);
    dP = Q;
    fprintf('Process noise coefficient: %.3e\n',Q(1,1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute Kalman Filter Equations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xNew = zeros(2*mp.Fend.corr.Npix,1);
    Pnew = zeros(2*mp.Fend.corr.Npix,2);
    Kall = zeros(2,Npairs,mp.Fend.corr.Npix); %--Store the Kalman gain matrices. For diagnostics only

    xMinusAll = xOld + dX; %State Estimate Extrapolation (Recall Phi = Identity)
    PminusAll = Pold + dP; %Covariance Estimate Extrapolation (Recall Phi = Identity)
    for ii=1:mp.Fend.corr.Npix
        %--Format several matrices
        z = zAll(:,ii); % z has dimensions [Npairs x Npix]
        H = Hall(:,:,ii); % H has dimensions [Npairs x 2]
        xMinus = xMinusAll(2*(ii-1)+1:2*(ii-1)+2); % x has dimensions [2 x 1]
        Pminus = PminusAll(2*(ii-1)+1:2*(ii-1)+2, 1:2); % P has dimensions [2 x 2]
        
        %--Compute the Kalman gain matrix.
        K = Pminus*H.'/(H*Pminus*H.' + Rmat); % dimensions [2 x Npairs]
        Kall(:,:,ii) = K; % Store Kalman gains for diagnostics only
        
        %--State Estimate Update
        xNew(2*(ii-1)+1:2*(ii-1)+2) = xMinus + K*(z-H*xMinus);
        Pnew(2*(ii-1)+1:2*(ii-1)+2,1:2) = (eye(2)-K*H)*Pminus;
    end

    %--For diagnostics only. Print the mean state covariance matrix to the command line
    P11 = mean(Pnew(1:2:end,1));
    P22 = mean(Pnew(2:2:end,2));
    P12 = mean(Pnew(1:2:end,2));
    P21 = mean(Pnew(2:2:end,1));
    fprintf('P11,P22,P12,P21: \t%.2e \t%.2e \t%.2e \t%.2e \n',P11,P22,P12,P21);

    % Re-order Xnew into Eest
    EestKF = zeros(mp.Fend.corr.Npix,1); %for re-ordering state for control
    for jj=1:mp.Fend.corr.Npix
        EestKF(jj) = xNew(2*(jj-1)+1) + 1i*xNew(2*(jj-1)+2);
    end
    Eest = EestKF;
    
    %--Save out new state covariance P as the old for the next correction iteration.
    ev.Pold_KF_array(:,:,si) = Pnew;

    fprintf(['mean(abs(K)) = ' num2str(mean(mean(mean(abs(Kall))))) '.\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %End Kalman Filter Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save out the estimates
ev.Eest(:,si) = Eest;
ev.IincoEst(:,si) =  I0vec-abs(Eest).^2; %--Compute the incoherent light

end %--End of loop over the wavelengths

%--Other data to save out
ev.ampSqMean = mean(ampSq(:)); %--Mean probe intensity
ev.ampNorm = amp/sqrt(InormProbe); %--Normalized probe amplitude maps

%--Calculate the mean normalized intensity over the whole dark hole at all
% wavelengths.
ev.Iest = abs(ev.Eest).^2;
ev.InormEst = mean(ev.Iest(:));

fprintf(' done. Time: %.3f\n',toc);

end %--END OF FUNCTION