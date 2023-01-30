% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Estimate the final focal plane electric field via pair-wise probing.
% The estimator itself can be a batch process or a Kalman filter.
%
% INPUTS
% ------
% mp : structure of model parameters
% ev : structure of estimation variables
% jacStruct : (optional) structure of control Jacobians
%
% RETURNS
% -------
% ev : structure of estimation variables
%
% REFERENCES
% ----------
% A. Give'on, B. Kern, and S. Shaklan, "Pair-wise, deformable mirror, 
% image plane-based diversity electric field estimation for high contrast 
% coronagraphy," in Proceedings of SPIE, vol. 8151, p. 815110, 2011.
%
% T. D. Groff and N. J. Kasdin, "Kalman filtering techniques for focal plane 
% electric field estimation," Journal of the Optical Society of America A, 
% vol. 30, no. 1, pp. 128-139, 2013.
%
% New variables for Kalman filter
%  - mp.est.ItrStartKF:  Which correction iteration to start recursive estimate
%  - mp.est.tExp
%  - mp.est.num_im
%  - mp.readNoiseStd
%  - mp.peakCountsPerPixPerSec
%  - mp.est.Qcoef
%  - mp.est.Rcoef

function ev = falco_est_pairwise_probing(mp, ev, varargin)

Itr = ev.Itr;
whichDM = mp.est.probe.whichDM;

% Number of elements to correct depending of regular pixels or fibers
if mp.flagFiber
    Ncorr = mp.Fend.Nfiber;
else
    Ncorr = mp.Fend.corr.Npix;
end

if strcmpi(mp.estimator, 'pwp-kf') || strcmpi(mp.estimator, 'pairwise-kf') || strcmpi(mp.estimator, 'pairwise-rect-kf')
    useKalmanFilter = true;
else
    useKalmanFilter = false;
end

%% Input checks
if ~isa(mp.est.probe, 'Probe')
    error('mp.est.probe must be an instance of class Probe')
end

% If scheduled, change some aspects of the probe.
% Empty values mean they are not scheduled.
if ~isempty(mp.est.probeSchedule.xOffsetVec)
    if length(mp.est.probeSchedule.xOffsetVec) < mp.Nitr
        error('mp.est.probeSchedule.xOffsetVec must have enough values for all WFSC iterations.')
    end
    mp.est.probe.xOffset = mp.est.probeSchedule.xOffsetVec(Itr);
end

if ~isempty(mp.est.probeSchedule.yOffsetVec)
    if length(mp.est.probeSchedule.yOffsetVec) < mp.Nitr
        error('mp.est.probeSchedule.yOffsetVec must have enough values for all WFSC iterations.')
    end
    mp.est.probe.yOffset = mp.est.probeSchedule.yOffsetVec(Itr);
end
fprintf('Probe offsets at the DM are (x=%.2f, y=%.2f) actuators.\n', mp.est.probe.xOffset, mp.est.probe.yOffset);

if ~isempty(mp.est.probeSchedule.rotationVec)
    if length(mp.est.probeSchedule.rotationVec) < mp.Nitr
        error('mp.est.probeSchedule.rotationVec must have enough values for all WFSC iterations.')
    end
    mp.est.probe.rotation = mp.est.probeSchedule.rotationVec(Itr);
end

if ~isempty(mp.est.probeSchedule.InormProbeVec)
    if length(mp.est.probeSchedule.InormProbeVec) < mp.Nitr
        error('mp.est.probeSchedule.InormProbeVec must have enough values for all WFSC iterations.')
    end
end

%--If there is a third input, it is the Jacobian structure
if size(varargin, 2) == 1
    jacStruct = varargin{1};
end

%%

%--"ev" is passed in only for the Kalman filter. Clear it for the batch
% process to avoid accidentally using old data.
switch lower(mp.estimator)
    case{'pairwise', 'pairwise-square', 'pwp-bp-square', 'pairwise-rect', 'pwp-bp', }
        clear ev
end

% Augment which DMs are used if the probing DM isn't used for control.
if whichDM == 1 && ~any(mp.dm_ind == 1)
    mp.dm_ind = [mp.dm_ind(:); 1];
elseif whichDM == 2 && ~any(mp.dm_ind == 2)
    mp.dm_ind = [mp.dm_ind(:); 2];
end

%--Select number of actuators across based on chosen DM for the probing
if whichDM == 1
    Nact = mp.dm1.Nact;
elseif whichDM == 2
    Nact = mp.dm2.Nact;
end

%--Store the initial DM commands
if any(mp.dm_ind == 1);  DM1Vnom = mp.dm1.V;  else; DM1Vnom = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2);  DM2Vnom = mp.dm2.V;  else; DM2Vnom = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1

% Initialize output arrays
Npairs = mp.est.probe.Npairs; % % Number of image PAIRS for DM Diversity or Kalman filter initialization
ev.imageArray = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1+2*Npairs, mp.Nsbp);
if mp.flagFiber
    ev.Ifiber = zeros(mp.Fend.Nfiber, 1);
end
ev.Eest = zeros(Ncorr, mp.Nsbp*mp.compact.star.count);
ev.IincoEst = zeros(Ncorr, mp.Nsbp*mp.compact.star.count);

ev.IprobedMean = 0;
ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
if whichDM == 1;  ev.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, 1+2*Npairs, mp.Nsbp);  end
if whichDM == 2;  ev.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, 1+2*Npairs, mp.Nsbp);  end

%--Generate evenly spaced probes along the complex unit circle
% NOTE: Nprobes=Npairs*2;   
probePhaseVec = [0 Npairs];
for k = 1:Npairs-1
    probePhaseVec = [probePhaseVec probePhaseVec(end)-(Npairs-1)];
    probePhaseVec = [probePhaseVec probePhaseVec(end)+(Npairs)];
end
probePhaseVec = probePhaseVec * pi / Npairs;

switch mp.estimator
    case{'pairwise', 'pairwise-square', 'pwp-bp-square', 'pwp-kf', 'pairwise-kf'}
        
        switch lower(mp.est.probe.axis)
            case 'y'
                badAxisVec = repmat('y', [2*Npairs, 1]);
            case 'x'
                badAxisVec = repmat('x', [2*Npairs, 1]);
            case{'alt', 'xy', 'alternate'}
                % Change probe ordering for odd- vs even-numbered
                % WFSC iterations.
                if mod(Itr, 2) == 1
                    badAxisVec = repmat('x', [2*Npairs, 1]);
                    badAxisVec(3:4:end) = 'y';
                    badAxisVec(4:4:end) = 'y';
                else
                    badAxisVec = repmat('y', [2*Npairs, 1]);
                    badAxisVec(3:4:end) = 'x';
                    badAxisVec(4:4:end) = 'x';
                end

            case 'multi'
                badAxisVec = repmat('m', [2*Npairs, 1]);
        end
end


%% Get images and perform estimates in each sub-bandpass

fprintf('Estimating electric field with batch process estimation ...\n'); tic;

for iStar = 1:mp.compact.star.count

    modvar = ModelVariables;
    modvar.starIndex = iStar;
    modvar.whichSource = 'star';

for iSubband = 1:mp.Nsbp

    modvar.sbpIndex = iSubband;
    fprintf('Wavelength: %u/%u ... ', iSubband, mp.Nsbp);
    modeIndex = (iStar-1)*mp.Nsbp + iSubband;
    fprintf('Mode: %u/%u ... ', modeIndex, mp.jac.Nmode);    

    %% Measure current contrast level average, and on each side of Image Plane
    % Reset DM commands to the unprobed state:
    if any(mp.dm_ind == 1); mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vnom); end
    if any(mp.dm_ind == 2); mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vnom); end
    
    %% Separate out values of images at dark hole pixels and delta DM voltage settings
    
    Iplus  = zeros([mp.Fend.corr.Npix, Npairs]); % Pixels of plus probes' intensities
    Iminus = zeros([mp.Fend.corr.Npix, Npairs]); % Pixels of minus probes' intensities
    if mp.flagFiber
        Ifiberplus  = zeros([mp.Fend.Nfiber, Npairs]); % Pixels of plus probes' intensities
        Ifiberminus = zeros([mp.Fend.Nfiber, Npairs]); % Pixels of minus probes' intensities
    end
    DM1Vplus  = zeros([Nact, Nact, Npairs]);
    DM1Vminus = zeros([Nact, Nact, Npairs]);
    DM2Vplus  = zeros([Nact, Nact, Npairs]);
    DM2Vminus = zeros([Nact, Nact, Npairs]);

    %% Compute probe shapes and take probed images:

    %--Take initial, unprobed image (for unprobed DM settings).
    whichImage = 1;
    mp.isProbing = false; % tells the camera whether to use the exposure time for either probed or unprobed images.
    if ~mp.flagFiber
        I0 = falco_get_sbp_image(mp, iSubband);
        ev.score.Inorm = mean(I0(mp.Fend.score.maskBool));
        ev.corr.Inorm  = mean(I0(mp.Fend.corr.maskBool));
    else
        [I0,I0fiber] = falco_get_sbp_image(mp, iSubband);
        I0fibervec = I0fiber;
        ev.score.Inorm = mean(I0fibervec);
        ev.corr.Inorm  = mean(I0fibervec);
    end
    I0vec = I0(mp.Fend.corr.maskBool); % Vectorize the correction region pixels
    
    if iStar == 1 % Image already includes all stars, so don't sum over star loop
        ev.Im = ev.Im + mp.sbp_weights(iSubband)*I0; % subband-averaged image for plotting
        if mp.flagFiber; ev.Ifiber = ev.Ifiber + mp.sbp_weights(iSubband)*I0fiber;end % 

        %--Store values for first image and its DM commands
        ev.imageArray(:, :, whichImage, iSubband) = I0;
        if any(mp.dm_ind == 1);  ev.dm1.Vall(:, :, whichImage, iSubband) = mp.dm1.V;  end
        if any(mp.dm_ind == 2);  ev.dm2.Vall(:, :, whichImage, iSubband) = mp.dm2.V;  end
    end
    
    %--Compute the average Inorm in the scoring and correction regions
    fprintf('Measured unprobed Inorm (Corr / Score): %.2e \t%.2e \n',ev.corr.Inorm,ev.score.Inorm);    

    % Set (approximate) probe intensity based on current measured Inorm
    if isempty(mp.est.probeSchedule.InormProbeVec)
        ev.InormProbeMax = mp.est.InormProbeMax;
        if mp.flagFiber
            InormProbe = min([sqrt(max(I0fibervec)*1e-8), ev.InormProbeMax]);
        else
            InormProbe = min([sqrt(max(I0vec)*1e-5), ev.InormProbeMax]);
        end
        fprintf('Chosen probe intensity: %.2e \n', InormProbe);
    else
        InormProbe = mp.est.probeSchedule.InormProbeVec(Itr);
        fprintf('Scheduled probe intensity: %.2e \n', InormProbe);
    end

    %--Perform the probing
    mp.isProbing = true; % tells the camera whether to use the exposure time for either probed or unprobed images.
    iOdd = 1; iEven = 1; % Initialize index counters
    for iProbe = 1:2*Npairs           

        %--Generate the DM command map for the probe
        switch lower(mp.estimator)
            case{'pairwise-rect', 'pwp-bp', 'pairwise-rect-kf'} 
                probeCmd = falco_gen_pairwise_probe(mp, InormProbe, probePhaseVec(iProbe), iStar, mp.est.probe.rotation);
            case{'pairwise', 'pairwise-square', 'pwp-bp-square', 'pwp-kf', 'pairwise-kf'}
                probeCmd = falco_gen_pairwise_probe_square(mp, InormProbe, probePhaseVec(iProbe), badAxisVec(iProbe), mp.est.probe.rotation);
        end
        %--Select which DM to use for probing. Allocate probe to that DM
        if whichDM == 1
            dDM1Vprobe = probeCmd ./ mp.dm1.VtoH; % Now in volts
            dDM2Vprobe = 0;
        elseif whichDM == 2
            dDM1Vprobe = 0;        
            dDM2Vprobe = probeCmd ./ mp.dm2.VtoH; % Now in volts
        end

        if any(mp.dm_ind == 1)
            mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vnom + dDM1Vprobe); 
        end
        if any(mp.dm_ind == 2)
            mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vnom + dDM2Vprobe);
        end

        %--Take probed image
        if mp.flagFiber
            [Im,Ifiber] = falco_get_sbp_image(mp, iSubband);
            ev.IprobedMean = ev.IprobedMean + Ifiber/(2*Npairs); %--Inorm averaged over all the probed images
        else
            Im = falco_get_sbp_image(mp, iSubband);
            ev.IprobedMean = ev.IprobedMean + mean(Im(mp.Fend.corr.maskBool))/(2*Npairs); %--Inorm averaged over all the probed images
        end
        whichImage = 1+iProbe; %--Increment image counter

        %--Store probed image and its DM settings
        ev.imageArray(:, :, whichImage, iSubband) = Im;
        if any(mp.dm_ind == 1);  ev.dm1.Vall(:, :, whichImage, iSubband) = mp.dm1.V;  end
        if any(mp.dm_ind == 2);  ev.dm2.Vall(:, :, whichImage, iSubband) = mp.dm2.V;  end

        %--Report results
        probeSign = ['-', '+'];
        if mp.flagFiber
            fprintf('Actual Probe %d%s Contrast is: %.2e \n', ceil(iProbe/2), probeSign(mod(iProbe, 2)+1), mean(Ifiber));
        else
            fprintf('Actual Probe %d%s Contrast is: %.2e \n', ceil(iProbe/2), probeSign(mod(iProbe, 2)+1), mean(Im(mp.Fend.corr.maskBool)));
        end
        
        %--Assign image to positive or negative probe collection:
        if mod(iProbe, 2) == 1  % Odd; for plus probes
            if whichDM == 1;  DM1Vplus(:, :, iOdd) = dDM1Vprobe + DM1Vnom;  end
            if whichDM == 2;  DM2Vplus(:, :, iOdd) = dDM2Vprobe + DM2Vnom;  end
            Iplus(:, iOdd) = Im(mp.Fend.corr.maskBool);
            if mp.flagFiber
                Ifiberplus(:, iOdd) = Ifiber;
            end

            iOdd = iOdd + 1;
        elseif mod(iProbe, 2) == 0  % Even; for minus probes
            if whichDM == 1;  DM1Vminus(:, :, iEven) = dDM1Vprobe + DM1Vnom;  end
            if whichDM == 2;  DM2Vminus(:, :, iEven) = dDM2Vprobe + DM2Vnom;  end 
            Iminus(:, iEven) = Im(mp.Fend.corr.maskBool);
            if mp.flagFiber
                Ifiberminus(:, iEven) = Ifiber;
            end
            iEven = iEven + 1;
        end
    end

    %% Calculate probe amplitudes and measurement vector. (Refer again to Give'on+ SPIE 2011 to undersand why.)
    ampSq = (Iplus+Iminus)/2 - repmat(I0vec, [1,Npairs]);  % square of probe E-field amplitudes
    ampSq(ampSq < 0) = 0;  % If probe amplitude is zero, amplitude is zero there.
    if ~mp.flagFiber
        amp = sqrt(ampSq);   % E-field amplitudes, dimensions: [mp.Fend.corr.Npix, Npairs]
        isnonzero = all(amp, 2);
        zAll = ((Iplus-Iminus)/4).';  % Measurement vector, dimensions: [Npairs, mp.Fend.corr.Npix]
    else
        ampSqFiber = (Ifiberplus+Ifiberminus)/2 - repmat(I0fibervec, [1,Npairs]);  % square of probe E-field amplitudes
        ampSqFiber(ampSqFiber < 0) = 0;  % If probe amplitude is zero, amplitude is zero there.
        amp = sqrt(ampSqFiber);   % E-field amplitudes, dimensions: [mp.Fend.corr.Npix, Npairs]
        isnonzero = all(amp, 2);
        zAll = ((Ifiberplus-Ifiberminus)/4).';  % Measurement vector, dimensions: [Npairs, mp.Fend.corr.Npix]
    end
    ampSq2Dcube = zeros(mp.Fend.Neta, mp.Fend.Nxi, mp.est.probe.Npairs);
    for iProbe=1:Npairs % Display the actual probe intensity
        ampSq2D = zeros(mp.Fend.Neta, mp.Fend.Nxi); ampSq2D(mp.Fend.corr.maskBool) = ampSq(:, iProbe); 
        ampSq2Dcube(:, :, iProbe) = ampSq2D;
        fprintf('*** Mean measured Inorm for probe #%d  =\t%.3e \n',iProbe,mean(ampSq2D(mp.Fend.corr.maskBool)));
    end

    %% Plot relevant data for all the probes
    ev.iStar = iStar;
    if whichDM == 1
        DMV4plot = DM1Vplus - repmat(DM1Vnom, [1, 1, size(DM1Vplus, 3)]);
    elseif whichDM == 2
        DMV4plot = DM2Vplus - repmat(DM2Vnom, [1, 1, size(DM2Vplus, 3)]);
    end
    falco_plot_pairwise_probes(mp, ev, DMV4plot, ampSq2Dcube, iSubband)

    %% Perform the estimation
    
    if mp.est.flagUseJac %--Use Jacobian for estimation. This is fully model-based if the Jacobian is purely model-based, or it is better if the Jacobian is adaptive based on empirical data.
        
        dEplus  = zeros(size(Iplus ));
        for iProbe=1:Npairs
            if whichDM == 1
                dV = DM1Vplus(:, :, iProbe) - DM1Vnom;
                dEplus(:, iProbe) = squeeze(jacStruct.G1(:, :, modeIndex))*dV(mp.dm1.act_ele);
            elseif whichDM == 2
                dV = DM2Vplus(:, :, iProbe) - DM2Vnom;
                dEplus(:, iProbe) = squeeze(jacStruct.G2(:, :, modeIndex))*dV(mp.dm2.act_ele);
            end
        end
        
    else %--Get the probe phase from the model and the probe amplitude from the measurements

        % For unprobed field based on model:
        if any(mp.dm_ind == 1); mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vnom); end
        if any(mp.dm_ind == 2); mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vnom); end

        if mp.flagFiber
            [~, E0] = model_compact(mp, modvar);
            E0vec = E0;
            %--For probed fields based on model:
            Eplus  = zeros(size(Ifiberplus ));
            Eminus = zeros(size(Ifiberminus));
        else
            E0 = model_compact(mp, modvar);
            E0vec = E0(mp.Fend.corr.maskBool);
             %--For probed fields based on model:
            Eplus  = zeros(size(Iplus ));
            Eminus = zeros(size(Iminus));
        end
        
        for iProbe=1:Npairs
            % For plus probes:
            if whichDM == 1
                mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vplus(:, :, iProbe));
            elseif whichDM == 2
                mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vplus(:, :, iProbe));
            end
            if(mp.flagFiber)
                [~, Etemp] = model_compact(mp, modvar);
                Eplus(:, iProbe) = Etemp;
            else
                Etemp = model_compact(mp, modvar);
                Eplus(:, iProbe) = Etemp(mp.Fend.corr.maskBool);
            end
            
            % For minus probes:
            if whichDM == 1
                mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vminus(:, :, iProbe));
            elseif whichDM == 2
                mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vminus(:, :, iProbe));
            end
            if mp.flagFiber
                [~, Etemp] = model_compact(mp, modvar);
                Eminus(:, iProbe) = Etemp;
            else
                Etemp = model_compact(mp, modvar);
                Eminus(:, iProbe) = Etemp(mp.Fend.corr.maskBool);
            end
            
        end

        %%--Create delta E-fields for each probe image. Then create Npairs phase angles.
        dEplus  = Eplus  - repmat(E0vec, [1, Npairs]);
        dEminus = Eminus - repmat(E0vec, [1, Npairs]);
        dEprobe = (dEplus - dEminus)/2; % Take the average to mitigate nonlinearity
        [dphdm, amp_model] = deal(zeros([Ncorr, Npairs])); %--phases and model amp of the probes
        for iProbe = 1:Npairs
            dphdm(:, iProbe) = atan2(imag(dEprobe(:, iProbe)), real(dEprobe(:, iProbe)));
            % model predicted probe amplitude, only for diagnostics
            amp_model(:, iProbe) = abs(dEprobe(:, iProbe));
        end
        
    end 
    
    %% Batch process the measurements to estimate the electric field in the dark hole. Done pixel by pixel.

    if useKalmanFilter
        Hall = zeros(Npairs, 2, mp.Fend.corr.Npix);
    end
    Eest = zeros(Ncorr, 1);
    zerosCounter = 0;
    partialCounter = 0;
    for ipix = 1:Ncorr
        
        H = zeros(Npairs, 2); % Initialize the observation matrix
        
        if mp.est.flagUseJac 
            dE = dEplus(ipix, :).';
            H = [real(dE), imag(dE)];
            Epix = pinv(H)*zAll(:, ipix); %--Batch process estimation
        else
            % H = zeros(Npairs, 2); % Observation matrix
            goodProbeInds = find(amp(ipix, :) > 0);
            NpairsGood = length(goodProbeInds);
            
            % If <2 probe pairs had good measurements, can't do pinv. Leave Eest as zero.
            if NpairsGood < 2
                zerosCounter = zerosCounter + 1;
                Epix = [0; 0];
            
            % Otherwise, use the 2+ good probe pair measurements for that pixel:
            else
                for iProbe = 1:Npairs
                    H(iProbe, :) = amp(ipix, iProbe)*[cos(dphdm(ipix, iProbe)), sin(dphdm(ipix, iProbe)) ];
                end
                % Epix = pinv(H)*zAll(:, ipix); %--Batch process estimation
                Epix = pinv(H(goodProbeInds, :)) * zAll(goodProbeInds, ipix); %--Batch process estimation
                
                % Record how many pixels didn't use all the probe pairs
                if NpairsGood < Npairs
                    partialCounter = partialCounter + 1;
                end
                
            end

        end
        
        Eest(ipix) = Epix(1) + 1i*Epix(2);
        if useKalmanFilter
            Hall(:, :, ipix) = H;
        end
        
    end

    Eest(abs(Eest).^2 > mp.est.Ithreshold) = 0;  % If estimate is too bright, the estimate was probably bad. !!!!!!!!!!!!!!BE VERY CAREFUL WITH THIS VALUE!!!!!!!!!!!!!!!
    fprintf('%d of %d pixels were given zero probe amplitude. \n', zerosCounter, mp.Fend.corr.Npix);
    if Npairs > 2 % Only possible for Npairs > 2
        fprintf('%d of %d pixels threw out at least one probe pair. \n', partialCounter, mp.Fend.corr.Npix); 
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Kalman Filter Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    if useKalmanFilter

        if ev.Itr < mp.est.kf.ItrStart

            % Re-organize the batch-processed E-field estimate into the first state estimate for the Kalman filter
            xOld = zeros(2*Ncorr, 1);
            for ii = 1:Ncorr
                xOld(2*(ii-1)+1:2*(ii-1)+2) = [real(Eest(ii)); imag(Eest(ii))];
            end
            ev.xOld = xOld;

            %--Initialize the state covariance matrix (2x2 for each dark hole pixel)
            ev.Pold_KF_array = repmat(mp.est.kf.Pcoef0*eye(2), [mp.Fend.corr.Npix, 1, mp.Nsbp]);

        else

            xOld = ev.xOld;
            Pold = ev.Pold_KF_array(:, :, modeIndex);

            %--Compute the change in E-field since last correction iteration
            if mp.est.flagUseJac %--Use Jacobian to compute delta E-field from previous correction step to now.
                if whichDM == 1
                    dE = squeeze(jacStruct.G1(:, :, modeIndex)) * mp.dm1.dV(mp.dm1.act_ele);
                elseif whichDM == 2
                    dE = squeeze(jacStruct.G2(:, :, modeIndex)) * mp.dm2.dV(mp.dm2.act_ele);
                end
            else
                % For Xminus, use nonlinear dynamics instead of Gamma. This means
                % difference the output of model_compact rather than using the
                % Jacobian.
                %--Previous unprobed field based on model:
                if whichDM == 1;  mp.dm1 = falco_set_constrained_voltage(mp.dm1, DM1Vnom - mp.dm1.dV); end
                if whichDM == 2;  mp.dm2 = falco_set_constrained_voltage(mp.dm2, DM2Vnom - mp.dm2.dV); end
                if mp.flagFiber
                    [~, Eprev] = model_compact(mp, modvar);
                else
                    Eprev = model_compact(mp, modvar);
                end
                EprevVec = Eprev(mp.Fend.corr.maskBool);
                dE = E0vec - EprevVec; % Change in unprobed E-field between correction iterations
            end

            %--Construct dX, the change in state, from dE
            dX = zeros(size(xOld));
            for ii = 1:mp.Fend.corr.Npix
               dX(2*(ii-1)+1:2*(ii-1)+2) = [real(dE(ii)); imag(dE(ii))];
            end

            %--Compute Sensor Noise, R. You can calculate this from the properties
            % of your optical system.
            % ncounts_readout = mp.readNoiseStd;  % Read noise in counts (ADU). Might need to be in photo-electrons
            % ncounts_shot = sqrt(ev.IprobedMean*mp.peakCountsPerPixPerSec);
            % Dark current not included here (yet).
            ncounts_peak = mp.detector.peakFluxVec(iSubband) * mp.detector.tExpProbedVec(iSubband); % Using exposure time for probes
            ncounts_std = sqrt( (sqrt(2)*ev.IprobedMeanSubband(iSubband)*ncounts_peak + mp.detector.readNoiseStd^2)/mp.detector.Nexp);
            Rvar = (ncounts_std/ncounts_peak)^2; % Don't forget to square it since R = E<n*n.'>. This is a variable scalar
            Rmat = mp.est.kf.Rcoef*Rvar*eye(Npairs);
            fprintf('Sensor noise coefficient: %.3e\n', Rmat(1, 1));

            %--Compute Process Noise, Q. This you just have to tune (via the scalar mp.est.kf.Qcoef).
            Q00 = ev.corr.Inorm; %--Set coefficient to the current measured contrast. Can change to the last estimated contrast.
            Q = Q00*mp.est.kf.Qcoef*repmat( eye(2), [mp.Fend.corr.Npix, 1]);
            dP = Q;
            fprintf('Process noise coefficient: %.3e\n', Q(1, 1));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute Kalman Filter Equations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xNew = zeros(2*mp.Fend.corr.Npix, 1);
            Pnew = zeros(2*mp.Fend.corr.Npix, 2);
            Kall = zeros(2, Npairs, mp.Fend.corr.Npix); %--Store the Kalman gain matrices. For diagnostics only

            xMinusAll = xOld + dX; %State Estimate Extrapolation (Recall Phi = Identity)
            PminusAll = Pold + dP; %Covariance Estimate Extrapolation (Recall Phi = Identity)
            for ii=1:mp.Fend.corr.Npix
                %--Format several matrices
                z = zAll(:, ii); % z has dimensions [Npairs x Npix]
                H = Hall(:, :, ii); % H has dimensions [Npairs x 2]
                xMinus = xMinusAll(2*(ii-1)+1:2*(ii-1)+2); % x has dimensions [2 x 1]
                Pminus = PminusAll(2*(ii-1)+1:2*(ii-1)+2, 1:2); % P has dimensions [2 x 2]

                %--Compute the Kalman gain matrix.
                K = Pminus*H.'/(H*Pminus*H.' + Rmat); % dimensions [2 x Npairs]
                Kall(:, :, ii) = K; % Store Kalman gains for diagnostics only

                %--State Estimate Update
                xNew(2*(ii-1)+1:2*(ii-1)+2) = xMinus + K*(z-H*xMinus);
                Pnew(2*(ii-1)+1:2*(ii-1)+2, 1:2) = (eye(2)-K*H)*Pminus;
            end

            %--For diagnostics only. Print the mean state covariance matrix to the command line
            P11 = mean(Pnew(1:2:end, 1));
            P22 = mean(Pnew(2:2:end, 2));
            P12 = mean(Pnew(1:2:end, 2));
            P21 = mean(Pnew(2:2:end, 1));
            fprintf('P11,P22,P12,P21: \t%.2e \t%.2e \t%.2e \t%.2e \n', P11, P22, P12, P21);

            % Re-order Xnew into Eest
            EestKF = zeros(mp.Fend.corr.Npix, 1); %for re-ordering state for control
            for jj=1:mp.Fend.corr.Npix
                EestKF(jj) = xNew(2*(jj-1)+1) + 1i*xNew(2*(jj-1)+2);
            end
            Eest = EestKF;

            %--Save out new state covariance P as the old for the next correction iteration.
            ev.Pold_KF_array(:, :, modeIndex) = Pnew;

            fprintf(['mean(abs(K)) = ' num2str(mean(mean(mean(abs(Kall))))) '.\n']);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % End Kalman Filter Computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save out the estimates
ev.Eest(:, modeIndex) = Eest;
if mp.flagFiber
    ev.IincoEst(:, modeIndex) =  I0fibervec - abs(Eest).^2; % incoherent light
else
    ev.IincoEst(:, modeIndex) =  I0vec - abs(Eest).^2; % incoherent light
end

if mp.flagPlot && ~mp.flagFiber
    Eest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Eest2D(mp.Fend.corr.maskBool) = Eest;
    %figure(701); imagesc(real(Eest2D)); title('real(Eest)', 'Fontsize', 18); set(gca, 'Fontsize', 18); axis xy equal tight; colorbar;
    %figure(702); imagesc(imag(Eest2D)); title('imag(Eest)', 'Fontsize', 18); set(gca, 'Fontsize', 18); axis xy equal tight; colorbar;
    %figure(703); imagesc(log10(abs(Eest2D).^2)); title('abs(Eest)^2', 'Fontsize', 18); set(gca, 'Fontsize', 18); axis xy equal tight; colorbar;
    %drawnow;
end

end %--End of loop over the wavelengths
end %--End of loop over stars

%--Other data to save out
ev.ampSqMean = mean(ampSq(:)); %--Mean probe intensity
ev.ampNorm = amp/sqrt(InormProbe); %--Normalized probe amplitude maps
ev.InormProbe = InormProbe;        
ev.maskBool = mp.Fend.corr.maskBool; %--for resizing Eest and IincoEst for plotting
ev.amp_model = amp_model;

mp.isProbing = false; % tells the camera whether to use the exposure time for either probed or unprobed images.

fprintf(' done. Time: %.3f\n',toc);

end %--END OF FUNCTION