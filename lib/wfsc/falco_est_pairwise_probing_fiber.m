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

function [ev] = falco_est_pairwise_probing_fiber(mp,varargin)

%--If there is a second input, it is the Jacobian structure
if(size(varargin, 2)==1)
    jacStruct = varargin{1};
end

%--Select number of actuators across based on chosen DM for the probing
if(mp.est.probe.whichDM==1)
    Nact = mp.dm1.Nact;
elseif(mp.est.probe.whichDM==2)
    Nact = mp.dm2.Nact;
end

%--Store the initial DM commands
if(any(mp.dm_ind==1))
    DM1Vnom = mp.dm1.V;  
end
if(any(mp.dm_ind==2))
    DM2Vnom = mp.dm2.V;  
else
    DM2Vnom = zeros(size(mp.dm1.V)); 
end

% Definitions:
Npairs = mp.est.probe.Npairs; % % Number of image PAIRS for DM Diversity or Kalman filter initialization
ev.Icube = zeros(mp.Fend.Neta,mp.Fend.Nxi,1+2*Npairs);
if(any(mp.dm_ind==1))
    ev.Vcube.dm1 = zeros(mp.dm1.Nact,mp.dm1.Nact,1+2*Npairs);
end
if(any(mp.dm_ind==2))
    ev.Vcube.dm2 = zeros(mp.dm2.Nact,mp.dm2.Nact,1+2*Npairs);
end

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
if(mp.flagLenslet)
    ev.Eest = zeros(mp.Fend.Nlens, mp.Nsbp); %Final estimated E-field
    ev.IincoEst = zeros(mp.Fend.Nlens, mp.Nsbp); %Estimated incoming intensity
else
    ev.Est = zeros(mp.Fend.Nfiber, mp.Nsbp);
    ev.IincoEst = zeros(mp.Fend.Nfiber, mp.Nsbp);
end
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
    if(mp.flagLenslet)
        Iplus  = zeros([mp.Fend.Nlens, Npairs]); % Pixels of plus probes' intensities
        Iminus = zeros([mp.Fend.Nlens, Npairs]); % Pixels of minus probes' intensities
    else
        Iplus = zeros([mp.Fend.Nfiber, Npairs]);
        Iminus = zeros([mp.Fend.Nfiber, Npairs]);
    end
    DM1Vplus  = zeros([Nact,Nact, Npairs]);
    DM1Vminus = zeros([Nact,Nact, Npairs]);
    DM2Vplus  = zeros([Nact,Nact, Npairs]);
    DM2Vminus = zeros([Nact,Nact, Npairs]);

    %% Compute probe shapes and take probed images:

    %--Take initial, unprobed image (for unprobed DM settings).
    whichImg = 1;
    I0 = max(max(falco_get_sbp_image_fiber(mp,si)));
    ev.I0mean = ev.I0mean+I0/mp.Nsbp; %--Getting the sub-bandpass-averaged Inorm

    %--Store values for first image and its DM commands
    ev.Icube(:,:, whichImg) = I0;
    if(any(mp.dm_ind==1))
        ev.Vcube.dm1(:,:,whichImg) = mp.dm1.V;
    end
    if(any(mp.dm_ind==2))
        ev.Vcube.dm2(:,:,whichImg) = mp.dm2.V; 
    end

    %--Compute the average Inorm in the scoring and correction regions
    ev.Inorm = mean(I0);
    fprintf('Measured unprobed Ifiber : %.2e \n',ev.Inorm);    

    % Set (approximate) probe intensity based on current measured Inorm
    ev.InormProbeMax = 1e-6;
    InormProbe = min([sqrt(max(I0)*1e-12), ev.InormProbeMax]); %--Change this to a high percentile value (e.g., 90%) instead of the max to avoid being tricked by noise
    fprintf('Chosen probe intensity: %.2e \n',InormProbe);
    
    %--Perform the probing
    iOdd=1; iEven=1; %--Initialize index counters
    for iProbe=1:2*Npairs

        %--Generate the command map for the probe
        probeCmd = falco_gen_pairwise_probe_fiber(mp,InormProbe,probePhaseVec(iProbe),badAxisVec(iProbe));

        figure(901);
        imagesc(probeCmd); axis equal tight;
        
        %--Select which DM to use for probing. Allocate probe to that DM
        if(mp.est.probe.whichDM == 1)
            dDM1Vprobe = probeCmd./mp.dm1.VtoH; % Now in volts
            dDM2Vprobe = 0;
        elseif(mp.est.probe.whichDM == 2)
            dDM1Vprobe = 0;        
            dDM2Vprobe = probeCmd./mp.dm1.VtoH; % Now in volts
        end
        if(any(mp.dm_ind==1))
            mp.dm1.V = DM1Vnom+dDM1Vprobe;
        end
        if(any(mp.dm_ind==2))
            mp.dm2.V = DM2Vnom+dDM2Vprobe;
        end

        %--Take probed image
        Im = max(max(falco_get_sbp_image_fiber(mp,si)));
        whichImg = 1+iProbe; %--Increment image counter
        ev.IprobedMean = ev.IprobedMean + mean(Im)/(2*Npairs); %--Inorm averaged over all the probed images

        %--Store probed image and its DM settings
        ev.Icube(:,whichImg) = Im;
        if(any(mp.dm_ind==1))
            ev.Vcube.dm1(:,:,whichImg) = mp.dm1.V;
        end
        if(any(mp.dm_ind==2))
            ev.Vcube.dm2(:,:,whichImg) = mp.dm2.V;
        end

        %--Report results
        probeSign = ['-','+'];
        fprintf('Actual Probe %d%s Contrast is: %.2e \n',ceil(iProbe/2),probeSign(mod(iProbe,2)+1),mean(Im));

        %--Assign image to positive or negative probe collection:
        if mod(iProbe,2)==1  % Odd; for plus probes
            if(any(mp.dm_ind==1));  DM1Vplus(:,:,iOdd) = dDM1Vprobe + DM1Vnom;  end
            if(any(mp.dm_ind==2));  DM2Vplus(:,:,iOdd) = dDM2Vprobe + DM2Vnom;  end
            Iplus(:,iOdd) = Im;
            iOdd=iOdd+1;
        elseif mod(iProbe,2)==0  % Even; for minus probes
            if(any(mp.dm_ind==1));  DM1Vminus(:,:,iEven) = dDM1Vprobe + DM1Vnom;  end
            if(any(mp.dm_ind==2));  DM2Vminus(:,:,iEven) = dDM2Vprobe + DM2Vnom;  end 
            Iminus(:,iEven) = Im;
            iEven=iEven+1;      
        end
    end

    %% Calculate probe amplitudes and measurement vector. (Refer again to Give'on+ SPIE 2011 to undersand why.)
    ampSq = (Iplus+Iminus)/2 - repmat(I0,[1,Npairs]);  % square of probe E-field amplitudes
    ampSq(ampSq<0) = 0;  % If probe amplitude is zero, amplitude is zero there.
    amp = sqrt(ampSq);   % E-field amplitudes, dimensions: [mp.Fend.corr.Npix, Npairs]
    isnonzero = all(amp,2);
    zAll = ((Iplus-Iminus)/4).';  % Measurement vector, dimensions: [Npairs,mp.Fend.Nlens]
    
    %% Perform the estimation
    
    if(mp.est.flagUseJac) %--Use Jacobian for estimation. This is fully model-based if the Jacobian is purely model-based, or it is better if the Jacobian is adaptive based on empirical data.
        
        dEplus  = zeros(size(Iplus));
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
        if(any(mp.dm_ind==1))
            mp.dm1.V = DM1Vnom;
        end
        if(any(mp.dm_ind==2))
            mp.dm2.V = DM2Vnom;
        end
        
        [~,E0] = max(max(model_compact(mp, modvar)));
        
        %--For probed fields based on model:
        Eplus  = zeros(size(Iplus));
        Eminus = zeros(size(Iminus));
        for iProbe=1:Npairs
            % For plus probes:
            if(any(mp.dm_ind==1))
                mp.dm1.V = squeeze(DM1Vplus(:,:,iProbe));
            end
            if(any(mp.dm_ind==2))
                mp.dm2.V = squeeze(DM2Vplus(:,:,iProbe));
            end
            [~,Etemp] = max(max(model_compact(mp, modvar)));
            Eplus(:,iProbe) = Etemp;
            % For minus probes:
            if(any(mp.dm_ind==1))
                mp.dm1.V = squeeze(DM1Vminus(:,:,iProbe));
            end
            if(any(mp.dm_ind==2))
                mp.dm2.V = squeeze(DM2Vminus(:,:,iProbe));
            end
            [~,Etemp] = max(max(model_compact(mp, modvar)));
            Eminus(:,iProbe) = Etemp;
        end

        %%--Create delta E-fields for each probe image. Then create Npairs phase angles.
        dEplus  = Eplus  - repmat(E0,[1,Npairs]);
        dEminus = Eminus - repmat(E0,[1,Npairs]);
        if(mp.flagLenslet)
            dphdm = zeros([mp.Fend.Nlens, Npairs]); %--phases of the probes
        else
            dphdm = zeros([mp.Fend.Nfiber, Npairs]);
        end
        for iProbe=1:Npairs
            dphdm(:,iProbe) = atan2(imag(dEplus(:,iProbe))-imag(dEminus(:,iProbe)),real(dEplus(:,iProbe))-real(dEminus(:,iProbe)));
        end
        
    end 
    
%% Batch process the measurements to estimate the electric field in the dark hole. Done pixel by pixel.

if(strcmpi(mp.estimator,'pwp-bp'))
    if(mp.flagLenslet)
        Eest = zeros(mp.Fend.Nlens,1);
        loopend = mp.Fend.Nlens;
    else
        Eest = zeros(mp.Fend.Nfiber,1);
        loopend = mp.Fend.Nfiber;
    end
    zerosCounter = 0;
    for ipix=1:loopend
        
        if(mp.est.flagUseJac) 
            dE = dEplus(ipix,:).';
            H = [real(dE),imag(dE)];
        else
            H = zeros(Npairs,2); % Observation matrix
            if(isnonzero(ipix)==1) % Leave Eest for a pixel as zero if any probe amplitude is zero there.
                for iProbe=1:Npairs
                    H(iProbe,:) = amp(ipix,iProbe)*[cos(dphdm(ipix,iProbe)), sin(dphdm(ipix,iProbe))];
                end
            else
                zerosCounter = zerosCounter+1;
            end
        end
        Epix = pinv(H)*zAll(:,ipix); %--Batch process estimation
        Eest(ipix) = Epix(1) + 1i*Epix(2);
    end
    Eest(abs(Eest).^2 > 1e-4) = 0;  % If estimate is too bright, the estimate was probably bad. !!!!!!!!!!!!!!BE VERY CAREFUL WITH THIS HARD-CODED VALUE!!!!!!!!!!!!!!!
    fprintf('%d of %d pixels were given zero probe amplitude. \n',zerosCounter,mp.Fend.corr.Npix); 
    Eest(abs(Eest).^2 < 0) = 1e-12;  % If estimate is too bright, the estimate was probably bad. !!!!!!!!!!!!!!BE VERY CAREFUL WITH THIS HARD-CODED VALUE!!!!!!!!!!!!!!!

end   

%% Save out the estimates
ev.Eest(:,si) = Eest;

ev.IincoEst(:,si) =  I0-abs(Eest).^2; %--Compute the incoherent light

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