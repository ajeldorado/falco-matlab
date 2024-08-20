function ev = falco_est_modal_ekf_maintenance(mp, ev, varargin)

%% This stuff has been copy-pasted

Itr = ev.Itr;

whichDM = mp.est.probe.whichDM;

if ~isa(mp.est.probe, 'Probe')
    error('mp.est.probe must be an instance of class Probe')
end
%--If there is a third input, it is the Jacobian structure
if size(varargin, 2) == 1
    jacStruct = varargin{1};
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

% Initialize output arrays
ev.imageArrayPlus = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);
ev.imageArrayMinus = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);
ev.Eest = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IincoEst = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IprobedMean = 0;
ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
if whichDM == 1;  ev.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, 1, mp.Nsbp);  end
if whichDM == 2;  ev.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, 1, mp.Nsbp);  end

%% Get dither command

if any(mp.dm_ind == 1);  DM1Vdither = normrnd(0,mp.est.dither,[mp.dm1.Nact mp.dm1.Nact]); else; DM1Vdither = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2);  DM2Vdither = normrnd(0,mp.est.dither,[mp.dm1.Nact mp.dm1.Nact]); else; DM2Vdither = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1

ev.dither = get_dm_command_vector(mp,DM1Vdither, DM2Vdither);

%% Set total command for estimator image
% TODO: need to save these commands for each iteration separately

if Itr > 1
    if ~isfield(mp.dm1,'dV'); mp.dm1.dV = zeros(mp.dm1.Nact);end
    if ~isfield(mp.dm2,'dV'); mp.dm2.dV = zeros(mp.dm2.Nact);end
    efc_command = get_dm_command_vector(mp,mp.dm1.dV, mp.dm2.dV);
else
    efc_command = 0*ev.dither;
    mp.dm1.dV = zeros(size(DM1Vdither));
    mp.dm2.dV = zeros(size(DM1Vdither));
end

% Generate command to apply to DMs Plus
% Note if dm_drift_ind ~= i, the command is set to zero in
% falco_drift_injection
if any(mp.dm_ind == 1)
    % note falco_set_constrained_voltage does not apply the command to the
    % DM
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift + mp.dm1.dV + DM1Vdither + mp.dm1.V_shift); 
elseif any(mp.dm_drift_ind == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift); 
elseif any(mp.dm_ind_static == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz);
end

if any(mp.dm_ind == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift + mp.dm2.dV + DM2Vdither + mp.dm2.V_shift); 
elseif any(mp.dm_drift_ind == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift); 
elseif any(mp.dm_ind_static == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz);
end

% Do safety check to make sure no actuators are pinned
ev = pinned_act_safety_check(mp,ev);

% closed_loop_command_plus = efc_command + get_dm_command_vector(mp,mp.dm1.V_shift, mp.dm2.V_shift) + dither;
% closed_loop_command_minus = efc_command + get_dm_command_vector(mp,mp.dm1.V_shift, mp.dm2.V_shift) - dither;

% Get images

for iSubband = 1:mp.Nsbp
    ev.imageArrayPlus(:,:,1,iSubband) = falco_get_sbp_image(mp, iSubband);
end

% Generate command to apply to DMs Minus
if any(mp.dm_ind == 1)
    % note falco_set_constrained_voltage does not apply the command to the
    % DM
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift + mp.dm1.dV - DM1Vdither + mp.dm1.V_shift); 
elseif any(mp.dm_drift_ind == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift); 
elseif any(mp.dm_ind_static == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz);
end

if any(mp.dm_ind == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift + mp.dm2.dV - DM2Vdither + mp.dm2.V_shift); 
elseif any(mp.dm_drift_ind == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift); 
elseif any(mp.dm_ind_static == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz);
end

% Do safety check to make sure no actuators are pinned
ev = pinned_act_safety_check(mp,ev);

closed_loop_command = efc_command + get_dm_command_vector(mp,mp.dm1.V_shift, mp.dm2.V_shift);
% Get images

y_measured = zeros(mp.Fend.corr.Npix,mp.Nsbp);
for iSubband = 1:mp.Nsbp
    ev.imageArrayMinus(:,:,1,iSubband) = falco_get_sbp_image(mp, iSubband);
end

for iSubband = 1:mp.Nsbp
    I0plus = ev.imageArrayPlus(:,:,1,iSubband) * ev.peak_psf_counts(iSubband);
    I0minus = ev.imageArrayMinus(:,:,1,iSubband) * ev.peak_psf_counts(iSubband);
    y_measured(:,iSubband) = I0plus(mp.Fend.corr.mask) - I0minus(mp.Fend.corr.mask);
end
%% Perform the estimation
ev = ekf_estimate(mp,ev,jacStruct,y_measured,closed_loop_command);


%% Save out the estimate
% TODO: add star and wavelength loop?
if mp.flagSim
    sbp_texp = mp.detector.tExpUnprobedVec; % exposure times for non-pairwise-probe images in each subband.
else
    sbp_texp  = mp.tb.info.sbp_texp;
end
ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);

for iSubband = 1:1:mp.Nsbp
    ev.Eest(:,iSubband) = (ev.x_hat(1:2:end,iSubband) + 1i*ev.x_hat(2:2:end,iSubband))/ (ev.e_scaling(iSubband) * sqrt(sbp_texp(iSubband)));
    if any(mp.dm_ind == 1);  ev.dm1.Vall(:, :, 1, iSubband) = mp.dm1.V;  end
    if any(mp.dm_ind == 2);  ev.dm2.Vall(:, :, 1, iSubband) = mp.dm2.V;  end
    if mp.unprobeImg
        if any(mp.dm_ind == 1)
            mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift + DM1Vdither + mp.dm1.V_shift);
        elseif any(mp.dm_ind_static == 1)
            mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz);
        end
        if any(mp.dm_ind == 2) 
            mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift + DM2Vdither + mp.dm2.V_shift);
        elseif any(mp.dm_ind_static == 2)
            mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz);
        end
        ev.imageArray(:,:,1,iSubband) = falco_get_sbp_image(mp, iSubband);
    else
        ev.imageArray(:, :, 1, iSubband) = ev.imageArrayMinus(:,:,1,iSubband);
    end
    ev.Im = ev.Im + mp.sbp_weights(iSubband) * ev.imageArray(:, :, 1, iSubband);
end
I0vec = y_measured./ev.peak_psf_counts;
ev.IincoEst = I0vec - abs(ev.Eest).^2; % incoherent light

%--Other data to save out
% TODO: not sure if this is returning the right thing? do I need to return
% ampNorm?
ev.ampSqMean = mean(I0vec(:)); %--Mean probe intensity
% ev.ampNorm = mean(I0vec(:)); %--Normalized probe amplitude maps

% ev.Im = ev.imageArray(:,:,1,mp.si_ref);
ev.IprobedMean = mean(mean(ev.imageArray)); 

mp.isProbing = false;

%% If itr = itr_OL get OL data. NOTE THIS SHOULD BE BEFORE THE 
% "Remove control from DM command so that controller images are correct" block


if any(mp.est.itr_ol==ev.Itr) == true
    [mp,ev] = get_open_loop_data(mp,ev);
else
    ev.IOLScoreHist(ev.Itr,:) = ev.IOLScoreHist(ev.Itr-1,:);
end

%% Remove control from DM command so that controller images are correct
if any(mp.dm_ind == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift + DM1Vdither + mp.dm1.V_shift);
elseif any(mp.dm_ind_static == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz);
end
if any(mp.dm_ind == 2) 
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift + DM2Vdither + mp.dm2.V_shift);
elseif any(mp.dm_ind_static == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz);
end

save_ekf_data(mp,ev,DM1Vdither, DM2Vdither)

fprintf(' done. Time: %.3f\n',toc);
end

function [ev] = ekf_estimate(mp, ev, jacStruct, y_measured, closed_loop_command)
%% Estimation part. All EKFs are avanced in parallel
if mp.flagSim
    sbp_texp = mp.detector.tExpUnprobedVec;
else
    sbp_texp = mp.tb.info.sbp_texp;
end
dy_hat = [];
H_T = [];
for iSubband = 1:1:mp.Nsbp
    %--Estimate of the closed loop electric field:
    G = (ev.G_tot_cont(:,:,iSubband)*ev.e_scaling(iSubband))*sqrt(sbp_texp(iSubband));
    u_hat_CL = ev.u_hat - closed_loop_command;
    %--Estimate of the measurement:
    dy_dx_hat_CL = 4 * G * ev.dither;
    dy_hat_real_imag = dy_dx_hat_CL + G * u_hat_CL;
    dy_hat = [dy_hat, dy_hat_real_imag(1:2:end) + dy_hat_real_imag(2:2:end)];
    H_T = [H_T, G(1:2:end, :)'*diag(dy_dx_hat_CL(1:2:end)) + G(2:2:end, :)'*diag(dy_dx_hat_CL(2:2:end))];
    H = H_T';
    ev.P(:,:,iSubband) = ev.P(:,:,iSubband) + ev.Q(:,:,iSubband);
    
    ev.R = eye(length(dy_hat))* (10^(mp.ctrl.sched_mat(ev.Itr, 2)));
   
    P_H_T = mypagemtimes(ev.P(:,:,:,iSubband), H');
    S = mypagemtimes(H, P_H_T) + ev.R;
    S_inv = mypageinv(S) ;%does this need to be a pinv?

    % S_inv = np.linalg.pinv(S)
    K = mypagemtimes(P_H_T, S_inv);
    ev.P(:,:,:,iSubband) = ev.P(:,:,:,iSubband) - mypagemtimes(P_H_T, permute(K,[2,1,3]));
    
    % EKF correction:
    ev.u_hat = K*(y_measured(:,iSubband) - dy_hat);
    
    ev.x_hat(:,iSubband) = G * ev.u_hat;

end

end


function comm_vector = get_dm_command_vector(mp,command1, command2)

if any(mp.dm_ind == 1); comm1 = command1(mp.dm1.act_ele);  else; comm1 = []; end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2); comm2 = command2(mp.dm2.act_ele);  else; comm2 = []; end
comm_vector = [comm1;comm2];

end


function ev = pinned_act_safety_check(mp,ev)
% Update new pinned actuators
if any(mp.dm_ind == 1) || any(mp.dm_drift_ind == 1)
    ev.dm1.new_pinned_actuators = setdiff(mp.dm1.pinned, ev.dm1.initial_pinned_actuators);
    ev.dm1.act_ele_pinned = mp.dm1.pinned(ismember(ev.dm1.new_pinned_actuators,mp.dm1.act_ele));
end
if any(mp.dm_ind == 2) || any(mp.dm_drift_ind == 2)
    ev.dm2.new_pinned_actuators = setdiff(mp.dm2.pinned, ev.dm2.initial_pinned_actuators);
    ev.dm2.act_ele_pinned = mp.dm2.pinned(ismember(ev.dm2.new_pinned_actuators,mp.dm2.act_ele));

end


%  Check that no new actuators have been pinned
if size(ev.dm1.new_pinned_actuators,2)>0 || size(ev.dm2.new_pinned_actuators,2)>0

    % Print error warning
    fprintf('New DM1 pinned: [%s]\n', join(string(ev.dm1.new_pinned_actuators), ','));
    fprintf('New DM2 pinned: [%s]\n', join(string(ev.dm2.new_pinned_actuators), ','));

    % If actuators are used in jacobian, quit.
    if size(ev.dm1.act_ele_pinned,2)>0 || size(ev.dm2.act_ele_pinned,2)>0
        save(fullfile([mp.path.config,'/','/ev_exit_',num2str(ev.Itr),'.mat']),'ev')
        save(fullfile([mp.path.config,'/','/mp_exit_',num2str(ev.Itr),'.mat']),"mp")

        error('New actuators in act_ele pinned, exiting loop');
    end
end

end


function out = mypageinv(in)

dim = size(in,3);
out = zeros(size(in));
for i = 1:dim
    out(:,:,i) = inv(in(:,:,i));
end

end

function out = mypagemtimes(X,Y) 

dim1 = size(X,3);
dim2 = size(Y,3);
if(dim1~=dim2); error('X and Y need to be the same size.'); end
out = zeros(size(X,1),size(Y,2),dim1);
for i = 1:dim1
    out(:,:,i) = X(:,:,i)*Y(:,:,i);
end

end

function [mp,ev] = get_open_loop_data(mp,ev)
%% Remove control and dither from DM command 

% If DM is used for drift and control, apply V_dz and Vdrift, if DM is only
% used for control, apply V_dz
if (any(mp.dm_drift_ind == 1) && any(mp.dm_ind == 1)) || any(mp.dm_drift_ind == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift);
elseif any(mp.dm_ind == 1) || any(mp.dm_ind_static == 1)
    mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz);
end

if (any(mp.dm_drift_ind == 2) && any(mp.dm_ind == 2)) || any(mp.dm_drift_ind == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift);
elseif any(mp.dm_ind == 2) || any(mp.dm_ind_static == 2)
    mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz);
end

% Do safety check for pinned actuators
disp('OL DM safety check.')
ev = pinned_act_safety_check(mp,ev);

if ev.Itr == 1
    ev.IOLScoreHist = zeros(mp.Nitr,mp.Nsbp);
end

I_OL = zeros(size(ev.imageArray(:,:,1,1),1),size(ev.imageArray(:,:,1,1),2),mp.Nsbp);
for iSubband = 1:mp.Nsbp
    I0 = falco_get_sbp_image(mp, iSubband);
    I_OL(:,:,iSubband) = I0;
    
    ev.IOLScoreHist(ev.Itr,iSubband) = mean(I0(mp.Fend.score.mask));
    
end


ev.normI_OL_sbp = I_OL;

disp(['mean OL contrast: ',num2str(mean(ev.IOLScoreHist(ev.Itr,:)))])
end

function save_ekf_data(mp,ev,DM1Vdither, DM2Vdither)
drift = zeros(mp.dm1.Nact,mp.dm1.Nact,length(mp.dm_drift_ind));
dither = zeros(mp.dm1.Nact,mp.dm1.Nact,length(mp.dm_ind));
efc = zeros(mp.dm1.Nact,mp.dm1.Nact,length(mp.dm_ind));


if mp.dm_drift_ind(1) == 1; drift(:,:,1) = mp.dm1.V_drift;end
if mp.dm_drift_ind(1) == 2; drift(:,:,1) = mp.dm2.V_drift ; else drift(:,:,2) = mp.dm2.V_drift; end


if mp.dm_ind(1) == 1; dither(:,:,1) = DM1Vdither;end
if mp.dm_ind(1) == 2; dither(:,:,1) = DM2Vdither ; else dither(:,:,2) = DM2Vdither; end

if mp.dm_ind(1) == 1; efc(:,:,1) = mp.dm1.dV;end
if mp.dm_ind(1) == 2; efc(:,:,1) = mp.dm2.dV ; else efc(:,:,2) = mp.dm2.dV; end

% TODO: move to plot_progress_iact
fitswrite(drift,fullfile([mp.path.config,'/','/drift_command_it',num2str(ev.Itr),'.fits']))
fitswrite(dither,fullfile([mp.path.config,'/','dither_command_it',num2str(ev.Itr),'.fits']))
fitswrite(efc,fullfile([mp.path.config,'/','efc_command_it',num2str(ev.Itr-1),'.fits']))

if ev.Itr == 1
    dz_init = zeros(mp.dm1.Nact,mp.dm1.Nact,length(mp.dm_ind));
    if mp.dm_ind(1) == 1; dz_init(:,:,1) = mp.dm1.V_dz;end
    if mp.dm_ind(1) == 2; dz_init(:,:,1) = mp.dm2.V_dz ; else dz_init(:,:,2) = mp.dm2.V_dz; end

    fitswrite(dz_init,fullfile([mp.path.config,'/','dark_zone_command_0_pwp.fits']))
end

end
