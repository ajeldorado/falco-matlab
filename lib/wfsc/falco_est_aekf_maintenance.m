function ev = falco_est_aekf_maintenance(mp, ev, varargin)
% Modified version of falco_est_ekf_maintenance to handle incoherent estimation
% Uses 3-element state vector per pixel: [real(E), imag(E), incoherent_intensity]

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
ev.imageArray = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);
ev.Eest = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IincoEst = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IprobedMean = 0;
ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
if whichDM == 1;  ev.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, 1, mp.Nsbp);  end
if whichDM == 2;  ev.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, 1, mp.Nsbp);  end

%% Get dither command
% Set random number generator seed
% Dither commands get re-used every dither_cycle_iters iterations
if mod(Itr-1, mp.est.dither_cycle_iters) == 0 || Itr == 1
    ev.dm1_seed_num = 0; 
    ev.dm2_seed_num = 1000; % Don't want same random commands on DM1 and DM2
    disp(['Dither random seed reset at iteration ', num2str(Itr)])
else
    ev.dm1_seed_num = ev.dm1_seed_num + 1; 
    ev.dm2_seed_num = ev.dm2_seed_num + 1;
end

% Generate random dither command
if any(mp.dm_ind == 1)  
    rng(ev.dm1_seed_num);
    dither_command1 = mp.est.dither.*randn(mp.dm1.Nact,1);%sqrt(iteration);
    dither_command1(mp.dm1.pinned) = 0; % Don't dither pinned actuators
else
    dither_command1 = zeros(mp.dm1.Nact,1);
end

if any(mp.dm_ind == 2)  
    rng(ev.dm2_seed_num);
    dither_command2 = mp.est.dither.*randn(mp.dm2.Nact,1);
    dither_command2(mp.dm2.pinned) = 0; % Don't dither pinned actuators
else
    dither_command2 = zeros(mp.dm2.Nact,1);
end

%% Injection of drift
%FALCO drift is applied to ev.dm1.V, while the estimator assumes drift command is dV.

% Get drift comand vector that estimator knows about
if ~isfield(ev,'drift_command_total')
   ev.drift_command_total = zeros(mp.dm1.Nact + mp.dm2.Nact,1);
end

drift_command = falco_drift_injection(mp, ev);
ev.drift_command_total = ev.drift_command_total + drift_command; 

%% Apply probing command
% Get the current probe command (this is probably a bad proxy function)

if mp.est.probe.whichDM == 1
    probe_command1 = dither_command1;
    probe_command2 = dither_command2;
elseif mp.est.probe.whichDM == 2
    probe_command1 = dither_command1;
    probe_command2 = dither_command2;
end

mp.dm1.V = mp.dm1.V + probe_command1;
mp.dm2.V = mp.dm2.V + probe_command2;

%% Take Image
if mp.flagFiber
    [ev.Im, ev.Ifiber] = falco_get_summed_image(mp);
else
    ev.Im = falco_get_summed_image(mp);
end

%% Perform the estimation
% Get the DM commands known to the estimator
closed_loop_command = get_dm_command_vector(mp,mp.dm1.V - mp.dm1.V_dz, mp.dm2.V - mp.dm2.V_dz);

if mp.flagFiber
    if mp.est.flagUseJac
        ev = ekf_estimate_incoherent(mp, ev, jacStruct, ev.Ifiber, closed_loop_command);
    else
        ev = ekf_estimate_incoherent(mp, ev, [], ev.Ifiber, closed_loop_command);
    end
else
    Im_meas = ev.Im(mp.Fend.corr.maskBool);
    if mp.est.flagUseJac
        ev = ekf_estimate_incoherent(mp, ev, jacStruct, Im_meas, closed_loop_command);
    else
        ev = ekf_estimate_incoherent(mp, ev, [], Im_meas, closed_loop_command);
    end
end

%% Pack up E-field estimates
% Get the closed-loop E-field estimate back to contrast units:
for iSubband = 1:mp.Nsbp
    ev.Eest(:,iSubband) = (ev.x_hat(1:ev.SS:end,iSubband) + 1i*ev.x_hat(2:ev.SS:end,iSubband)) ./ (ev.e_scaling(iSubband) * sqrt(ev.peak_psf_counts(iSubband)/ev.peak_psf_counts(iSubband)));
    ev.IincoEst(:,iSubband) = ev.x_hat(3:ev.SS:end,iSubband) ./ ev.peak_psf_counts(iSubband);
end

%% Remove the probing command
mp.dm1.V = mp.dm1.V - probe_command1;
mp.dm2.V = mp.dm2.V - probe_command2;

%% Safety check on pinned actuators (remove once this is tested)
ev = pinned_act_safety_check(mp,ev);

%% Optional output showing image
if mp.flagDisplayProgress
    if mp.flagFiber
        ev.Im = falco_get_summed_image(mp);
    end
    
    %--Plot the raw image in the dark hole
    figure(300);
    imagesc(mp.Fend.xisDL, mp.Fend.etasDL, log10(ev.Im), [-10, -4]); axis xy equal tight; 
    colorbar; colormap parula; drawnow;
    pause(0.1);
end

%% Control inputs
% Add the dithering command to the DM drift command as the
% estimator does not know about it

efc_command = get_dm_command_vector(mp,probe_command1,probe_command2);

% Estimate the closed loop E-field 
if mp.flagSim
    sbp_texp  = mp.detector.tExpUnprobedVec;
else
    sbp_texp  = mp.tb.info.sbp_texp;
end

for iSubband = 1:1:mp.Nsbp
    ev.x_hat(:,iSubband) = ev.x_hat(:,iSubband) + (ev.G_tot_cont(:,:,iSubband)*ev.e_scaling(iSubband))*sqrt(sbp_texp(iSubband))*efc_command;
end
mp.dm1.V_shift = mp.dm1.dV;
mp.dm2.V_shift = mp.dm2.dV;

mp.dm1.dV = zeros(size(mp.dm1.V_dz));
mp.dm2.dV = zeros(size(mp.dm2.V_dz));
    % TODO: save each drift command

%% =========================================================================
function [ev] = ekf_estimate_incoherent(mp, ev, jacStruct, y_measured, closed_loop_command)
%% Estimation part. All EKFs are advanced in parallel - MODIFIED FOR 3-ELEMENT STATE

if mp.flagSim
    sbp_texp = mp.detector.tExpUnprobedVec;
else
    sbp_texp = mp.tb.info.sbp_texp;
end

for iSubband = 1:1:mp.Nsbp

    %--Estimate of the closed loop electric field (3-element state):
    x_hat_CL = ev.x_hat(:,iSubband) + (ev.G_tot_cont(:,:,iSubband)*ev.e_scaling(iSubband))*sqrt(sbp_texp(iSubband))*closed_loop_command;

    %--Estimate of the measurement (modified for 3-element state):
    % Total intensity = coherent + incoherent + dark current
    y_hat = x_hat_CL(1:3:end).^2 + x_hat_CL(2:3:end).^2 + x_hat_CL(3:3:end) + (mp.est.dark_current*sbp_texp(iSubband));

    ev.R(ev.R_indices) = reshape(y_hat + (mp.est.read_noise)^2,size(ev.R(ev.R_indices)));

    %--CORRECTED observation matrix H for 3-element state:
    % Use same pattern as standard EKF but adapted for 3-element state
    
    % Debug information - ADD THIS FOR DEBUGGING
    fprintf('=== AEKF Debug Info ===\n');
    fprintf('ev.SS = %d, ev.BS = %d, ev.SL = %d\n', ev.SS, ev.BS, ev.SL);
    fprintf('x_hat_CL size: [%d, %d]\n', size(x_hat_CL));
    fprintf('mp.Fend.corr.Npix = %d\n', mp.Fend.corr.Npix);
    fprintf('y_measured size: [%d, %d]\n', size(y_measured));
    fprintf('Before H construction - ev.H size: [%d, %d, %d]\n', size(ev.H));
    
    % CORRECTED H matrix construction
    dh = zeros(size(ev.x_hat(:,iSubband)));
    dh(1:ev.SS:end) = 2 * x_hat_CL(1:3:end);  % Real components: 2*real(E)
    dh(2:ev.SS:end) = 2 * x_hat_CL(2:3:end);  % Imaginary components: 2*imag(E)  
    dh(3:ev.SS:end) = ones(mp.Fend.corr.Npix, 1);  % Incoherent components: 1

    ev.H(ev.H_indices) = dh;
    H_T = permute(ev.H,[2,1,3]);
    
    % More debug information
    fprintf('After H construction - H size: [%d, %d, %d]\n', size(ev.H));
    fprintf('H_T size: [%d, %d, %d]\n', size(H_T));
    fprintf('P size: [%d, %d, %d]\n', size(ev.P(:,:,:,iSubband)));
    fprintf('=== End AEKF Debug ===\n');

    ev.P(:,:,:,iSubband) = ev.P(:,:,:,iSubband) + ev.Q(:,:,:,iSubband);
   
    P_H_T = mypagemtimes(ev.P(:,:,:,iSubband), H_T);
    S = mypagemtimes(ev.H, P_H_T) + ev.R;
    S_inv = mypageinv(S) ;%does this need to be a pinv?

    % S_inv = np.linalg.pinv(S)
    K = mypagemtimes(P_H_T, S_inv);
    ev.P(:,:,:,iSubband) = ev.P(:,:,:,iSubband) - mypagemtimes(P_H_T, permute(K,[2,1,3]));
    
    % EKF correction:
    dy = (y_measured(:,iSubband) - y_hat);
    
    % Modified for 3-element state
    dy_hat_stacked  = zeros(size(K));
    dy_hat_stacked(1,:,:) = dy.';
    dy_hat_stacked(2,:,:) = dy.';
    dy_hat_stacked(3,:,:) = dy.'; % Third element for incoherent
    
    dx_hat_stacked = K.*dy_hat_stacked;

    dx_hat = zeros(size(x_hat_CL));
    dx_hat(1:3:end) = dx_hat_stacked(1,:,:); % Real part
    dx_hat(2:3:end) = dx_hat_stacked(2,:,:); % Imaginary part  
    dx_hat(3:3:end) = dx_hat_stacked(3,:,:); % Incoherent intensity

    ev.x_hat(:,iSubband) = ev.x_hat(:,iSubband) + dx_hat;

end

%% =========================================================================
function comm_vector = get_dm_command_vector(mp,command1, command2)

if any(mp.dm_ind == 1); comm1 = command1(mp.dm1.act_ele);  else; comm1 = []; end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2); comm2 = command2(mp.dm2.act_ele);  else; comm2 = []; end
comm_vector = [comm1;comm2];

%% =========================================================================
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
    if ~isempty(ev.dm1.new_pinned_actuators); fprintf('New DM1 pinned: [%s]\n', join(string(ev.dm1.new_pinned_actuators), ',')); end
    if ~isempty(ev.dm2.new_pinned_actuators); fprintf('New DM2 pinned: [%s]\n', join(string(ev.dm2.new_pinned_actuators), ',')); end

    % If actuators are used in jacobian, quit.
    if size(ev.dm1.act_ele_pinned,2)>0 || size(ev.dm2.act_ele_pinned,2)>0
        save(fullfile(mp.path.config,['ev_exit_',num2str(ev.Itr),'.mat']),'ev')
        save(fullfile(mp.path.config,['mp_exit_',num2str(ev.Itr),'.mat']),'mp')

        error('New actuators in act_ele pinned, exiting loop');
    end
end

%% =========================================================================
function out = mypageinv(in)

dim = size(in,3);
out = zeros(size(in));
for i = 1:dim
    out(:,:,i) = inv(in(:,:,i));
end

%% =========================================================================
function out = mypagemtimes(X,Y) 

dim1 = size(X,3);
dim2 = size(Y,3);
if(dim1~=dim2); 
    % Enhanced error message for debugging
    fprintf('ERROR in mypagemtimes: Matrix dimension mismatch!\n');
    fprintf('  X dimensions: [%d, %d, %d]\n', size(X));
    fprintf('  Y dimensions: [%d, %d, %d]\n', size(Y));
    error('X and Y need to have the same third dimension size. X has %d pages, Y has %d pages.', dim1, dim2); 
end

dim = size(X,3);
out = zeros(size(X,1),size(Y,2),dim);
for i = 1:dim
    out(:,:,i) = X(:,:,i)*Y(:,:,i);
end