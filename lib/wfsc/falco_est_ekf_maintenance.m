function ev = falco_est_ekf_maintenance(mp, ev, varargin)

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

% Initialize EKF
if Itr == 1
    
    [mp, ev, jacStruct] = initialize_ekf_maintenance(mp, ev, jacStruct);
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
% sfr: what is this doing
% if any(mp.dm_ind == 1);  DM1Vnom = mp.dm1.V;  else; DM1Vnom = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
% if any(mp.dm_ind == 2);  DM2Vnom = mp.dm2.V;  else; DM2Vnom = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1

% Initialize output arrays
ev.imageArray = zeros(mp.FendNeta, mp.Fend.Nxi, 1, mp.Nsbp);
ev.Eest = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IincoEst = zeros(mp.Fend.corr.Npix, mp.Nsbp*mp.compact.star.count);
ev.IprobedMean = 0;
ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
% TODO: do I need to update Vall?
if whichDM == 1;  ev.dm1.Vall = zeros(mp.dm1.Nact, mp.dm1.Nact, 1, mp.Nsbp);  end
if whichDM == 2;  ev.dm2.Vall = zeros(mp.dm2.Nact, mp.dm2.Nact, 1, mp.Nsbp);  end


%% Get Drift Command
mp = falco_drift_injection(mp);

%% Get dither command
% Is this if else necessary?

if any(mp.dm_ind == 1);  DM1Vdither = normrnd(0,ev.dither,[mp.dm1.Nact mp.dm1.Nact]); else; DM1Vdither = zeros(size(mp.dm1.V)); end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2);  DM2Vdither = normrnd(0,ev.dither,[mp.dm1.Nact mp.dm1.Nact]); else; DM2Vdither = zeros(size(mp.dm2.V)); end % The 'else' block would mean we're only using DM1

%% Set total command for estimator image
% TODO: need to save these commands for each iteration separately
% TODO: what is controller command? - mp.dm1.dV (make sure sign is right)

mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift + DM1Vdither + mp.dm1.dV);
mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift + DM2Vdither + mp.dm2.dV);

dither = get_dm_command_vector(mp,DM1Vdither, DM2Vdither);
if Itr > 1
    efc_command = get_dm_command_vector(mp,mp.dm1.dV, mp.dm2.dV);
else
    efc_command = 0*dither;
    mp.dm1.dV = zeros(size(DM1Vdither));
    mp.dm2.dV = zeros(size(DM1Vdither));
end
% TODO: check sign on efc command
closed_loop_command = dither + efc_command;

%% Get images

y_measured = zeros(length(mp.Fend.score.mask),mp.Nsbp);
for iSubband = 1:mp.Nsbp
    ev.imageArray(:,:,1,iSubband) = falco_get_sbp_image(mp, iSubband);
    I0 = ev.imageArray(:,:,1,iSubband) * ev.peak_psf_counts(iSubband);
    y_measured(:,iSubband) = I0(mp.Fend.score.mask);

end

%% Perform the estimation
ev = ekf_estimate(mp,ev,y_measured,closed_loop_command);


%% Save out the estimate
% TODO: add star and wavelength loop?
for si = 1:1:mp.Nsbp
    ev.Eest(:,si) = ev.x_hat(:,si) / (ev.e_scaling(si) * sqrt(mp.tb.info.sbp_texp(si)));
    if any(mp.dm_ind == 1);  ev.dm1.Vall(:, :, 1, si) = mp.dm1.V;  end
    if any(mp.dm_ind == 2);  ev.dm2.Vall(:, :, 1, si) = mp.dm2.V;  end
end
I0vec = y_measured.*ev.peak_psf_counts;
ev.IincoEst = I0vec - abs(ev.Eest).^2; % incoherent light

%--Other data to save out
% TODO: not sure if this is returning the right thing? do I need to return
% ampNorm?
ev.ampSqMean = mean(I0vec(:)); %--Mean probe intensity
% ev.ampNorm = mean(I0vec(:)); %--Normalized probe amplitude maps

ev.Im = ev.imageArray(:,:,1,mp.si_ref);
ev.IprobedMean = mean(mean(ev.imageArray)); 

mp.isProbing = false;

fprintf(' done. Time: %.3f\n',toc);

%% Remove control from DM command so that controller images are correct
mp.dm1 = falco_set_constrained_voltage(mp.dm1, mp.dm1.V_dz + mp.dm1.V_drift + DM1Vdither);
mp.dm2 = falco_set_constrained_voltage(mp.dm2, mp.dm2.V_dz + mp.dm2.V_drift + DM2Vdither);


end




function [ev] = ekf_estimate(mp, ev, y_measured, closed_loop_command)


% if ev.Itr == 0
%     closed_loop_command = dither + periodic_dm_shift;
% else
%     closed_loop_command = -controller.corrections(:,end) + dither + periodic_dm_shift;
% end

%% Estimation part. All EKFs are avanced in parallel
% e_scaling = zeros(mp.Nsbp,1);


for si = 1:1:mp.Nsbp

    %--Estimate of the closed loop electric field:
    x_hat_CL = ev.x_hat(:,si) + (jacStruct.G_tot(:,:,si).'*ev.e_scaling(si))*sqrt(ev.exposure_time_coron)*closed_loop_command;

    %--Estimate of the measurement:
    y_hat = x_hat_CL(1:ev.SS:end).^2 + x_hat_CL(2:ev.SS:end).^2 + (mp.dark_current*ev.exposure_time_coron);

    ev.R(ev.R_indices) = reshape(y_hat + (ev.read_noise)^2,size(ev.R(ev.R_indices)));

    ev.H(ev.H_indices) = 2 * x_hat_CL;
%     H_T = H.transpose(0,2,1);
    % TODO: check that this transpose is correct
    H_T = permute(ev.H,[2,1,3]);

    ev.P(:,:,:,si) = ev.P(:,:,:,si) + ev.Q(:,:,:,si);
   
    P_H_T = pagemtimes(ev.P(:,:,:,si), H_T);
    S = pagemtimes(ev.H, P_H_T) + ev.R;
    S_inv = pageinv(S) ;%does this need to be a pinv?

    % S_inv = np.linalg.pinv(S)
    K = pagemtimes(P_H_T, S_inv);
    ev.P(:,:,:,si) = ev.P(:,:,:,si) - pagemtimes(P_H_T, permute(K,[2,1,3]));
    
%     EKF correction:
%     dx_hat = pagemtimes(K, (y_measured(:,si) - y_hat).reshape((-1,ev.BS/ev.SS,1))).reshape(-1);
%     dx_hat = permute(K,[1,3,2])*(y_measured(:,si) - y_hat);%.reshape((-1,ev.BS/ev.SS,1))).reshape(-1);
    dy = (y_measured(:,si) - y_hat);
    
    dy_hat_stacked  = zeros(size(K));
    dy_hat_stacked(1,:,:) = dy.';
    dy_hat_stacked(2,:,:) = dy.';
    
    dx_hat_stacked = K.*dy_hat_stacked;

    dx_hat = zeros(size(x_hat_CL));
    dx_hat(1:ev.SS:end) = dx_hat_stacked(1,:,:);
    dx_hat(2:ev.SS:end) = dx_hat_stacked(2,:,:);


    ev.x_hat(:,si) = ev.x_hat(:,si) + dx_hat;
%     
%     % Convert E_hat to contrast units:
%     E_hat = (ev.x_hat(:,si)(::SS) + complex(0,1)*ev.x_hat(wavelength)(1::SS))/(ev.e_scaling(wavelength) * np.sqrt(ev.exposure_time_coron)); %Estimate of the electric field from EKF state estimate
%              
%     
%     % Update based on scaling factors
% %             x_hat(wavelength) *= e_field_scale_factors(wavelength)(-1)
% %             E_hat *= e_field_scale_factors(wavelength)(-1)
%     
%     E = hcipy.Field(np.zeros(wfsc_utils.focal_grid.size, dtype='complex'), wfsc_utils.focal_grid);
%     E(dark_zone) = E_hat;
%     E_estimates(wavelength) = E;
%     

end



end


function comm_vector = get_dm_command_vector(mp,command1, command2)

if any(mp.dm_ind == 1); comm1 = command1(mp.dm1.act_ele);  else; comm1 = []; end % The 'else' block would mean we're only using DM2
if any(mp.dm_ind == 2); comm2 = command2(mp.dm2.act_ele);  else; comm2 = []; end
comm_vector = [comm1;comm2];

end













