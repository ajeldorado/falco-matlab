function [ev, mp] = falco_est_ekf_maintenance_modal(mp, ev, varargin)

end %--END OF FUNCTION


function ev = modal_ekf_estimate(mp,ev,jacStruct,y_measured,dither,cont_command)

    y_plus = y_measured(:,:,1);
    y_minus = y_measured(:,:,end);

for iSubband = 1:1:mp.Nsbp

     %--Jacobian (convert from sqrt(contrast)/nm to sqrt(photons/s)/nm)
    G = mp.G_tot_cont(:, :, iSubband) * mp.est.e_scaling(iSubband) * sqrt(ev.sbp_texp); % Jacobian for the full DM
    
    %--Dither-modulated Jacobian
    J_du = G' * dither;
    
    %--Measurement noise covariance matrix (R)
    % if mp.est.low_photon_regime
        % dark_current_photons = mp.est.dark_current * mp.est.quantum_efficiency * mp.exp.t_coron_sbp;
        % read_noise_photons = mp.est.read_noise * mp.est.quantum_efficiency;
        % ev.dm_R = diag(y_prev + y_measured + 2 * dark_current_photons^2 + 2 * read_noise_photons^2);
    % else
    dark_current_photons = mp.est.dark_current * mp.est.quantum_efficiency * ev.sbp_texp;
    read_noise_photons = mp.est.read_noise * mp.est.quantum_efficiency;

    ev.R = diag(y_plus + y_minus + 2 * dark_current_photons^2 + 2 * read_noise_photons^2);
    % end

    %--Measurement operator (H)
    JJdu = G * diag(J_du);
    ev.H = 4 * (JJdu(:, 1:2:end) + JJdu(:, 2:2:end));
    
    %--Prediction step for the covariance matrix
    ev.P(:, :, iSubband) = ev.P(:, :, iSubband) + ev.Q(:, :, iSubband);
    
    %--Kalman Gain Calculation
    P_H = ev.P(:, :, iSubband) * ev.H';
    S = ev.H * P_H + ev.R;
    
    if rcond(S) < 1e-12
        K = P_H * pinv(S);
        fprintf('Warning: S matrix is singular. Using pseudo-inverse.\n');
    else
        K = P_H / S;
    end
    
    %--Update step for the covariance matrix
    ev.P(:, :, iSubband) = ev.P(:, :, iSubband) - K * P_H';
    
    %--Measurement residual (dy)
    dy = y_plus - y_minus;
    
    %--Modelled intensity difference
    controlled_command = cont_command + get_dm_command_vector(mp,mp.dm1.V_dz, mp.dm2.V_dz);

    E_hat_plus = G' * (controlled_command + dither);
    E_hat_minus = G' * (controlled_command - dither);
    I_hat_plus = E_hat_plus(1:2:end).^2 + E_hat_plus(2:2:end).^2;
    I_hat_minus = E_hat_minus(1:2:end).^2 + E_hat_minus(2:2:end).^2;
    dy_hat = I_hat_plus - I_hat_minus;
    
    %--Update DM command estimate
    residual = dy - dy_hat;
    ev.x_hat(:, iSubband) = ev.x_hat(:, iSubband) + K * residual;
    
    ev.Eest(:,iSubband) = (E_hat_plus(1:2:end) + 1i*E_hat_plus(2:2:end));

end

% % Prep estimate to save out
% for iSubband = 1:1:mp.Nsbp
%     ev.Eest(:,iSubband) = (ev.x_hat(1:2:end,iSubband) + 1i*ev.x_hat(2:2:end,iSubband))/ (ev.e_scaling(iSubband) * sqrt(ev.sbp_texp(iSubband)));
% end
end