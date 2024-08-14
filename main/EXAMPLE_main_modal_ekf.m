%% Section 1: Initialization of Parameters
n_iter = 5;
num_pix_dark_zone = 10;  % reduced for simplicity in testing
num_actuators = 6;  % reduced for simplicity in testing
wavelengths = 640;

%% Section 2: Dummy Data Generation
direct_image_data = ones(num_pix_dark_zone, num_pix_dark_zone) * 1000;
jacobian_data = normrnd(0, 1, [num_actuators, 2 * num_pix_dark_zone]);
dark_zone_mask = int16(ones(num_pix_dark_zone, num_pix_dark_zone));  % Convert bool to int16

% Simulate the creation of dummy FITS files
fitswrite(direct_image_data, 'dummy_direct_image.fits');
fitswrite(jacobian_data, 'dummy_jacobian.fits');
fitswrite(dark_zone_mask, 'dummy_dark_zone_mask.fits');

%% Section 3: Load and Check Dummy FITS Files
direct_images = fitsread('dummy_direct_image.fits');
direct_maxes_binned = max(direct_images(:));  % counts/s
e_scaling = sqrt(direct_maxes_binned);  % sqrt(counts/s)

%% Section 4: Rearrange Jacobian
jacobian = fitsread('dummy_jacobian.fits');
G_reordered = zeros(size(jacobian));
G_reordered(:, 1:2:end) = jacobian(:, 1:size(jacobian, 2) / 2);
G_reordered(:, 2:2:end) = jacobian(:, size(jacobian, 2) / 2 + 1:end);

jacobians_boston_pixel = G_reordered;

%% Section 5: Function to Get DM Modes and Drift
function [M, QM] = get_DM_modes_and_drift(G, r)
    [U, S, V] = svd(G, 'econ');

    if isempty(r)
        r = length(diag(S));
    end

    M = V(:, 1:r);
    QM = diag(diag(S(1:r, 1:r))) * (U(:, 1:r)' * eye(size(G, 1)) * U(:, 1:r)) * diag(diag(S(1:r, 1:r)));
end

[M, Q] = get_DM_modes_and_drift(jacobians_boston_pixel, 3);

%% Section 6: Initialize Variables for Testing
r = size(M, 2);
x_hat = zeros(r, 1);
F = eye(r);
P = Q * 10;
E_floor = zeros(2 * num_pix_dark_zone, 1);
closed_loop_correction = zeros(num_actuators, 1);

%% Section 7: Function to Take Images
function images = take_images(cmd, num_pix_dark_zone, wavelengths)
    images = struct();
    extended_cmd = repmat(cmd, ceil(num_pix_dark_zone / length(cmd)), 1);
    extended_cmd = extended_cmd(1:num_pix_dark_zone);  % Resize to match image size

    for wl = wavelengths
        images.(sprintf('wl%d', wl)) = extended_cmd * 0.1 + normrnd(0, 0.01, [num_pix_dark_zone, 1]);
    end
end

x_hat_history = zeros(n_iter, r);

%% Section 8: Test Loop with Adjusted Kalman Filter
try
    for i = 1:n_iter
        probe = normrnd(0, 0.1, [num_actuators, 1]);

        closed_loop_command_plus = closed_loop_correction + probe;
        closed_loop_command_minus = closed_loop_correction - probe;

        images_plus = take_images(closed_loop_command_plus, num_pix_dark_zone, wavelengths);
        images_minus = take_images(closed_loop_command_minus, num_pix_dark_zone, wavelengths);

        y = images_plus.wl640 - images_minus.wl640;  % Result is a column vector

        Gu = jacobians_boston_pixel' * closed_loop_correction;  % Gu is now a column vector
        Gprobe = jacobians_boston_pixel' * probe;  % Gprobe is now a column vector

        E_hat = Gu + M * x_hat + E_floor;  % E_hat is a column vector

        y_hat = 4 * (Gprobe(1:2:end) .* E_hat(1:2:end) + Gprobe(2:2:end) .* E_hat(2:2:end));  % y_hat is a column vector
        
        R = diag(2 * (E_hat(1:2:end).^2 + E_hat(2:2:end).^2 + Gprobe(1:2:end).^2 + Gprobe(2:2:end).^2));

        % Reshape H to produce the same number of outputs as y_hat
        H_full = zeros(length(y_hat), r);  % Adjust H to match y_hat dimensions
        H_full(:, :) = 4 * (Gprobe(1:2:end) .* M(1:2:end, :) + Gprobe(2:2:end) .* M(2:2:end, :));

        % Check dimensions before S computation
        fprintf('Dimensions before S computation:\n');
        fprintf('H_full: [%d, %d], P: [%d, %d], R: [%d, %d]\n', size(H_full), size(P), size(R));
        
        S = H_full * P * H_full' + R;  % S should now match the dimensions of R

        % Check S dimensions
        fprintf('S: [%d, %d]\n', size(S));

        K = P * H_full' / S;  % K should match the state dimension

        % Check K dimensions
        fprintf('K: [%d, %d]\n', size(K));

        x_hat = x_hat + K * (y - y_hat);  % Ensure dimensions align

        x_hat_history(i, :) = x_hat';

        % Print intermediate values for comparison
        fprintf('Iteration %d\n', i);
        fprintf('x_hat: %s\n', mat2str(round(x_hat', 4)));
        fprintf('P: \n%s\n\n', mat2str(round(P, 4)));
        fprintf('y_hat: %s\n', mat2str(round(y_hat, 4)));
        fprintf('y: %s\n', mat2str(round(y, 4)));
        fprintf('Gu: %s\n', mat2str(round(Gu', 4)));
        fprintf('Gprobe: %s\n', mat2str(round(Gprobe', 4)));
        fprintf('E_hat: %s\n', mat2str(round(E_hat', 4)));
        fprintf('H_full: %s\n', mat2str(round(H_full, 4)));
        fprintf('S: %s\n', mat2str(round(S, 4)));
        fprintf('K: %s\n', mat2str(round(K, 4)));
        fprintf('--------------------------------------------------\n');
    end
catch ME
    fprintf('Error during iteration %d: %s\n', i, ME.message);
end
