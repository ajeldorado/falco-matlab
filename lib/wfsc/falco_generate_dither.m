function [mp,ev, DM1Vdither, DM2Vdither] = falco_generate_dither(mp,ev)

%% Get dither command
% Set random number generator seed
% Dither commands get re-used every dither_cycle_iters iterations
Itr = ev.Itr;
DM1Vdither = zeros(size(mp.dm1.V));
DM2Vdither = zeros(size(mp.dm2.V));
if mod(Itr-1, mp.est.dither_cycle_iters) == 0 || Itr == 1
    ev.dm1_seed_num = 0; 
    ev.dm2_seed_num = 1000; % Don't want same random commands on DM1 and DM2
    disp(['Dither random seed reset at iteration ', num2str(Itr)])
else
    ev.dm1_seed_num = ev.dm1_seed_num + 1; 
    ev.dm2_seed_num = ev.dm2_seed_num + 1;
end

if isfield(mp.est, 'dither_type')
    dither_type = mp.est.dither_type;
else
    dither_type = 'random';
end

switch lower(dither_type)
    case{'random'}
        if any(mp.dm_ind_est == 1)  
            rng(ev.dm1_seed_num); 
            DM1Vdither(mp.dm1.act_ele) = normrnd(0,mp.est.dither,[mp.dm1.Nele, 1]);  
        end % The 'else' block would mean we're only using DM2
        
        if any(mp.dm_ind_est == 2)  
            rng(ev.dm2_seed_num); 
            DM2Vdither(mp.dm2.act_ele) = normrnd(0,mp.est.dither,[mp.dm2.Nele, 1]);  
        end % The 'else' block would mean we're only using DM1

    case{'generated'}
        %% Trying a new optimal dither test
        if any(mp.dm_ind_est == 1)  
            rng(ev.dm1_seed_num); 
            mask_opt_dm1 = mp.est.dither_opt(1:mp.dm1.Nele);
            % 1. Generate Standard Normal Noise (Mean=0, STD=1)
            random_noise_dm1 = randn([mp.dm1.Nele, 1]);
            
            % 2. Apply the Spatial Importance Map (Hadamard Mask)
            masked_dither_dm1 = random_noise_dm1 .* abs(mask_opt_dm1);
            
            % 3. Hardware Safety Clamp (Force STD to exactly mp.est.dither)
            actual_std_dm1 = std(masked_dither_dm1);
            scaling_factor_dm1 = mp.est.dither / max(actual_std_dm1, 1e-8);
            final_dither_dm1 = masked_dither_dm1 * scaling_factor_dm1;
            
            DM1Vdither(mp.dm1.act_ele) = final_dither_dm1; 
        end
        
        % === GENERATE MDZM DM2 DITHER ===
        if any(mp.dm_ind_est == 2)  
            rng(ev.dm2_seed_num); 

            mask_opt_dm2 = mp.est.dither_opt(mp.dm1.Nele+1:end);
            % 1. Generate Standard Normal Noise (Mean=0, STD=1)
            random_noise_dm2 = randn([mp.dm2.Nele, 1]);
            
            % 2. Apply the Spatial Importance Map (Hadamard Mask)
            masked_dither_dm2 = random_noise_dm2 .* abs(mask_opt_dm2);
            
            % 3. Hardware Safety Clamp (Force STD to exactly mp.est.dither)
            actual_std_dm2 = std(masked_dither_dm2);
            scaling_factor_dm2 = mp.est.dither / max(actual_std_dm2, 1e-8);
            final_dither_dm2 = masked_dither_dm2 * scaling_factor_dm2;
            DM2Vdither(mp.dm2.act_ele) = final_dither_dm2; 
        end

        case{'tophat'}
            nact = mp.dm1.Nact;
            act_x = linspace(-nact/2, nact/2 - 1, nact);
            act_y = linspace(-nact/2, nact/2 - 1, nact);
            [XA, YA] = meshgrid(act_x, act_y);
            R = sqrt(XA.^2 + YA.^2);
            cycle_order = reshape([1:(mp.est.dither_cycle_iters/2); (mp.est.dither_cycle_iters/2):(mp.est.dither_cycle_iters-1)], 1, []);
            i = cycle_order(mod(Itr - 1, length(cycle_order)) + 1);
            tp_vector = linspace(3, 24, mp.est.dither_cycle_iters);
            if any(mp.dm_ind_est == 1)
                temp2 = double(R < tp_vector(i));
                dither_dm1 = (1 - temp2) * mp.est.dither;
                DM1Vdither(mp.dm1.act_ele) = dither_dm1(mp.dm1.act_ele);        
            end
            
            if any(mp.dm_ind_est == 2)
                temp2 = double(R < tp_vector(i));
                dither_dm2 = (1 - temp2) * mp.est.dither;
                DM2Vdither(mp.dm2.act_ele) = dither_dm2(mp.dm2.act_ele);
            end

        case{'disk'}
            nact = mp.dm1.Nact;
            act_x = linspace(-nact/2, nact/2 - 1, nact);
            act_y = linspace(-nact/2, nact/2 - 1, nact);
            [XA, YA] = meshgrid(act_x, act_y);
            tp_vector = linspace(2, 21, mp.est.dither_cycle_iters);
            cycle_order = reshape([1:(mp.est.dither_cycle_iters/2); (mp.est.dither_cycle_iters/2):(mp.est.dither_cycle_iters-1)], 1, []);
            i = cycle_order(mod(Itr - 1, length(cycle_order)) + 1);
            diskstart = tp_vector(i);
            
            temp1 = double(sqrt(XA.^2 + YA.^2) > diskstart);
            temp2 = double(sqrt(XA.^2 + YA.^2) < diskstart + mp.est.dither_diskthick);
            if any(mp.dm_ind_est == 1)
                dither_dm1 = temp1.*temp2.*mp.est.dither;
                DM1Vdither(mp.dm1.act_ele) = dither_dm1(mp.dm1.act_ele);        
            end
            
            if any(mp.dm_ind_est == 2)
                dither_dm2 = temp1.*temp2.*mp.est.dither;
                DM2Vdither(mp.dm2.act_ele) = dither_dm2(mp.dm2.act_ele);
            end
         
end 

end




