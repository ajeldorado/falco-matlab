% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to train an adaptive model to improve the Jacobian based on 
% measured data and the model.
%
% INPUTS:
% -mp = structure of model parameters
% -ev = structure containing variables for the estimator
%
% OUTPUTS:
% -mp = structure of model parameters
%
% REVISION HISTORY:
% --------------
% Moved on 2019-09-30 by A.J. Riggs to be a separate file. Had tried
%   nesting in falco_wfsc_loop to save RAM but that did not work.

function mp = falco_train_model(mp,ev)
    Itr = ev.Itr;
    n_batch = mp.NitrTrain;
    % n_batch = 5; %% INITIALIZE THE DATA STRUCTURE THAT SAVES THE TRAINING DATA

    if(Itr==1)
        data_train.u1 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u2 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u1p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.u2p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.I = zeros(size(mp.Fend.corr.mask, 1), size(mp.Fend.corr.mask, 2), 2*mp.est.probe.Npairs+1, n_batch);
    else
        data_train = mp.data_train;
    end

    data_train.u1(:, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = mp.dm1.dV;
    data_train.u2(:, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = mp.dm2.dV;
    data_train.u1p(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = ev.Vcube.dm1;
    data_train.u2p(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = ev.Vcube.dm2;
    data_train.I(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = ev.Icube;

    if rem(Itr, n_batch) == 0
        % convert the WFSC data to standard input to the system ID function
        n_activeAct1 = length(mp.dm1.act_ele);
        n_activeAct2 = length(mp.dm2.act_ele);
        n_pairs = mp.est.probe.Npairs;
        n_pix = sum(sum(mp.Fend.corr.mask));

        uAll = zeros(n_activeAct1+n_activeAct2, n_batch); % control commands of all the iterations, including both DM1 and DM2
        uProbeAll = zeros(n_activeAct1+n_activeAct2, 2*n_pairs, n_batch); % probe commands of all the iterations, including both DM1 and DM2
        IAll = zeros(n_pix, 2*n_pairs+1, n_batch); % difference image of all the iterations


        for kc = 1 : n_batch % convert the 2D images and DM commands to vectors
            u1_2D = data_train.u1(:, :, kc);
            u2_2D = data_train.u2(:, :, kc);
            uAll(1:n_activeAct1, kc) = u1_2D(mp.dm1.act_ele);
            uAll(n_activeAct1+1:end, kc) = u2_2D(mp.dm2.act_ele);
            I_2D = data_train.I(:, :,1, kc);
            IAll(:, 1, kc) = I_2D(mp.Fend.corr.mask);
            for kp = 1 : 2*mp.est.probe.Npairs
                u1p_2D = data_train.u1p(:, :, kp+1, kc) - data_train.u1p(:, :, 1, kc);
                u2p_2D = data_train.u2p(:, :, kp+1, kc) - data_train.u2p(:, :, 1, kc);
                I_2D = data_train.I(:, :,kp+1, kc);

                uProbeAll(1:n_activeAct1, kp, kc) = u1p_2D(mp.dm1.act_ele);
                uProbeAll(n_activeAct1+1:end, kp, kc) = u2p_2D(mp.dm2.act_ele);
                IAll(:, kp+1, kc) = I_2D(mp.Fend.corr.mask);
            end
        end

        data_train.u1 = uAll(1:n_activeAct1, :);
        data_train.u2 = uAll(n_activeAct1+1:end, :);
        data_train.u1p = uProbeAll(1:n_activeAct1, :, :);
        data_train.u2p = uProbeAll(n_activeAct1+1:end, :, :);
        data_train.I = IAll;

        save([mp.path.jac, 'data_train.mat'],'data_train') %    save data_train data_train
        save([mp.path.jac, 'jacStruct.mat'], 'jacStruct'); %save jacStruct jacStruct

        if Itr == n_batch %--Call System ID after final iteration of training
            py.falco_systemID.linear_vl() %--First training
        else %--All later trainings
            Q0 = exp(jacStructLearned.noise_coef(1));
            Q1 = exp(jacStructLearned.noise_coef(2));
            R0 = exp(jacStructLearned.noise_coef(3));
            R1 = exp(jacStructLearned.noise_coef(4));
            R2 = exp(jacStructLearned.noise_coef(5));
            print_flag = false;
            path2data = mp.path.jac;
            lr = mp.est.lr;
            lr2 = mp.est.lr2;
            epoch = mp.est.epoch;
            py.falco_systemID.linear_vl(Q0, Q1, R0, R1, R2, lr, lr2, epoch, print_flag,path2data);
        end
        mp.flagUseLearnedJac = 1;
        data_train.u1 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u2 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u1p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.u2p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.I = zeros(size(mp.Fend.corr.mask, 1), size(mp.Fend.corr.mask, 2), 2*mp.est.probe.Npairs+1, n_batch);

    end 
    
    mp.data_train = data_train;

end %--END OF FUNCTION
