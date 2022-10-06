% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Run a wavefront estimation and control loop.

function [mp, out] = falco_wfsc_loop(mp, out)

fprintf('\nBeginning Trial %d of Series %d.\n', mp.TrialNum, mp.SeriesNum);
mp.thput_vec = zeros(mp.Nitr+1, 1);

flagBreak = false;

for Itr = 1:mp.Nitr
    
    
    %% Bookkeeping
    fprintf(['WFSC Iteration: ' num2str(Itr) '/' num2str(mp.Nitr) ', ' datestr(now) '\n' ]);
    
    % user-defined bookkeeping updates for each iteration
    if isfield(mp, 'funTopofloopBookkeeping') && ~isempty(mp.funTopofloopBookkeeping)
        mp = mp.funTopofloopBookkeeping(mp, Itr);
    end
    
    if mp.flagSim
        fprintf('Zernike modes used in this Jacobian:\t');
        fprintf('%d ', mp.jac.zerns);
        fprintf('\n');
    end
    
    ev.Itr = Itr;
    cvar.Itr = Itr;
    out.Itr = Itr;
    
    % Updated DM info
    if strcmpi(mp.controller, 'plannedefc')
        mp.dm_ind = mp.dm_ind_sched{Itr}; % Change which DMs are used
    end
    disp(['DMs to be used in this iteration = [' num2str(mp.dm_ind(:)') ']']);
    
    out.serialDateVec(Itr) = now;
    out.datetimeArray = datetime(out.serialDateVec,'ConvertFrom','datenum');
    out = store_dm_command_history(mp, out, Itr);

    %% Normalization and throughput calculations
    
    mp = falco_compute_psf_norm_factor(mp);
    
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
    out.thput(Itr) = max(thput);
    mp.thput_vec(Itr) = max(thput); % note: max() needed when mp.flagFiber==true
    
    %% Control Jacobian

    mp = falco_set_jacobian_modal_weights(mp); 
    
    if any(mp.relinItrVec == Itr)
        cvar.flagRelin = true;        
    else
        cvar.flagRelin = false;
    end
    
    if ((Itr == 1) && ~cvar.flagRelin && (strcmpi(mp.estimator, 'scc') || strcmpi(mp.estimator, 'iefc'))) % load jacStruct from file
        load([mp.path.jac filesep mp.jac.fn], 'jacStruct');
    elseif cvar.flagRelin % recompute jacStruct
        out.ctrl.relinHist(Itr) = true;
        jacStruct =  model_Jacobian(mp);
    end
    
    [mp, jacStruct] = falco_ctrl_cull_weak_actuators(mp, cvar, jacStruct);

    % Load the improved Jacobian if using the E-M technique
    if mp.flagUseLearnedJac
        jacStructLearned = load('jacStructLearned.mat');
        if any(mp.dm_ind == 1);  jacStruct.G1 = jacStructLearned.G1;  end
        if any(mp.dm_ind == 1);  jacStruct.G2 = jacStructLearned.G2;  end
    end
    
    %% Inject drift for (Only) Dark Zone Maintenance
    % Get Drift Command
    if strcmpi(mp.estimator,'ekf_maintenance')
        [mp, ev] = falco_drift_injection(mp, ev);
    end

    %% Wavefront Estimation
    
    if (Itr > 1); EestPrev = ev.Eest; end % save previous estimate for Delta E plot
    ev = falco_est(mp, ev, jacStruct);
    
    out = falco_store_intensities(mp, out, ev, Itr);

    %% Plot the expected and measured delta E-fields
    
    if (Itr > 1); EsimPrev = Esim; end % save previous value for Delta E plot
    Esim = compute_simulated_efield_for_delta_efield_plot(mp);
    if Itr > 1
        out = falco_plot_DeltaE(mp, out, ev.Eest, EestPrev, Esim, EsimPrev, Itr);
    end
    % Add model E-field to ev for saving
    ev.Esim = Esim;

    %% Progress plots (PSF, NI, and DM surfaces)
    % plot_wfsc_progress also saves images and ev probe data
    
    if Itr == 1; hProgress.master = 1; end % initialize the handle
    [out, hProgress] = plot_wfsc_progress(mp, out, ev, hProgress, Itr, ImSimOffaxis);
        
    %% Compute and Plot the Singular Mode Spectrum of the Electric Field
    if mp.flagSVD
        out = falco_plot_singular_mode_spectrum_of_Efield(mp, out, jacStruct, ev.Eest, Itr);
    end
        
    %% Wavefront Control
    
    % control strategy
    if isfield(mp, 'funCtrlStrategy') && ~isempty(mp.funCtrlStrategy)
        mp = mp.funCtrlStrategy(mp, out, Itr);
    end
    
    cvar.Eest = ev.Eest;
    [mp, cvar] = falco_ctrl(mp, cvar, jacStruct);
    
    % Save data to 'out'
    out = falco_store_controller_data(mp, out, cvar, Itr);
        
    %--Enforce constraints on DM commands 
    if any(mp.dm_ind == 1); mp.dm1 = falco_enforce_dm_constraints(mp.dm1); end
    if any(mp.dm_ind == 2); mp.dm2 = falco_enforce_dm_constraints(mp.dm2); end

    % Update the dynamic map of pinned actuators and tied actuators.
    % Actuators can be tied electrically (have same voltage) or by the
    % neighbor rule (have a constant offset).
    if any(mp.dm_ind == 1); mp.dm1 = falco_update_dm_constraints(mp.dm1); end
    if any(mp.dm_ind == 2); mp.dm2 = falco_update_dm_constraints(mp.dm2); end
    
    %--Update DM actuator gains for new voltages (stays same if 'fitType' == linear)
    if any(mp.dm_ind == 1); mp.dm1 = falco_update_dm_gain_map(mp.dm1); end
    if any(mp.dm_ind == 2); mp.dm2 = falco_update_dm_gain_map(mp.dm2); end
    
    %% Report Various Stats

    out = falco_compute_dm_stats(mp, out, Itr);

    %--Calculate sensitivities to small Zernike phase aberrations at entrance pupil.
    if ~isempty(mp.eval.Rsens) && ~isempty(mp.eval.indsZnoll)
        out.Zsens(:, :, Itr) = falco_get_Zernike_sensitivities(mp);
    end

    %--REPORTING NORMALIZED INTENSITY
    if out.InormHist(Itr+1) ~= 0
        if mp.flagFiber
            fprintf('Prev and New Measured NI (SMF):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
                out.InormFiberHist(Itr), out.InormFiberHist(Itr+1), out.InormFiberHist(Itr)/out.InormFiberHist(Itr+1) );
            fprintf('Prev and New Measured NI (pixels):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
                out.InormHist(Itr), out.InormHist(Itr+1), out.InormHist(Itr)/out.InormHist(Itr+1) );
        else
            fprintf('Prev and New Measured NI:\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
                out.InormHist(Itr), out.InormHist(Itr+1), out.InormHist(Itr)/out.InormHist(Itr+1) );
        end
        if ~mp.flagSim
            fprintf('\n');
        end
    else
        fprintf('Previous Measured NI:\t\t\t %.2e \n', out.InormHist(Itr))
    end

    %% Save 'out' structure after each iteration in case the trial terminates early.
    fnSnippet = [mp.path.config filesep mp.runLabel,'_snippet.mat'];
    fprintf('Saving data snippet to \n%s\n', fnSnippet)
    save(fnSnippet, 'out');
    fprintf('...done.\n\n')

    %% SAVE THE TRAINING DATA OR RUN THE E-M Algorithm
    if mp.flagTrainModel; mp = falco_train_model(mp,ev); end
    
    %% End early? You can change the value of bEndEarly in debugger mode, but you cannot change mp.Nitr or Itr
    if flagBreak
        break;
    end

end %--END OF ESTIMATION + CONTROL LOOP

Itr = mp.Nitr;

%% Update 'out' structure and progress plot one last time
Itr = Itr + 1;

out = store_dm_command_history(mp, out, Itr);

[mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
out.thput(Itr) = thput;
mp.thput_vec(Itr) = max(thput); % max() used for if mp.flagFiber==true

% Update progress plot using image from controller (if new image was taken)
if isfield(cvar, 'Im') && ~mp.ctrl.flagUseModel
    ev.Im = cvar.Im;
    [out, hProgress] = plot_wfsc_progress(mp, out, ev, hProgress, Itr, ImSimOffaxis);
end

if strcmpi(mp.estimator,'ekf_maintenance')  % sfr
   out.IOLScoreHist = ev.IOLScoreHist;
end


%% Save out an abridged workspace

fnSnippet = [mp.path.config filesep mp.runLabel,'_snippet.mat'];
fprintf('\nSaving data snippet to:\n\t%s\n', fnSnippet)
save(fnSnippet, 'out');
fprintf('...done.\n')

%% Save out the data from the workspace
if mp.flagSaveWS
    clear ev cvar G* h* jacStruct; % Save a ton of space when storing the workspace

    fnAll = fullfile(mp.path.ws,[mp.runLabel, '_all.mat']);
    disp(['Saving entire workspace to file ' fnAll '...'])
    save(fnAll);
    fprintf('done.\n\n')
else
    disp('Entire workspace NOT saved because mp.flagSaveWS==false')
end

end %--END OF main FUNCTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = store_dm_command_history(mp, out, Itr)

    if isfield(mp, 'dm1'); if(isfield(mp.dm1,'V'));  out.dm1.Vall(:,:,Itr) = mp.dm1.V;  end;  end
    if isfield(mp, 'dm2'); if(isfield(mp.dm2,'V'));  out.dm2.Vall(:,:,Itr) = mp.dm2.V;  end;  end
    if isfield(mp, 'dm8'); if(isfield(mp.dm8,'V'));  out.dm8.Vall(:,Itr) = mp.dm8.V(:);  end;  end
    if isfield(mp, 'dm9'); if(isfield(mp.dm9,'V'));  out.dm9.Vall(:,Itr) = mp.dm9.V(:);  end;  end

end


function out = falco_store_controller_data(mp, out, cvar, Itr)
    
    if any(strcmp({'plannedefc', 'gridsearchefc'}, mp.controller))
        out.ctrl.dmfacHist(Itr) = cvar.dmfacUsed;
        out.ctrl.log10regHist(Itr) = cvar.log10regUsed;
        out.log10regHist(Itr) = out.ctrl.log10regHist(Itr); % kept for backwards compatibility
    end

    % If the unprobed image for the next WFSC iteration was already taken,
    % then use it to compute the NI for the next iteration.
    if isfield(cvar, 'Im') && ~mp.ctrl.flagUseModel
        out.IrawScoreHist(Itr+1) = mean(cvar.Im(mp.Fend.score.maskBool));
        out.IrawCorrHist(Itr+1) = mean(cvar.Im(mp.Fend.corr.maskBool));
        out.InormHist(Itr+1) = out.IrawCorrHist(Itr+1);
        if mp.flagFiber
            out.InormFiberHist(Itr+1) = cvar.Ifiber;
        end
    end
    
end


function Esim = compute_simulated_efield_for_delta_efield_plot(mp)

    %--Model-based estimate for comparing Delta E (1st star only)
    modvar = ModelVariables;
    modvar.whichSource = 'star';
    modvar.starIndex = 1; % 1ST STAR ONLY
    for si = mp.Nsbp:-1:1
        modvar.sbpIndex = si;
        Etemp = model_compact(mp, modvar);
        Esim(:, si) = Etemp(mp.Fend.corr.maskBool);
    end
end


function [out, hProgress] = plot_wfsc_progress(mp, out, ev, hProgress, Itr, ImSimOffaxis)

    if strcmpi(mp.estimator, 'iefc')
        ev.Eest = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
    end

    Im = ev.Im;
    if any(mp.dm_ind == 1); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm); else; DM1surf = zeros(mp.dm1.compact.Ndm); end
    if any(mp.dm_ind == 2); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm); else; DM2surf = zeros(mp.dm2.compact.Ndm); end
    
    % Add open loop contrast history to out variable.
    if strcmpi(mp.estimator,'ekf_maintenance') 
        out.IOLScoreHist = ev.IOLScoreHist;
    end
    
    if isfield(mp, 'testbed')
        out.InormHist_tb.total = out.InormHist; 
        Im_tb.Im = Im;
        Im_tb.E = zeros([size(Im), mp.Nsbp]);
        Im_tb.Iinco = zeros([size(Im), mp.Nsbp]);
        if ~strcmpi(mp.estimator, 'perfect')
            for si = 1:mp.Nsbp
                tmp = zeros(size(Im));
                tmp(mp.Fend.corr.maskBool) = ev.Eest(:, si);
                Im_tb.E(:, :, si) = tmp; % modulated component 
 
                tmp = zeros(size(Im));
                tmp(mp.Fend.corr.maskBool) = ev.IincoEst(:, si);
                Im_tb.Iinco(:, :, si) = tmp; % unmodulated component 

                out.InormHist_tb.mod(Itr, si) = mean(abs(ev.Eest(:, si)).^2);
                out.InormHist_tb.unmod(Itr, si) = mean(ev.IincoEst(:, si));

                Im_tb.ev = ev; % Passing the probing structure so I can save it
            end
            clear tmp;
        else
            out.InormHist_tb.mod = NaN(Itr, mp.Nsbp);
            out.InormHist_tb.unmod = NaN(Itr, mp.Nsbp);
        end
        hProgress = falco_plot_progress_testbed(hProgress, mp, Itr, out.InormHist_tb, Im_tb, DM1surf, DM2surf);

    else
        if mp.flagFiber
            hProgress = falco_plot_progress(hProgress, mp, Itr, out.InormFiberHist, Im, DM1surf, DM2surf, ImSimOffaxis);
        else
            hProgress = falco_plot_progress(hProgress, mp, Itr, out.InormHist, Im, DM1surf, DM2surf, ImSimOffaxis);
        end

        if ~mp.flagFiber
            out.InormHist_tb.total = out.InormHist; 
            Im_tb.Im = Im;
            Im_tb.E = zeros([size(Im), mp.Nsbp]);
            Im_tb.Iinco = zeros([size(Im), mp.Nsbp]);

            for si = 1:mp.Nsbp
                    tmp = zeros(size(Im));
                    tmp(mp.Fend.corr.maskBool) = ev.Eest(:, si);
                    Im_tb.E(:, :, si) = tmp; % modulated component 

                    tmp = zeros(size(Im));
                    tmp(mp.Fend.corr.maskBool) = ev.IincoEst(:, si);
                    Im_tb.Iinco(:, :, si) = tmp; % unmodulated component 

                    out.InormHist_tb.mod(Itr, si) = mean(abs(ev.Eest(:, si)).^2);
                    out.InormHist_tb.unmod(Itr, si) = mean(ev.IincoEst(:, si));

                    Im_tb.ev = ev; % Passing the probing structure so I can save it
            end
        end
        
    end
end
