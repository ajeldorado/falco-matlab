% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Run a wavefront estimation and control loop.

function [mp, out] = falco_wfsc_loop(mp, out)

for Itr = 1:mp.Nitr
    
    %% Bookkeeping
    fprintf(['WFSC Iteration: ' num2str(Itr) '/' num2str(mp.Nitr) '\n' ]);
    
    if mp.flagSim
        fprintf('Zernike modes used in this Jacobian:\t');fprintf('%d ',mp.jac.zerns); fprintf('\n');
    end
    
    ev.Itr = Itr;
    cvar.Itr = Itr;
    
    % Updated DM info
    if strcmpi(mp.controller, 'plannedefc')
        mp.dm_ind = mp.dm_ind_sched{Itr}; % Change which DMs are used
    end
    fprintf('DMs to be used in this iteration = [')
    for jj = 1:length(mp.dm_ind)
        fprintf(' %d', mp.dm_ind(jj)); 
    end
    fprintf(' ]\n')
    
    out = store_dm_command_history(mp, out, Itr);

    %% Normalization and throughput calculations
    
    mp = falco_compute_PSF_norm_factor(mp);
    
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
    out.thput(Itr) = thput;
    mp.thput_vec(Itr) = max(thput); % note: max() needed when mp.flagFiber==true

    %% Control Jacobian

    mp = falco_set_jacobian_weights(mp); 
    
    relinearizeThisIter = any(mp.relinItrVec == Itr);
    if (Itr == 1) || relinearizeThisIter
        jacStruct =  model_Jacobian(mp);
    end
    
    [mp, jacStruct] = falco_ctrl_cull_weak_actuators(mp, cvar, jacStruct);

    % Load the improved Jacobian if using the E-M technique
    if mp.flagUseLearnedJac
        jacStructLearned = load('jacStructLearned.mat');
        if any(mp.dm_ind == 1);  jacStruct.G1 = jacStructLearned.G1;  end
        if any(mp.dm_ind == 1);  jacStruct.G2 = jacStructLearned.G2;  end
    end

    %% Wavefront Estimation
    
    if (Itr > 1); EestPrev = ev.Eest; end % save previous estimate for Delta E plot
    ev = falco_est(mp, ev, jacStruct);
    
    out = store_intensities(mp, out, ev, Itr);
    
    %% Progress plots (PSF, NI, and DM surfaces)

    if Itr == 1; hProgress.master = 1; end % initialize the handle
    [out, hProgress] = plot_wfsc_progress(mp, out, ev, hProgress, Itr, ImSimOffaxis);
    
    %% Plot the expected and measured delta E-fields
    
    if (Itr > 1); EsimPrev = Esim; end % save previous value for Delta E plot
    Esim = compute_simulated_efield_for_delta_efield_plot(mp);
    if Itr > 1
        out = falco_plot_DeltaE(mp, out, ev.Eest, EestPrev, Esim, EsimPrev, Itr);
    end
    
    %% Compute and Plot the Singular Mode Spectrum of the Electric Field
    if mp.flagSVD
        out = falco_plot_singular_mode_spectrum_of_Efield(mp, out, jacStruct, ev.Eest, Itr);
    end
        
    %% Wavefront Control
    
    jacStruct = falco_apply_spatial_weighting_to_Jacobian(mp, jacStruct);

    cvar.Eest = ev.Eest;
    cvar.NeleAll = mp.dm1.Nele + mp.dm2.Nele + mp.dm3.Nele + mp.dm4.Nele + mp.dm5.Nele + mp.dm6.Nele + mp.dm7.Nele + mp.dm8.Nele + mp.dm9.Nele; %--Number of total actuators used 
    [mp, cvar] = falco_ctrl(mp, cvar, jacStruct);
    
    % Save data to 'out'
    out.log10regHist(Itr) = cvar.log10regUsed;
    if isfield(cvar, 'Im') && ~mp.ctrl.flagUseModel
        out.InormHist(Itr+1) = mean(cvar.Im(mp.Fend.corr.maskBool));
        out.IrawCorrHist(Itr+1) = mean(cvar.Im(mp.Fend.corr.maskBool));
        out.IrawScoreHist(Itr+1) = mean(cvar.Im(mp.Fend.score.maskBool));
    end
    
    %--Enforce constraints on DM commands 
    if any(mp.dm_ind == 1); mp.dm1 = falco_enforce_dm_constraints(mp.dm1); end
    if any(mp.dm_ind == 2); mp.dm2 = falco_enforce_dm_constraints(mp.dm2); end
    
    %--Update DM actuator gains for new voltages
    if any(mp.dm_ind == 1); mp.dm1 = falco_update_dm_gain_map(mp.dm1); end
    if any(mp.dm_ind == 2); mp.dm2 = falco_update_dm_gain_map(mp.dm2); end
    
    %% Report Various Stats

    out = falco_compute_dm_stats(mp, out, Itr);

    %--Calculate sensitivities to small Zernike phase aberrations at entrance pupil.
    if ~isempty(mp.eval.Rsens) && ~isempty(mp.eval.indsZnoll)
        out.Zsens(:, :, Itr) = falco_get_Zernike_sensitivities(mp);
    end

    %--REPORTING NORMALIZED INTENSITY
    if isfield(cvar, 'cMin') && mp.ctrl.flagUseModel == false
        fprintf('Prev and New Measured NI:\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
            out.InormHist(Itr), out.InormHist(Itr+1), out.InormHist(Itr)/out.InormHist(Itr+1) );
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

end %--END OF ESTIMATION + CONTROL LOOP


%% Update 'out' structure and progress plot one last time
Itr = Itr + 1;

out = store_dm_command_history(mp, out, Itr);

[mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
out.thput(Itr) = thput;
mp.thput_vec(Itr) = max(thput); % max() used for if mp.flagFiber==true

if isfield(cvar, 'Im')
    ev.Im = cvar.Im;
    [out, hProgress] = plot_wfsc_progress(mp, out, ev, hProgress, Itr, ImSimOffaxis);
end
%% Save out an abridged workspace

fnSnippet = [mp.path.config filesep mp.runLabel,'_snippet.mat'];
fprintf('\nSaving data snippet to:\n\t%s\n', fnSnippet)
save(fnSnippet, 'out');
fprintf('...done.\n')

%% Save out the data from the workspace
if mp.flagSaveWS
    clear ev cvar G* h* jacStruct; % Save a ton of space when storing the workspace

    % Don't bother saving the large 2-D, floating point maps in the workspace (they take up too much space)
    mp.P1.full.mask=1; mp.P1.compact.mask=1;
    mp.P3.full.mask=1; mp.P3.compact.mask=1;
    mp.P4.full.mask=1; mp.P4.compact.mask=1;
    mp.F3.full.mask=1; mp.F3.compact.mask=1;

    mp.P1.full.E = 1; mp.P1.compact.E = 1; mp.Eplanet = 1; 
    mp.dm1.full.mask = 1; mp.dm1.compact.mask = 1; mp.dm2.full.mask = 1; mp.dm2.compact.mask = 1;
    mp.complexTransFull = 1; mp.complexTransCompact = 1;

    mp.dm1.compact.inf_datacube = 0;
    mp.dm2.compact.inf_datacube = 0;
    mp.dm8.compact.inf_datacube = 0;
    mp.dm9.compact.inf_datacube = 0;
    mp.dm8.inf_datacube = 0;
    mp.dm9.inf_datacube = 0;

    fnAll = [mp.path.ws mp.runLabel, '_all.mat'];
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


function Esim = compute_simulated_efield_for_delta_efield_plot(mp)

    %--Model-based estimate for comparing Delta E (1st star only)
    modvar.whichSource = 'star';
    modvar.starIndex = 1; % 1ST STAR ONLY
    for si = mp.Nsbp:-1:1
        modvar.sbpIndex = si;
        Etemp = model_compact(mp, modvar);
        Esim(:, si) = Etemp(mp.Fend.corr.maskBool);
    end
end


function out = store_intensities(mp, out, ev, Itr)

    %% Calculate the average measured, coherent, and incoherent intensities
    
    % Apply subband weights and then sum over subbands
    Iest = abs(ev.Eest).^2;
    Iinco = ev.IincoEst;
    for iMode = 1:mp.jac.Nmode
        iSubband = mp.jac.sbp_inds(iMode);
        Iest(:, iMode) = mp.sbp_weights(iSubband) * Iest(:, iMode);
        Iinco(:, iMode) = mp.sbp_weights(iSubband) * Iinco(:, iMode);
    end
    IestAllBands = sum(Iest, 2);
    IincoAllBands = sum(Iinco, 2);
    
    
    % Put intensities back into 2-D arrays to use correct indexing of scoring region.
    % Modulated
    Iest2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Iest2D(mp.Fend.corr.maskBool) = IestAllBands(:);
    out.IestScoreHist(Itr) = mean(Iest2D(mp.Fend.score.maskBool));
    out.IestCorrHist(Itr) = mean(Iest2D(mp.Fend.corr.maskBool));
    % Unmodulated
    Iinco2D = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    Iinco2D(mp.Fend.corr.maskBool) = IincoAllBands(:);
    out.IincoScoreHist(Itr) = mean(Iinco2D(mp.Fend.score.maskBool));
    out.IincoCorrHist(Itr) = mean(Iinco2D(mp.Fend.corr.maskBool));

    % Measured
    out.IrawScoreHist(Itr) = mean(ev.Im(mp.Fend.score.maskBool));
    out.IrawCorrHist(Itr) = mean(ev.Im(mp.Fend.corr.maskBool));
    out.InormHist(Itr) = out.IrawCorrHist(Itr); % InormHist is a vestigial variable
    
    %% Calculate the measured, coherent, and incoherent intensities by subband
    
    % measured intensities
    for iSubband = 1:mp.Nsbp
        imageMeas = squeeze(ev.imageArray(:, :, 1, iSubband));
        out.normIntMeasCorr(Itr, iSubband) = mean(imageMeas(mp.Fend.corr.maskBool));
        out.normIntMeasScore(Itr, iSubband) = mean(imageMeas(mp.Fend.score.maskBool));
        clear im        
    end
    
    % estimated
    for iMode = 1:(mp.Nsbp*mp.compact.star.count)
        
        imageModVec = abs(ev.Eest(:, iMode)).^2;
        imageUnmodVec = ev.IincoEst(:, iMode);
        
        out.normIntModCorr(Itr, iMode) = mean(imageModVec);
        out.normIntModScore(Itr, iMode) = mean(imageModVec(mp.Fend.scoreInCorr));
        
        out.normIntUnmodCorr(Itr, iMode) = mean(imageUnmodVec);
        out.normIntModScore(Itr, iMode) = mean(imageUnmodVec(mp.Fend.scoreInCorr));
    end
    
end


function [out, hProgress] = plot_wfsc_progress(mp, out, ev, hProgress, Itr, ImSimOffaxis)

    Im = ev.Im;
    if any(mp.dm_ind == 1); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm); else; DM1surf = zeros(mp.dm1.compact.Ndm); end
    if any(mp.dm_ind == 2); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm); else; DM2surf = zeros(mp.dm2.compact.Ndm); end
    
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
        hProgress = falco_plot_progress(hProgress, mp, Itr, out.InormHist, Im, DM1surf, DM2surf, ImSimOffaxis);
    end
end
