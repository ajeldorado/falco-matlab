% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Run a wavefront estimation and control loop.

function [mp, out] = falco_wfsc_loop(mp, out)

% Take initial broadband image 
Im = falco_get_summed_image(mp);

% Begin the Estimation+Control Iterations
for Itr = 1:mp.Nitr

    %% Apply DM constraints now. Can't do within DM surface generator if calling a PROPER model. 
    if(any(mp.dm_ind==1));  mp.dm1 = falco_enforce_dm_constraints(mp.dm1);  end
    if(any(mp.dm_ind==2));  mp.dm2 = falco_enforce_dm_constraints(mp.dm2);  end
    
    %%
    %--Start of new estimation+control iteration
    fprintf(['Iteration: ' num2str(Itr) '/' num2str(mp.Nitr) '\n' ]);

    %--Re-compute the starlight normalization factor for the compact and full models (to convert images to normalized intensity). No tip/tilt necessary.
    mp = falco_get_PSF_norm_factor(mp);
    
    %% Updated DM data
    %--Change the selected DMs if using the scheduled EFC controller
    switch lower(mp.controller)
        case{'plannedefc'} 
            mp.dm_ind = mp.dm_ind_sched{Itr};
    end
    %--Report which DMs are used in this iteration
    fprintf('DMs to be used in this iteration = ['); for jj = 1:length(mp.dm_ind); fprintf(' %d',mp.dm_ind(jj)); end; fprintf(' ]\n');

    %--Fill in History of DM commands to Store
    if(isfield(mp,'dm1')); if(isfield(mp.dm1,'V'));  out.dm1.Vall(:,:,Itr) = mp.dm1.V;  end;  end
    if(isfield(mp,'dm2')); if(isfield(mp.dm2,'V'));  out.dm2.Vall(:,:,Itr) = mp.dm2.V;  end;  end
    if(isfield(mp,'dm8')); if(isfield(mp.dm8,'V'));  out.dm8.Vall(:,Itr) = mp.dm8.V(:);  end;  end
    if(isfield(mp,'dm9')); if(isfield(mp.dm9,'V'));  out.dm9.Vall(:,Itr) = mp.dm9.V(:);  end;  end

    %--Compute the DM surfaces
    if(any(mp.dm_ind==1)); DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm); else; DM1surf = zeros(mp.dm1.compact.Ndm); end
    if(any(mp.dm_ind==2)); DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm); else; DM2surf = zeros(mp.dm2.compact.Ndm); end

    %% Throughput and Normalized Intensity
    %--Calculate the core throughput (at higher resolution to be more accurate)
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp);
    if(mp.flagFiber)
        mp.thput_vec(Itr) = max(thput);
    else
        mp.thput_vec(Itr) = thput;
    end
    
    %--Compute the current normalized intensity level
    out.InormHist(Itr) = mean(Im(mp.Fend.corr.maskBool));
    
    %% Updated selection of Zernike modes targeted by the controller
    %--Decide with Zernike modes to include in the Jacobian
    if(Itr==1); mp.jac.zerns0 = mp.jac.zerns; end
    fprintf('Zernike modes used in this Jacobian:\t'); fprintf('%d ',mp.jac.zerns); fprintf('\n');
    
    %--Re-compute the Jacobian weights
    mp = falco_set_jacobian_weights(mp); 

    %% Actuator Culling: Initialization of Flag and Which Actuators

    %--If new actuators are added, perform a new cull of actuators.
    if(Itr==1)
        cvar.flagCullAct = true;
    else
        if(isfield(mp,'dm_ind_sched'))
            cvar.flagCullAct = ~isequal(mp.dm_ind_sched{Itr}, mp.dm_ind_sched{Itr-1});
        else
            cvar.flagCullAct = false;
        end
    end
    mp.flagCullActHist(Itr) = cvar.flagCullAct;

    %--Before performing new cull, include all actuators again
    if(cvar.flagCullAct)
        %--Re-include all actuators in the basis set. Need act_ele to be a column vector.
        if(any(mp.dm_ind==1)); mp.dm1.act_ele = (1:mp.dm1.NactTotal).'; end
        if(any(mp.dm_ind==2)); mp.dm2.act_ele = (1:mp.dm2.NactTotal).'; end
        if(any(mp.dm_ind==8)); mp.dm8.act_ele = (1:mp.dm8.NactTotal).'; end
        if(any(mp.dm_ind==9)); mp.dm9.act_ele = (1:mp.dm9.NactTotal).'; end
        %--Update the number of elements used per DM
        if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); else; mp.dm1.Nele = 0; end
        if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); else; mp.dm2.Nele = 0; end
        if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); else; mp.dm8.Nele = 0; end
        if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); else; mp.dm9.Nele = 0; end
    end

    %% Compute the control Jacobians for each DM
    
    %--Relinearize about the DMs only at the iteration numbers in mp.relinItrVec.
    if any(mp.relinItrVec == Itr)
        cvar.flagRelin = true;
    else
        cvar.flagRelin = false;
    end
    
    if (Itr == 1) || cvar.flagRelin
        jacStruct =  model_Jacobian(mp); %--Get structure containing Jacobians
    end

    %% Cull weak actuators, but only if(cvar.flagCullAct && cvar.flagRelin)
    [mp, jacStruct] = falco_ctrl_cull(mp, cvar, jacStruct);

    %% Load the improved Jacobian if using the E-M technique
    if mp.flagUseLearnedJac
        jacStructLearned = load('jacStructLearned.mat');
        if(any(mp.dm_ind==1));  jacStruct.G1 = jacStructLearned.G1;  end
        if(any(mp.dm_ind==1));  jacStruct.G2 = jacStructLearned.G2;  end
    end

    %% Wavefront Estimation
    
    %--Previous iteration's E-field for computing and plotting Delta E
    if(Itr > 1)
        EestPrev = Eest;
        EsimPrev = Esim;
    end

    %--Model-based estimate for comparing Delta E (1st star only)
    modvar.whichSource = 'star';
    modvar.starIndex = 1; % 1ST STAR ONLY
    for si = mp.Nsbp:-1:1
        modvar.sbpIndex = si;
        Etemp = model_compact(mp, modvar);
        Esim(:, si) = Etemp(mp.Fend.corr.maskBool);
    end
    clear modvar Etemp
    
    ev.Itr = Itr;
    
    ev = falco_est(mp, ev, jacStruct);
    
    Eest = ev.Eest;
    IincoEst = ev.IincoEst;
    Im = ev.Im;
    out.IestCorrHist(Itr) = ev.IestCorrMean;
    out.IestScoreHist(Itr) = ev.IestScoreMean;
    out.IincoCorrHist(Itr) = ev.IincoCorrMean;
    out.IincoScoreHist(Itr) = ev.IincoScoreMean;
    
    %% Plot the updates to the DMs and PSF
    if Itr == 1; hProgress.master = 1; end %--dummy value to intialize the progress plot's handle
    if isfield(mp, 'testbed')
        out.InormHist_tb.total = out.InormHist; 
        Im_tb.Im = Im;
        Im_tb.E = zeros([size(Im), mp.Nsbp]);
        Im_tb.Iinco = zeros([size(Im), mp.Nsbp]);
        if ~strcmpi(mp.estimator, 'perfect')
            for si = 1:mp.Nsbp
                tmp = zeros(size(Im));
                tmp(mp.Fend.corr.mask) = Eest(:, si);
                Im_tb.E(:, :, si) = tmp; % modulated component 
 
                tmp = zeros(size(Im));
                tmp(mp.Fend.corr.mask) = IincoEst(:, si);
                Im_tb.Iinco(:, :, si) = tmp; % unmodulated component 

                out.InormHist_tb.mod(Itr, si) = mean(abs(Eest(:, si)).^2);
                out.InormHist_tb.unmod(Itr, si) = mean(IincoEst(:, si));

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
    
    %% Plot the expected and measured delta E-fields
    if Itr > 1
        out = falco_plot_DeltaE(mp, out, Eest, EestPrev, Esim, EsimPrev, Itr);
    end
    
    %% Compute and Plot the Singular Mode Spectrum of the Electric Field
    if mp.flagSVD
        out = falco_plot_singular_mode_spectrum_of_Efield(mp, out, jacStruct, Eest, Itr);
    end
    
    %% Add spatially-dependent (and star-dependent) weighting to the control Jacobians

    for iStar = 1:mp.compact.star.count
        if(any(mp.dm_ind==1))
            jacStruct.G1(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G1(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm1.Nele, mp.jac.NmodePerStar]);
        end
        if(any(mp.dm_ind==2))
            jacStruct.G2(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G2(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm2.Nele, mp.jac.NmodePerStar]);
        end
        if(any(mp.dm_ind==8))
            jacStruct.G8(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G8(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm8.Nele, mp.jac.NmodePerStar]);
        end
        if(any(mp.dm_ind==9))
            jacStruct.G9(:, :, mp.jac.star_inds == iStar) = ...
                jacStruct.G9(:, :, mp.jac.star_inds == iStar) .* repmat(mp.WspatialVec(:, iStar), [1, mp.dm9.Nele, mp.jac.NmodePerStar]);
        end
    end

    %--Compute the number of total actuators for all DMs used. 
    cvar.NeleAll = mp.dm1.Nele + mp.dm2.Nele + mp.dm3.Nele + mp.dm4.Nele + mp.dm5.Nele + mp.dm6.Nele + mp.dm7.Nele + mp.dm8.Nele + mp.dm9.Nele; %--Number of total actuators used 
    
    %% Wavefront Control
    
    cvar.Eest = Eest;
    cvar.Itr = Itr;
    cvar.out.InormHist = out.InormHist(Itr);
    [mp, cvar] = falco_ctrl(mp, cvar, jacStruct);
    if isfield(cvar, 'Im') && ~mp.ctrl.flagUseModel
        Im = cvar.Im; 
    end
    
    %--Enforce constraints on DM commands 
    % (not needed here--just done here for stats and plotting)
    if(any(mp.dm_ind==1)); mp.dm1 = falco_enforce_dm_constraints(mp.dm1); end
    if(any(mp.dm_ind==2)); mp.dm2 = falco_enforce_dm_constraints(mp.dm2); end
    
    %--Update DM actuator gains for new voltages
    if(any(mp.dm_ind==1)); mp.dm1 = falco_update_dm_gain_map(mp.dm1); end
    if(any(mp.dm_ind==2)); mp.dm2 = falco_update_dm_gain_map(mp.dm2); end
    
    %--Save out regularization used.
    out.log10regHist(Itr) = cvar.log10regUsed; 

    %% Report and Store Various Stats

    out = falco_compute_dm_stats(mp, out, Itr);

    %--Calculate sensitivities to 1nm RMS of Zernike phase aberrations at entrance pupil.
    if( isempty(mp.eval.Rsens)==false || isempty(mp.eval.indsZnoll)==false )
        out.Zsens(:,:,Itr) = falco_get_Zernike_sensitivities(mp);
    end

    %--REPORTING NORMALIZED INTENSITY
    if isfield(cvar, 'cMin') && mp.ctrl.flagUseModel == false
        out.InormHist(Itr+1) = cvar.cMin;
        fprintf('Prev and New Measured NI:\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n\n',...
            out.InormHist(Itr), out.InormHist(Itr+1), out.InormHist(Itr)/out.InormHist(Itr+1) );
        if ~mp.flagSim
            fprintf('\n\n');
        end
    else
        fprintf('Previous Measured NI:\t\t\t %.2e \n\n', out.InormHist(Itr))
    end

    %--Save out DM commands after each iteration in case the trial crashes part way through.
    if(mp.flagSaveEachItr)
        fprintf('Saving DM commands for this iteration...')
        if(any(mp.dm_ind==1)); DM1V = mp.dm1.V; else; DM1V = 0; end
        if(any(mp.dm_ind==2)); DM2V = mp.dm2.V; else; DM2V = 0; end
        if(any(mp.dm_ind==8)); DM8V = mp.dm8.V; else; DM8V = 0; end
        if(any(mp.dm_ind==9)); DM9V = mp.dm9.V; else; DM9V = 0; end
        Nitr = mp.Nitr;
        thput_vec = mp.thput_vec;
        fnWS = sprintf('%sws_%s_Iter%dof%d.mat',mp.path.wsInProgress, mp.runLabel, Itr, mp.Nitr);
        save(fnWS,'Nitr','Itr','DM1V','DM2V','DM8V','DM9V','out.InormHist','thput_vec','out')
        fprintf('done.\n\n')
    end

    %% SAVE THE TRAINING DATA OR RUN THE E-M Algorithm
    if(mp.flagTrainModel)
        ev.Itr = Itr;
        mp = falco_train_model(mp,ev);
    end

end %--END OF ESTIMATION + CONTROL LOOP

%% Update progress plot one last time
Itr = Itr + 1;

%--Compute the DM surfaces
if(any(mp.dm_ind==1)); DM1surf =  falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  end
if(any(mp.dm_ind==2)); DM2surf =  falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  end

%--Data to store
if(any(mp.dm_ind==1)); out.dm1.Vall(:,:,Itr) = mp.dm1.V; end
if(any(mp.dm_ind==2)); out.dm2.Vall(:,:,Itr) = mp.dm2.V; end
if(any(mp.dm_ind==5)); out.dm5.Vall(:,:,Itr) = mp.dm5.V; end
if(any(mp.dm_ind==8)); out.dm8.Vall(:,Itr) = mp.dm8.V; end
if(any(mp.dm_ind==9)); out.dm9.Vall(:,Itr) = mp.dm9.V; end

%--Calculate the core throughput (at higher resolution to be more accurate)
[mp,thput,ImSimOffaxis] = falco_compute_thput(mp);
if(mp.flagFiber)
    mp.thput_vec(Itr) = max(thput);
else
    mp.thput_vec(Itr) = thput; %--record keeping
end

if(isfield(mp,'testbed') )
    out.InormHist_tb.total = out.InormHist; 
    Im_tb.Im = Im;
    Im_tb.E = zeros([size(Im),mp.Nsbp]);
    Im_tb.Iinco = zeros([size(Im),mp.Nsbp]);
    if(~strcmpi(mp.estimator,'perfect') )
        for si = 1:mp.Nsbp
            tmp = zeros(size(Im));
            tmp(mp.Fend.corr.mask) = Eest(:,si);
            Im_tb.E(:,:,si) = tmp; % modulated component 

            tmp = zeros(size(Im));
            tmp(mp.Fend.corr.mask) = IincoEst(:,si);
            Im_tb.Iinco(:,:,si) = tmp; % unmodulated component 

            out.InormHist_tb.mod(Itr,si) = mean(abs(Eest(:,si)).^2);
            out.InormHist_tb.unmod(Itr,si) = mean(IincoEst(:,si));

            Im_tb.ev = ev; % Passing the probing structure so I can save it
        end
        clear tmp;
    else
        out.InormHist_tb.mod = NaN(Itr,mp.Nsbp);
        out.InormHist_tb.unmod = NaN(Itr,mp.Nsbp);
    end
    hProgress = falco_plot_progress_testbed(hProgress,mp,Itr,out.InormHist_tb,Im_tb,DM1surf,DM2surf);
else
    hProgress = falco_plot_progress(hProgress,mp,Itr,out.InormHist,Im,DM1surf,DM2surf,ImSimOffaxis);
end

%% Save the final DM commands separately for faster reference
if(isfield(mp,'dm1')); if(isfield(mp.dm1,'V')); out.DM1V = mp.dm1.V; end; end
if(isfield(mp,'dm2')); if(isfield(mp.dm2,'V')); out.DM2V = mp.dm2.V; end; end
switch upper(mp.coro)
    case{'HLC','EHLC'}
        if(isfield(mp.dm8,'V')); out.DM8V = mp.dm8.V;  end
        if(isfield(mp.dm9,'V')); out.DM9V = mp.dm9.V;  end
end

%% Save out an abridged workspace

%--Variables to save out:
% contrast vs iter
% regularization history
%  DM1surf,DM1V, DM2surf,DM2V, DM8surf,DM9surf, fpm sampling, base pmgi thickness, base nickel thickness, dm_tilts, aoi, ...
% to reproduce your basic design results of NI, throughput, etc..

out.thput = mp.thput_vec;
out.Nitr = mp.Nitr;

fnOut = [mp.path.config filesep mp.runLabel,'_snippet.mat'];

fprintf('\nSaving abridged workspace to file:\n\t%s\n',fnOut)
save(fnOut,'out');
fprintf('...done.\n\n')

%% Save out the data from the workspace
if(mp.flagSaveWS)
    clear cvar G* h* jacStruct; % Save a ton of space when storing the workspace
    clear Gcomplex Gall Eall Eri Gri U S EriPrime IriPrime 

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

    fnAll = [mp.path.ws mp.runLabel,'_all.mat'];
    disp(['Saving entire workspace to file ' fnAll '...'])
    save(fnAll);
    fprintf('done.\n\n')
else
    disp('Entire workspace NOT saved because mp.flagSaveWS==false')
end

end %--END OF main FUNCTION
