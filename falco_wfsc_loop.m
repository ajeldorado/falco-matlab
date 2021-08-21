% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to run wavefront a estimation and control loop for various 
% types of coronagraphs.

function [mp, out] = falco_wfsc_loop(mp, out)

if(isfield(mp,'testbed') )
    bench = mp.bench;
end
%% Initializations of Arrays for Data Storage 
InormHist = zeros(mp.Nitr,1); % Measured, mean raw contrast in scoring region of dark hole.
InormSMFHist = zeros(mp.Nitr,1);
%--Raw contrast (broadband)

%% Plot the pupil masks
% if(mp.flagPlot); figure(101); imagesc(mp.P1.full.mask);axis image; colorbar; title('pupil');drawnow; end
% if(mp.flagPlot && (length(mp.P4.full.mask)==length(mp.P1.full.mask))); figure(102); imagesc(mp.P4.full.mask);axis image; colorbar; title('Lyot stop');drawnow; end
% if(mp.flagPlot && isfield(mp,'P3.full.mask')); figure(103); imagesc(padOrCropEven(mp.P1.full.mask,mp.P3.full.Narr).*mp.P3.full.mask);axis image; colorbar; drawnow; end

%% Take initial broadband image 
if mp.flagFiber
    [Im, Ifiber,Ifiber_arr] = falco_get_summed_image(mp);
else
    [Im,~,Imean_arr] = falco_get_summed_image(mp);
end

%%
sz = size(Im);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin the Correction Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Itr=1:mp.Nitr

    %--Start of new estimation+control iteration
    fprintf(['Iteration: ' num2str(Itr) '/' num2str(mp.Nitr) '\n' ]);
    
    if(~mp.flagSim && any(mp.search4OffaxisSMF_it==Itr))
        [ x_smf_out, y_smf_out ] = hcst_fiu_findOffaxisSMFwProbe(bench, mp);
        mp.est.probe.Xloc = x_smf_out;
        mp.est.probe.Yloc = y_smf_out;
        mp.Fend.x_fiber = mp.est.probe.Xloc;
        mp.Fend.y_fiber = mp.est.probe.Yloc;
        % We have to re-compute the Jacobian
        mp = falco_configure_fiber_dark_hole(mp);
        mp.relinItrVec = [mp.relinItrVec,Itr];
    end
    
    % new evaluation of jitter instance (1/2 pre-estimation)
    if mp.flagSim && mp.flagJitter && Itr~=1
        mp_newFiberMode = mp;
        mp_newFiberMode.Fend.x_fiber = mp.Fend.x_fiber + randn(1,1)*mp.Fend.jitt_amp;
        mp_newFiberMode.Fend.y_fiber = [mp.Fend.y_fiber] + randn(1,1)*mp.Fend.jitt_amp;
        mp = falco_configure_fiber_dark_hole(mp_newFiberMode);
    end
    
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
    if(isfield(mp,'dm5')); if(isfield(mp.dm5,'V'));  out.dm5.Vall(:,:,Itr) = mp.dm5.V;  end;  end
    if(isfield(mp,'dm8')); if(isfield(mp.dm8,'V'));  out.dm8.Vall(:,Itr) = mp.dm8.V(:);  end;  end
    if(isfield(mp,'dm9')); if(isfield(mp.dm9,'V'));  out.dm9.Vall(:,Itr) = mp.dm9.V(:);  end;  end

    %--Compute the DM surfaces
    if(any(mp.dm_ind==1)); DM1surf =  falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  else; DM1surf = zeros(mp.dm1.compact.Ndm);  end
    if(any(mp.dm_ind==2)); DM2surf =  falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  else; DM2surf = zeros(mp.dm2.compact.Ndm);    end
    if(any(mp.dm_ind==5))
        DM5surf =  falco_gen_dm_surf(mp.dm5, mp.dm5.compact.dx, mp.dm5.compact.Ndm); 
        figure(325); imagesc(DM5surf); axis xy equal tight; colorbar; drawnow;
        figure(326); imagesc(mp.dm5.V); axis xy equal tight; colorbar; drawnow;
    end

    switch upper(mp.coro)
        case{'EHLC'}
            mp.DM8surf = falco_gen_EHLC_FPM_surf_from_cube(mp.dm8,'compact'); %--Metal layer profile [m]
            mp.DM9surf = falco_gen_EHLC_FPM_surf_from_cube(mp.dm9,'compact'); %--Dielectric layer profile [m]
        case{'HLC','APHLC','SPHLC'}
            switch lower(mp.layout)
                case 'fourier'
                    mp.DM8surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm8,'compact'); %--Metal layer profile [m]
                    mp.DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9,'compact'); %--Dielectric layer profile [m]
            end
        case{'FOHLC'}
            mp.DM8amp = falco_gen_HLC_FPM_amplitude_from_cube(mp.dm8,'compact'); %--FPM amplitude transmission [amplitude]
            mp.DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9,'compact'); %--FPM phase shift in transmission [m]    
    end

    %% Updated plot and reporting
    %--Calculate the core throughput (at higher resolution to be more accurate)
    [mp,thput,ImSimOffaxis] = falco_compute_thput(mp);
    if(mp.flagFiber)
        mp.thput_vec(Itr) = max(thput);
    else
        mp.thput_vec(Itr) = thput; %--record keeping
    end
    
    %--Compute the current contrast level
    InormHist(Itr) = mean(Im(mp.Fend.corr.maskBool));
    if(mp.flagFiber)
        InormSMFHist(Itr) = Ifiber;
    end
    %--Plot the updates to the DMs and PSF
    if(Itr==1); hProgress.master = 1; end %--dummy value to intialize the handle variable
    if(isfield(mp,'testbed') )
        if(mp.flagFiber)
            InormHist_tb.total = InormSMFHist; 
            InormHist_tb.Inorm_arr = Ifiber_arr;
        else
            InormHist_tb.total = InormHist; 
            InormHist_tb.Inorm_arr = Imean_arr;
        end
        
        Im_tb.Im = Im;
        Im_tb.E = zeros(size(Im));
        if(Itr>1)
            InormHist_tb.mod(Itr-1) = mean(abs(EfieldVec(:)).^2);
            InormHist_tb.unmod(Itr-1) = mean(IincoVec(:));
            Im_tb.E(mp.Fend.corr.mask) = EfieldVec(:,ceil(mp.Nsbp/2));
            InormHist_tb.beta(Itr-1) = out.log10regHist(Itr-1);
        else
            InormHist_tb.mod = NaN;
            InormHist_tb.unmod = NaN;
        end
        hProgress = falco_plot_progress_hcst(hProgress,mp,Itr,InormHist_tb,Im_tb,DM1surf,DM2surf);
        Im_tb.Im_prev = zeros(size(Im));
        Im_tb.Im_prev(mp.Fend.corr.mask) = Im(mp.Fend.corr.mask);
    else
        hProgress = falco_plot_progress(hProgress,mp,Itr,InormHist,Im,DM1surf,DM2surf,ImSimOffaxis);
    end
    
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
        if(any(mp.dm_ind==5)); mp.dm5.act_ele = (1:mp.dm5.NactTotal).'; end
        if(any(mp.dm_ind==8)); mp.dm8.act_ele = (1:mp.dm8.NactTotal).'; end
        if(any(mp.dm_ind==9)); mp.dm9.act_ele = (1:mp.dm9.NactTotal).'; end
        %--Update the number of elements used per DM
        if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); else; mp.dm1.Nele = 0; end
        if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); else; mp.dm2.Nele = 0; end
        if(any(mp.dm_ind==5)); mp.dm5.Nele = length(mp.dm5.act_ele); else; mp.dm5.Nele = 0; end
        if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); else; mp.dm8.Nele = 0; end
        if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); else; mp.dm9.Nele = 0; end
    end

    %% Compute the control Jacobians for each DM
    
    %--Relinearize about the DMs only at the iteration numbers in mp.relinItrVec.
    if(any(mp.relinItrVec==Itr))
        cvar.flagRelin=true;
    else
        cvar.flagRelin=false;
    end
    
%     if( (Itr==1) || cvar.flagRelin ) %% && (~mp.flagUseLearnedJac)
    if((Itr==1) && mp.use_lastJacStruc)
        load([mp.path.falco,'/data/jac/lastJacStruc.mat']);
    else
        if(cvar.flagRelin ) %% && (~mp.flagUseLearnedJac)
            jacStruct =  model_Jacobian(mp); %--Get structure containing Jacobians
            save([mp.path.falco,'/data/jac/lastJacStruc.mat'],'jacStruct');
        end
    end

    %% Cull weak actuators, but only if(cvar.flagCullAct && cvar.flagRelin)
    [mp, jacStruct] = falco_ctrl_cull(mp,cvar,jacStruct);

    %% Load the improved Jacobian if using the E-M technique
    if(mp.flagUseLearnedJac)
        if mp.est.flagUseTensorflow
            jacStructLearned = load('jacStructLearned.mat');
        end
        if mp.est.flagEM
            load('jacStructLearned.mat');
        end
        if(any(mp.dm_ind==1));  jacStruct.G1 = jacStructLearned.G1;  end
        if(any(mp.dm_ind==2));  jacStruct.G2 = jacStructLearned.G2;  end
    end

    %% Wavefront Estimation
    % First do a fine alignment of the FPM if necesary
    if(any(mp.fineAlignment_it==Itr))
        bench = hcstr_realignFPMAndRecenter(bench,mp);
    end

    switch lower(mp.estimator)
        case{'perfect'}
            EfieldVec  = falco_est_perfect_Efield_with_Zernikes(mp);
        case{'pwp-bp','pwp-kf'}
% 			if(mp.flagFiber && mp.flagLenslet)
			if(mp.flagFiber)
				if(mp.est.flagUseJac) %--Send in the Jacobian if true
					ev = falco_est_pairwise_probing_fiber(mp,jacStruct);
				else %--Otherwise don't pass the Jacobian
					ev = falco_est_pairwise_probing_fiber(mp);
				end
			else
				if(mp.est.flagUseJac) %--Send in the Jacobian if true
					ev = falco_est_pairwise_probing(mp,jacStruct);
				else %--Otherwise don't pass the Jacobian
					ev = falco_est_pairwise_probing(mp);
				end
            end
            
            EfieldVec = ev.Eest;
            IincoVec = ev.IincoEst;
    end
    
    % new evaluation of jitter instance (2/2 post-estimation)
    if mp.flagSim && mp.flagJitter && Itr~=1
        mp_newFiberMode = mp;
        mp_newFiberMode.Fend.x_fiber = mp.Fend.x_fiber + randn(1,1)*mp.Fend.jitt_amp;
        mp_newFiberMode.Fend.y_fiber = [mp.Fend.y_fiber] + randn(1,1)*mp.Fend.jitt_amp;
        mp = falco_configure_fiber_dark_hole(mp_newFiberMode);
    end

    %% Compute and Plot the Singular Mode Spectrum of the Control Jacobian

    if(mp.flagSVD)
        
        if(cvar.flagRelin)
            
            ii=1;
            Gcomplex = [jacStruct.G1(:,:,ii), jacStruct.G2(:,:,ii), jacStruct.G3(:,:,ii), jacStruct.G4(:,:,ii), jacStruct.G5(:,:,ii), jacStruct.G6(:,:,ii), jacStruct.G7(:,:,ii), jacStruct.G8(:,:,ii), jacStruct.G9(:,:,ii)];
            Gall = zeros(mp.jac.Nmode*size(Gcomplex,1),size(Gcomplex,2));
            Eall = zeros(mp.jac.Nmode*size(EfieldVec,1),1);

            for ii=1:mp.jac.Nmode
                N = size(Gcomplex,1);
                inds = (ii-1)*N+1:ii*N;
                Gcomplex = [jacStruct.G1(:,:,ii), jacStruct.G2(:,:,ii), jacStruct.G3(:,:,ii), jacStruct.G4(:,:,ii), jacStruct.G5(:,:,ii), jacStruct.G6(:,:,ii), jacStruct.G7(:,:,ii), jacStruct.G8(:,:,ii), jacStruct.G9(:,:,ii)];
                Gall(inds,:) = Gcomplex;
                Eall(inds) = EfieldVec(:,ii);
            end

            Eri = [real(Eall); imag(Eall)];
            alpha2 = max(diag(real(Gall'*Gall)));
            Gri = [real(Gall); imag(Gall)];
            [U,S,~] = svd(Gri,'econ');
            s = diag(S);
        else
           
            for ii=1:mp.jac.Nmode
                N = size(Gcomplex,1);
                inds = (ii-1)*N+1:ii*N;
                Eall(inds) = EfieldVec(:,ii);
            end
            Eri = [real(Eall); imag(Eall)];
        
        end
        
        EriPrime = U'*Eri;
        IriPrime = abs(EriPrime).^2;
    
        %--Save out for later analysis
        out.EforSpectra{Itr} = EriPrime;
        out.smspectra{Itr} = IriPrime;
        out.sm{Itr} = s;
        out.alpha2{Itr} = alpha2;
        
        if(mp.flagPlot)
            figure(401); 
            loglog(out.sm{Itr}.^2/out.alpha2{Itr},smooth(out.smspectra{Itr},31),'Linewidth',3,'Color',[0.3, 1-(0.2+Itr/mp.Nitr)/(1.3),1 ]);
            set(gca,'Fontsize',20); grid on; 
            set(gcf,'Color',[1 1 1]);
            title('Singular Mode Spectrum','Fontsize',20)
            xlim([1e-10, 2*max(s.^2/alpha2)])
            ylim([1e-12, 1e-0]) 
            drawnow;
            hold on;
        end
        
    end
    
    %% Add spatially-dependent weighting to the control Jacobians
    if(any(mp.dm_ind==1)); jacStruct.G1 = jacStruct.G1.*repmat(mp.WspatialVec,[1,mp.dm1.Nele,mp.jac.Nmode]); end
    if(any(mp.dm_ind==2)); jacStruct.G2 = jacStruct.G2.*repmat(mp.WspatialVec,[1,mp.dm2.Nele,mp.jac.Nmode]); end
    if(any(mp.dm_ind==5)); jacStruct.G5 = jacStruct.G5.*repmat(mp.WspatialVec,[1,mp.dm5.Nele,mp.jac.Nmode]); end
    if(any(mp.dm_ind==8)); jacStruct.G8 = jacStruct.G8.*repmat(mp.WspatialVec,[1,mp.dm8.Nele,mp.jac.Nmode]); end 
    if(any(mp.dm_ind==9)); jacStruct.G9 = jacStruct.G9.*repmat(mp.WspatialVec,[1,mp.dm9.Nele,mp.jac.Nmode]); end
    
    %fprintf('Total Jacobian Calcuation Time: %.2f\n',toc);

    %--Compute the number of total actuators for all DMs used. 
    cvar.NeleAll = mp.dm1.Nele + mp.dm2.Nele + mp.dm3.Nele + mp.dm4.Nele + mp.dm5.Nele + mp.dm6.Nele + mp.dm7.Nele + mp.dm8.Nele + mp.dm9.Nele; %--Number of total actuators used 
    
    %% Wavefront Control

    cvar.Itr = Itr;
    cvar.EfieldVec = EfieldVec;
    cvar.InormHist = InormHist(Itr);
    if mp.flagFiber
        cvar.InormSMFHist = InormSMFHist(Itr);
    end
    [mp,cvar] = falco_ctrl(mp,cvar,jacStruct);
    
    %--Enforce constraints on DM commands 
    % (not needed here--just done here for stats and plotting)
    if(any(mp.dm_ind==1)); mp.dm1 = falco_enforce_dm_constraints(mp.dm1); end
    if(any(mp.dm_ind==2)); mp.dm2 = falco_enforce_dm_constraints(mp.dm2); end
    
    %--Save out regularization used.
    out.log10regHist(Itr) = cvar.log10regUsed; 

%-----------------------------------------------------------------------------------------
%% DM Stats

%--ID and OD of pupil in units of pupil diameters
% ID_pup = 0.303; % for WFIRST, mp.P1.IDnorm
OD_pup = 1.0;

%--Compute the DM surfaces
if(any(mp.dm_ind==1)); DM1surf =  falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  end
if(any(mp.dm_ind==2)); DM2surf =  falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  end

Ndm = length(DM1surf);
dx_dm = mp.P2.compact.dx/mp.P2.D; %--Normalized dx [Units of pupil diameters]
switch mp.centering
    case 'interpixel'
        xs = ( -(Ndm-1)/2:(Ndm-1)/2 )*dx_dm;
    otherwise
        xs = ( -(Ndm/2):(Ndm/2-1) )*dx_dm;
end
[XS,YS] = meshgrid(xs);
RS = sqrt(XS.^2 + YS.^2);
rmsSurf_ele = find(RS>=mp.P1.IDnorm/2 & RS<=OD_pup/2);

%--Compute the RMS stroke
if(any(mp.dm_ind==1))
    Nact = mp.dm1.Nact;
    pitch = mp.dm1.dm_spacing;
    dx_dm = pitch/mp.P2.D; %--Normalized dx [Units of pupil diameters]
    xs = ( -(Nact-1)/2:(Nact-1)/2 )*dx_dm;
    [XS,YS] = meshgrid(xs);
    RS = sqrt(XS.^2 + YS.^2);
    rmsStroke1_ele = find(RS>=mp.P1.IDnorm/2 & RS<=OD_pup/2);
end

%--Calculate and report updated P-V DM voltages.
if(any(mp.dm_ind==1))
    out.dm1.Vpv(Itr) = (max(max(mp.dm1.V))-min(min(mp.dm1.V)));
    Nrail1 = length(find( (mp.dm1.V <= -mp.dm1.maxAbsV) | (mp.dm1.V >= mp.dm1.maxAbsV) ));
    fprintf(' DM1 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm1.Vpv(Itr), Nrail1, mp.dm1.NactTotal, 100*Nrail1/mp.dm1.NactTotal); 
    if(size(mp.dm1.tied,1)>0);  fprintf(' DM1 has %d pairs of tied actuators.\n',size(mp.dm1.tied,1));  end  
end
if(any(mp.dm_ind==2))
    out.dm2.Vpv(Itr) = (max(max(mp.dm2.V))-min(min(mp.dm2.V)));
    Nrail2 = length(find( (mp.dm2.V <= -mp.dm2.maxAbsV) | (mp.dm2.V >= mp.dm2.maxAbsV) ));
    fprintf(' DM2 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm2.Vpv(Itr), Nrail2, mp.dm2.NactTotal, 100*Nrail2/mp.dm2.NactTotal); 
    if(size(mp.dm2.tied,1)>0);  fprintf(' DM2 has %d pairs of tied actuators.\n',size(mp.dm2.tied,1));  end 
end
if(any(mp.dm_ind==8))
    out.dm8.Vpv(Itr) = (max(max(mp.dm8.V))-min(min(mp.dm8.V)));
    Nrail8 = length(find( (mp.dm8.V <= mp.dm8.Vmin) | (mp.dm8.V >= mp.dm8.Vmax) ));
    fprintf(' DM8 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm8.Vpv(Itr), Nrail8,mp.dm8.NactTotal,100*Nrail8/mp.dm8.NactTotal); 
end
if(any(mp.dm_ind==9))
    out.dm9.Vpv(Itr) = (max(max(mp.dm9.V))-min(min(mp.dm9.V)));
    Nrail9 = length(find( (mp.dm9.V <= mp.dm9.Vmin) | (mp.dm9.V >= mp.dm9.Vmax) ));
    fprintf(' DM9 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm9.Vpv(Itr), Nrail9,mp.dm9.NactTotal,100*Nrail9/mp.dm9.NactTotal); 
end

%--Calculate and report updated RMS DM surfaces.
if(any(mp.dm_ind==1))
    out.dm1.Spv(Itr) = max(DM1surf(:))-min(DM1surf(:));
    out.dm1.Srms(Itr) = falco_rms(DM1surf(rmsSurf_ele));
    fprintf('RMS surface of DM1 = %.1f nm\n', 1e9*out.dm1.Srms(Itr))
end
if(any(mp.dm_ind==2))
    out.dm2.Spv(Itr) = max(DM2surf(:))-min(DM2surf(:));
    out.dm2.Srms(Itr) = falco_rms(DM2surf(rmsSurf_ele));
    fprintf('RMS surface of DM2 = %.1f nm\n', 1e9*out.dm2.Srms(Itr))
end

%--Calculate sensitivities to 1nm RMS of Zernike phase aberrations at entrance pupil.
if( isempty(mp.eval.Rsens)==false || isempty(mp.eval.indsZnoll)==false )
    out.Zsens(:,:,Itr) = falco_get_Zernike_sensitivities(mp);
end

% Take the next image to check the contrast level (in simulation only)
tic; fprintf('Getting updated summed image... ');
if mp.flagFiber
    [Im, Ifiber,Ifiber_arr] = falco_get_summed_image(mp);
    if(~mp.flagSim);bench.Femto.averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,Ifiber);end
else
    [Im,~,Imean_arr] = falco_get_summed_image(mp);
end
fprintf('done. Time = %.1f s\n',toc);

%--REPORTING NORMALIZED INTENSITY
InormHist(Itr+1) = mean(Im(mp.Fend.corr.maskBool));
fprintf('Prev and New Measured Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
    InormHist(Itr), InormHist(Itr+1), InormHist(Itr)/InormHist(Itr+1) ); 

fprintf('\n\n');

if mp.flagFiber
    InormSMFHist(Itr+1) = Ifiber;%mean(Im(mp.Fend.corr.maskBool));
    if false%~mp.flagSim
        % Increase the exp time if needed
        peakInCounts = Ifiber*(bench.info.SMFInt0s(1)/mp.peakPSFtint*bench.andor.tint);
        if mp.flagUseCamera4EFCSMF && peakInCounts<2e3
            hcst_andor_setExposureTime(bench,bench.andor.tint+1);
        end
    end
    fprintf('Prev and New Measured SMF Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
        InormSMFHist(Itr), InormSMFHist(Itr+1), InormSMFHist(Itr)/InormSMFHist(Itr+1) ); 

    fprintf('\n\n');
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
    fnWS = sprintf('%sws_%s_Iter%dof%d.mat',mp.path.ws_inprogress,mp.runLabel,Itr,mp.Nitr);
    save(fnWS,'Nitr','Itr','DM1V','DM2V','DM8V','DM9V','InormHist','InormSMFHist','thput_vec','Im')
    fprintf('done.\n\n')
end

% Save cmds for each iteration
if(~mp.flagSim)
    datacmds = hcst_DM_2Dto1D(bench,fliplr(mp.dm1.V'));
    cmds = datacmds+bench.DM.flatvec;
    save([bench.info.outDir,'cdms',datestr(now,'yyyymmddTHHMMSS'),'.mat'],'cmds');
end
%% SAVE THE TRAINING DATA OR RUN THE E-M Algorithm
if(mp.flagTrainModel)
    ev.Itr = Itr;
    if mp.est.flagUseTensorflow
        mp = falco_train_model(mp,ev);
    else
        mp = falco_train_model_matlab(mp,ev);
    end
end
end %--END OF ESTIMATION + CONTROL LOOP
%% ------------------------------------------------------------------------

%% Update progress plot one last time
Itr = Itr + 1;

%--Compute the DM surfaces
if(any(mp.dm_ind==1)); DM1surf =  falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  end
if(any(mp.dm_ind==2)); DM2surf =  falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  end

%--Data to store
if(any(mp.dm_ind==1)); out.dm1.Vall(:,:,Itr) = mp.dm1.V; end
if(any(mp.dm_ind==2)); out.dm2.Vall(:,:,Itr) = mp.dm2.V; end
if(any(mp.dm_ind==8)); out.dm8.Vall(:,Itr) = mp.dm8.V; end
if(any(mp.dm_ind==9)); out.dm9.Vall(:,Itr) = mp.dm9.V; end

%--Calculate the core throughput (at higher resolution to be more accurate)
[mp,thput,ImSimOffaxis] = falco_compute_thput(mp);
if(mp.flagFiber)
    mp.thput_vec(Itr) = max(thput);
else
    mp.thput_vec(Itr) = thput; %--record keeping
end

if(isfield(mp,'testbed'))
    if(mp.flagFiber)
        InormHist_tb.total = InormSMFHist; 
    else
        InormHist_tb.total = InormHist; 
    end
    InormHist_tb.mod(Itr-1) = mean(abs(EfieldVec(:)).^2);
    InormHist_tb.unmod(Itr-1) = mean(IincoVec(:));
    Im_tb.Im = Im;
    Im_tb.E = zeros(size(Im));
    Im_tb.E(mp.Fend.corr.mask) = EfieldVec(:,ceil(mp.Nsbp/2));
    InormHist_tb.beta(Itr-1) = out.log10regHist(Itr-1);
    hProgress = falco_plot_progress_hcst(hProgress,mp,Itr,InormHist_tb,Im_tb,DM1surf,DM2surf);
    
    if isfield(mp,'bench');mp = rmfield(mp,'bench');end

else
    if(mp.flagFiber)
        InormHist_show = InormSMFHist; 
    else
        InormHist_show = InormHist; 
    end
    hProgress = falco_plot_progress(hProgress,mp,Itr,InormHist_show,Im,DM1surf,DM2surf,ImSimOffaxis);
end
% %% Optional output variable: mp
% varargout{1} = mp;

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
out.InormHist = InormHist;
if(mp.flagFiber);out.InormSMFHist = InormSMFHist; end
if(~mp.flagSim); out.InormMod = InormHist_tb.mod;end

fnOut = [mp.path.ws_inprogress filesep mp.runLabel,'_snippet.mat'];

fprintf('\nSaving abridged workspace to file:\n\t%s\n',fnOut)
save(fnOut,'out');
fprintf('...done.\n\n')

%% Save out the data from the workspace
if(mp.flagSaveWS)
    clear cvar G* h* jacStruct; % Save a ton of space when storing the workspace

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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% falco_train_model.m is a nested function in order to save RAM since the 
% output of the Jacobian structure is large and I do not want it copied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% -mp = structure of model parameters
% -cvar = structure containing variables for the controller
% -jacStruct = structure containing the Jacobians
%
% OUTPUTS:
% -mp = structure of model parameters
function mp = falco_train_model(mp,ev)
    Itr = ev.Itr;
    n_batch = mp.NitrTrain;
    % n_batch = 5; %% INITIALIZE THE DATA STRUCTURE THAT SAVES THE TRAINING DATA

    if(Itr==1)
%     if(rem(Itr, n_batch)==1)
        data_train.u1 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u2 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u1p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.u2p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.I = zeros(size(mp.Fend.corr.mask, 1), size(mp.Fend.corr.mask, 2), 2*mp.est.probe.Npairs+1, n_batch);
        save([mp.path.jac, 'jacStruct.mat'], 'jacStruct'); %save jacStruct jacStruct
    else
        data_train = mp.data_train;
    end

    data_train.u1(:, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = mp.dm1.dV;
    if any(mp.dm_ind==2)
        data_train.u2(:, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = mp.dm2.dV;
    else
        data_train.u2(:, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = zeros(size(mp.dm1.dV));
    end
    data_train.u1p(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = ev.Vcube.dm1;
    if any(mp.dm_ind==2)
        data_train.u2p(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = ev.Vcube.dm2;
    else
        data_train.u2p(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = zeros(size(ev.Vcube.dm1));
    end
    data_train.I(:, :, :, Itr - n_batch*floor(Itr/n_batch-1e-3)) = ev.Icube;

    if rem(Itr, n_batch) == 0
        % convert the WFSC data to standard input to the system ID function
        n_activeAct1 = length(mp.dm1.act_ele);
        if any(mp.dm_ind==2)
            n_activeAct2 = length(mp.dm2.act_ele);
        else
            n_activeAct2 = length(mp.dm1.act_ele);
        end
        n_pairs = mp.est.probe.Npairs;
        n_pix = sum(sum(mp.Fend.corr.mask));

        uAll = zeros(n_activeAct1+n_activeAct2, n_batch); % control commands of all the iterations, including both DM1 and DM2
        uProbeAll = zeros(n_activeAct1+n_activeAct2, 2*n_pairs, n_batch); % probe commands of all the iterations, including both DM1 and DM2
        IAll = zeros(n_pix, 2*n_pairs+1, n_batch); % difference image of all the iterations


        for kc = 1 : n_batch % convert the 2D images and DM commands to vectors
            u1_2D = data_train.u1(:, :, kc);
            u2_2D = data_train.u2(:, :, kc);
            uAll(1:n_activeAct1, kc) = u1_2D(mp.dm1.act_ele);
            if any(mp.dm_ind==2)
                uAll(n_activeAct1+1:end, kc) = u2_2D(mp.dm2.act_ele);
            else
                uAll(n_activeAct1+1:end, kc) = u2_2D(mp.dm1.act_ele);
            end
            I_2D = data_train.I(:, :,1, kc);
            IAll(:, 1, kc) = I_2D(mp.Fend.corr.mask);
            for kp = 1 : 2*mp.est.probe.Npairs
                u1p_2D = data_train.u1p(:, :, kp+1, kc) - data_train.u1p(:, :, 1, kc);
                u2p_2D = data_train.u2p(:, :, kp+1, kc) - data_train.u2p(:, :, 1, kc);
                I_2D = data_train.I(:, :,kp+1, kc);

                uProbeAll(1:n_activeAct1, kp, kc) = u1p_2D(mp.dm1.act_ele);
                if any(mp.dm_ind==2)
                    uProbeAll(n_activeAct1+1:end, kp, kc) = u2p_2D(mp.dm2.act_ele);
                else
                    uProbeAll(n_activeAct1+1:end, kp, kc) = u2p_2D(mp.dm1.act_ele);
                end
                IAll(:, kp+1, kc) = I_2D(mp.Fend.corr.mask);
            end
        end

        data_train.u1 = uAll(1:n_activeAct1, :);
        data_train.u2 = uAll(n_activeAct1+1:end, :);
        data_train.u1p = uProbeAll(1:n_activeAct1, :, :);
        data_train.u2p = uProbeAll(n_activeAct1+1:end, :, :);
        data_train.I = IAll;

        save([mp.path.jac, 'data_train.mat'],'data_train') %    save data_train data_train
%         
%         if Itr == n_batch && mp.est.klearning == 1
%             lr = mp.est.lr;
%             lr2 = mp.est.lr2;
%             epoch = mp.est.epoch;
%             Q0 = mp.est.Q0;
%             Q1 = mp.est.Q1;
%             R0 = mp.est.R0;
%             R1 = mp.est.R1;
%             R2 = mp.est.R2;
%         else
%             lr = mp.est.lr;
%             lr2 = mp.est.lr2;
%             epoch = mp.est.epoch;
%             jacLearned = load([mp.path.jac, 'jacStructLearned.mat']);
%             Q0 = jacLearned.noise_coef(1);
%             Q1 = jacLearned.noise_coef(2);
%             R0 = jacLearned.noise_coef(3);
%             R1 = jacLearned.noise_coef(4);
%             R2 = jacLearned.noise_coef(5);
%         end
%         pycommand = ['python falco_systemID_main.py ', mp.path.jac, ' ', ...
%                     num2str(lr), ' ', num2str(lr2), ' ', num2str(epoch), ' ', 'True ', ...
%                     num2str(Q0), ' ', num2str(Q1), ' ', num2str(R0), ' ', num2str(R1), ' ', ...
%                     num2str(R2)];
%         system(pycommand);
        mp.flagUseLearnedJac = true;

        data_train.u1 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u2 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u1p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.u2p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.I = zeros(size(mp.Fend.corr.mask, 1), size(mp.Fend.corr.mask, 2), 2*mp.est.probe.Npairs+1, n_batch);

    end 
    
    mp.data_train = data_train;

end %--END OF FUNCTION


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% falco_train_model_EM.m is a nested function in order to save RAM since the 
% output of the Jacobian structure is large and I do not want it copied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% -mp = structure of model parameters
% -cvar = structure containing variables for the controller
% -jacStruct = structure containing the Jacobians
%
% OUTPUTS:
% -mp = structure of model parameters
function mp = falco_train_model_matlab(mp,ev)
    % a MATLAB version of system ID using E-M algorithm, in contrast the
    % tensorflow based python version
    % 
    Itr = ev.Itr;
    n_batch = mp.NitrTrain;
    % n_batch = 5; %% INITIALIZE THE DATA STRUCTURE THAT SAVES THE TRAINING DATA

    if(Itr==1)
        data_train.u1 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u2 = zeros(mp.dm1.Nact, mp.dm1.Nact, n_batch);
        data_train.u1p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.u2p = zeros(mp.dm1.Nact, mp.dm1.Nact, 2*mp.est.probe.Npairs+1, n_batch);
        data_train.I = zeros(size(mp.Fend.corr.mask, 1), size(mp.Fend.corr.mask, 2), 2*mp.est.probe.Npairs+1, n_batch);
        save([mp.path.jac, 'jacStruct.mat'], 'jacStruct'); %save jacStruct jacStruct
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
%         save([mp.path.jac, 'jacStruct.mat'], 'jacStruct'); %save jacStruct jacStruct
        path2data = mp.path.jac;

        if Itr == n_batch %--Call System ID after final iteration of training
            jacStructLearned = jacStruct;
            jacStructLearned.Q = mp.est.Q;
            jacStructLearned.R = mp.est.R;
            save([mp.path.jac, 'jacStructLearned.mat'],'jacStructLearned')
            LDSlearning(mp)
        else %--All later trainings
            LDSlearning(mp)
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% falco_ctrl_cull.m is a nested function in order to save RAM since the 
% output of the Jacobian structure is large and I do not want it copied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% -mp = structure of model parameters
% -cvar = structure containing variables for the controller
% -jacStruct = structure containing the Jacobians
%
% OUTPUTS:
% -mp = structure of model parameters

function [mp,jacStruct] = falco_ctrl_cull(mp,cvar,jacStruct)

        %% Cull Weak Actuators
        %--MOVE TO A FUNCTION
        %--Reduce the number of actuators used based on their relative strength in the Jacobian
        if(cvar.flagCullAct && cvar.flagRelin)
            fprintf('Weeding out weak actuators from the control Jacobian...\n'); 
            if(any(mp.dm_ind==1))
                %--Crop out very weak-effect actuators
                G1intNorm = zeros(mp.dm1.Nact);
                G1intNorm(1:end) = sum( mean(abs(jacStruct.G1).^2,3), 1);
                G1intNorm = G1intNorm/max(max(G1intNorm));
                mp.dm1.act_ele = find(G1intNorm>=10^(mp.logGmin));
                clear G1intNorm
            end
            if(any(mp.dm_ind==2))
                G2intNorm = zeros(mp.dm2.Nact);
                G2intNorm(1:end) = sum( mean(abs(jacStruct.G2).^2,3),1);
                G2intNorm = G2intNorm/max(max(G2intNorm));
                mp.dm2.act_ele = find(G2intNorm>=10^(mp.logGmin));
                clear G2intNorm
            end
            
            if(any(mp.dm_ind==5))
                G5intNorm = zeros(mp.dm5.Nact);
                G5intNorm(1:end) = sum( mean(abs(jacStruct.G5).^2,3),1);
                G5intNorm = G5intNorm/max(max(G5intNorm));
                mp.dm5.act_ele = find(G5intNorm>=10^(mp.logGmin));
                clear G5intNorm
            end

            if(any(mp.dm_ind==8))
                G8intNorm = zeros(mp.dm8.NactTotal,1);
                G8intNorm(1:end) = sum( mean(abs(jacStruct.G8).^2,3),1);
                G8intNorm = G8intNorm/max(max(G8intNorm));
                mp.dm8.act_ele = find(G8intNorm>=10^(mp.logGmin));
                clear G8intNorm
            end    
            if(any(mp.dm_ind==9))
                G9intNorm = zeros(mp.dm9.NactTotal,1);
                G9intNorm(1:end) = sum( mean(abs(jacStruct.G9).^2,3),1);
                G9intNorm = G9intNorm/max(max(G9intNorm));
                mp.dm9.act_ele = find(G9intNorm>=10^(mp.logGmin));
                clear G9intNorm
            end

            %--Add back in all actuators that are tied (to make the tied actuator logic easier)
            if(any(mp.dm_ind==1))
                for ti=1:size(mp.dm1.tied,1)
                    if(any(mp.dm1.act_ele==mp.dm1.tied(ti,1))==false);  mp.dm1.act_ele = [mp.dm1.act_ele; mp.dm1.tied(ti,1)];  end
                    if(any(mp.dm1.act_ele==mp.dm1.tied(ti,2))==false);  mp.dm1.act_ele = [mp.dm1.act_ele; mp.dm1.tied(ti,2)];  end
                end
                mp.dm1.act_ele = sort(mp.dm1.act_ele); %--Need to sort for the logic in model_Jacobian.m
            end
            if(any(mp.dm_ind==2))
                for ti=1:size(mp.dm2.tied,1)
                    if(any(mp.dm2.act_ele==mp.dm2.tied(ti,1))==false);  mp.dm2.act_ele = [mp.dm2.act_ele; mp.dm2.tied(ti,1)];  end
                    if(any(mp.dm2.act_ele==mp.dm2.tied(ti,2))==false);  mp.dm2.act_ele = [mp.dm2.act_ele; mp.dm2.tied(ti,2)];  end
                end
                mp.dm2.act_ele = sort(mp.dm2.act_ele);
            end
            if(any(mp.dm_ind==8))
                for ti=1:size(mp.dm8.tied,1)
                    if(any(mp.dm8.act_ele==mp.dm8.tied(ti,1))==false);  mp.dm8.act_ele = [mp.dm8.act_ele; mp.dm8.tied(ti,1)];  end
                    if(any(mp.dm8.act_ele==mp.dm8.tied(ti,2))==false);  mp.dm8.act_ele = [mp.dm8.act_ele; mp.dm8.tied(ti,2)];  end
                end
                mp.dm8.act_ele = sort(mp.dm8.act_ele);
            end
            if(any(mp.dm_ind==9))
                for ti=1:size(mp.dm9.tied,1)
                    if(any(mp.dm9.act_ele==mp.dm9.tied(ti,1))==false);  mp.dm9.act_ele = [mp.dm9.act_ele; mp.dm9.tied(ti,1)];  end
                    if(any(mp.dm9.act_ele==mp.dm9.tied(ti,2))==false);  mp.dm9.act_ele = [mp.dm9.act_ele; mp.dm9.tied(ti,2)];  end
                end
                mp.dm9.act_ele = sort(mp.dm9.act_ele);
            end
            
            %--Update the number of elements used per DM
            if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
            if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
            if(any(mp.dm_ind==5)); mp.dm5.Nele = length(mp.dm5.act_ele); end
            if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
            if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end
            %mp.NelePerDMvec = [length(mp.dm1.Nele), length(mp.dm2.Nele), length(mp.dm3.Nele), length(mp.dm4.Nele), length(mp.dm5.Nele), length(mp.dm6.Nele), length(mp.dm7.Nele), length(mp.dm8.Nele), length(mp.dm9.Nele) ];

            if(any(mp.dm_ind==1)); fprintf('  DM1: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm1.Nele, mp.dm1.NactTotal,100*mp.dm1.Nele/mp.dm1.NactTotal); end
            if(any(mp.dm_ind==2)); fprintf('  DM2: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm2.Nele, mp.dm2.NactTotal,100*mp.dm2.Nele/mp.dm2.NactTotal); end
            if(any(mp.dm_ind==5)); fprintf('  DM5: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm5.Nele, mp.dm5.NactTotal,100*mp.dm5.Nele/mp.dm5.NactTotal); end
            if(any(mp.dm_ind==8)); fprintf('  DM8: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm8.Nele, mp.dm8.NactTotal,100*mp.dm8.Nele/mp.dm8.NactTotal); end
            if(any(mp.dm_ind==9)); fprintf('  DM9: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm9.Nele, mp.dm9.NactTotal,100*mp.dm9.Nele/mp.dm9.NactTotal); end
            
            %--Crop out unused actuators from the control Jacobian
            if(any(mp.dm_ind==1)); jacStruct.G1 = jacStruct.G1(:,mp.dm1.act_ele,:); end
            if(any(mp.dm_ind==2)); jacStruct.G2 = jacStruct.G2(:,mp.dm2.act_ele,:); end
            if(any(mp.dm_ind==5)); jacStruct.G5 = jacStruct.G5(:,mp.dm5.act_ele,:); end
            if(any(mp.dm_ind==8)); jacStruct.G8 = jacStruct.G8(:,mp.dm8.act_ele,:); end
            if(any(mp.dm_ind==9)); jacStruct.G9 = jacStruct.G9(:,mp.dm9.act_ele,:); end
        end  

end %--END OF FUNCTION falco_ctrl_cull.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I placed the function falco_ctrl.m here as a nested
% function in order to save RAM since the Jacobian structure
% is large and I do not want it copied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% -mp = structure of model parameters
% -cvar = structure containing variables for the controller
% -jacStruct = structure containing the Jacobians
%
% OUTPUTS:
% -mp = structure of model parameters
% REVISION HISTORY:
% --------------
% Created by A.J. Riggs on 2018-10-04 by extracting material from falco_wfsc_loop.m.
% ---------------

function [mp,cvar] = falco_ctrl(mp,cvar,jacStruct)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Control Algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    fprintf('Using the Jacobian to make other matrices...'); tic;

    %--Compute matrices for linear control with regular EFC
    cvar.GstarG_wsum  = zeros(cvar.NeleAll,cvar.NeleAll); 
    cvar.RealGstarEab_wsum = zeros(cvar.NeleAll, 1);

    for im=1:mp.jac.Nmode

        Gstack = [jacStruct.G1(:,:,im), jacStruct.G2(:,:,im), jacStruct.G3(:,:,im), jacStruct.G4(:,:,im), jacStruct.G5(:,:,im), jacStruct.G6(:,:,im), jacStruct.G7(:,:,im), jacStruct.G8(:,:,im), jacStruct.G9(:,:,im)];

        %--Square matrix part stays the same if no re-linearization has occurrred. 
        cvar.GstarG_wsum  = cvar.GstarG_wsum  + mp.jac.weights(im)*real(Gstack'*Gstack); 

        %--The G^*E part changes each iteration because the E-field changes.
        Eweighted = mp.WspatialVec.*cvar.EfieldVec(:,im); %--Apply 2-D spatial weighting to E-field in dark hole pixels.
        cvar.RealGstarEab_wsum = cvar.RealGstarEab_wsum + mp.jac.weights(im)*real(Gstack'*Eweighted); %--Apply the Jacobian weights and add to the total.

    end
    clear GallCell Gstack Eweighted % save RAM
    
    %--Make the regularization matrix. (Define only the diagonal here to save RAM.)
    cvar.EyeGstarGdiag = max(diag(cvar.GstarG_wsum ))*ones(cvar.NeleAll,1);
    fprintf(' done. Time: %.3f\n',toc);

    %--Call the Controller Function
    fprintf('Control beginning ...\n'); tic
    switch lower(mp.controller)

        case{'plannedefc'} %--EFC regularization is scheduled ahead of time
            [dDM,cvar] = falco_ctrl_planned_EFC(mp,cvar);

        case{'gridsearchefc'}  %--Empirical grid search of EFC. Scaling factor for DM commands too.
            [dDM,cvar] = falco_ctrl_grid_search_EFC(mp,cvar);
            
            
        case{'plannedefccon'} %--Constrained-EFC regularization is scheduled ahead of time
            [dDM,cvar] = falco_ctrl_planned_EFCcon(mp,cvar);
            
        case{'sm-cvx'} %--Constrained & bounded stroke minimization using CVX. The quadratic cost function is solved directly CVX.
            cvar.dummy = 1;
            [dDM,cvar] = falco_ctrl_SM_CVX(mp,cvar);
            
        case{'tsm'}
            cvar.dummy = 1;
            [dDM,cvar] = falco_ctrl_total_stroke_minimization(mp,cvar); 
            
    end
    fprintf(' done. Time: %.3f sec\n',toc);
    
    %% Updates to DM commands

    %--Update the DM commands by adding the delta control signal
    if(any(mp.dm_ind==1));  mp.dm1.V = mp.dm1.V + dDM.dDM1V;  end
    if(any(mp.dm_ind==2));  mp.dm2.V = mp.dm2.V + dDM.dDM2V;  end
    if(any(mp.dm_ind==3));  mp.dm3.V = mp.dm3.V + dDM.dDM3V;  end
    if(any(mp.dm_ind==4));  mp.dm4.V = mp.dm4.V + dDM.dDM4V;  end
    if(any(mp.dm_ind==5));  mp.dm5.V = mp.dm5.V + dDM.dDM5V;  end
    if(any(mp.dm_ind==6));  mp.dm6.V = mp.dm6.V + dDM.dDM6V;  end
    if(any(mp.dm_ind==7));  mp.dm7.V = mp.dm7.V + dDM.dDM7V;  end
    if(any(mp.dm_ind==8));  mp.dm8.V = mp.dm8.V + dDM.dDM8V;  end
    if(any(mp.dm_ind==9));  mp.dm9.V = mp.dm9.V + dDM.dDM9V;  end

    %%--Save the delta from the previous command
    if(any(mp.dm_ind==1));  mp.dm1.dV = dDM.dDM1V;  end
    if(any(mp.dm_ind==2));  mp.dm2.dV = dDM.dDM2V;  end
    if(any(mp.dm_ind==3));  mp.dm3.dV = dDM.dDM3V;  end
    if(any(mp.dm_ind==4));  mp.dm4.dV = dDM.dDM4V;  end
    if(any(mp.dm_ind==5));  mp.dm5.dV = dDM.dDM5V;  end
    if(any(mp.dm_ind==6));  mp.dm6.dV = dDM.dDM6V;  end
    if(any(mp.dm_ind==7));  mp.dm7.dV = dDM.dDM7V;  end
    if(any(mp.dm_ind==8));  mp.dm8.dV = dDM.dDM8V;  end
    if(any(mp.dm_ind==9));  mp.dm9.dV = dDM.dDM9V;  end
    
%     %--Update the tied actuator pairs
%     if(any(mp.dm_ind==1));  mp.dm1.tied = dDM.dm1tied;  end
%     if(any(mp.dm_ind==2));  mp.dm2.tied = dDM.dm2tied;  end

end %--END OF NESTED FUNCTION

end %--END OF main FUNCTION
