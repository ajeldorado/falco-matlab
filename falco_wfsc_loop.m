% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% High-level function to run design simulations of various types of Lyot coronagraphs.
%
% Data are mostly passed in structures.
%
% Modified again by A.J. Riggs on May 23, 2017 to eliminate a lot of unnecessary
%   variables for full-knowledge design work in the main script,
%   in the config function, in the EFC controller, and in
%   the image generating function.
% Modified by A.J. Riggs on May 10, 2017 to eliminate a lot of unnecessary
%   variables for full-knowledge design work in the main script, falco_3DM_main.m .
% Modified by A.J. Riggs in April 2016 to include Babinet's principle for a 3DMLC
%   propagation option.
% Modified by A.J. Riggs from A.J.'s general coronagraphic WFSC code for design only 
%   in April 2016.
% Adapted by A.J. Riggs from A.J.'s Princeton HCIL code on August 31, 2016.

function [out] = falco_wfsc_loop(mp)


% function [out] = falco_wfsc_loop(fn_config,varargin)
% % Set default values of input parameters
% mp.flagPlot = false; % flag to plot PSF correction in real time
% %--Enable different arguments values by using varargin
% icav = 0;                     % index in cell array varargin
% while icav < size(varargin, 2)
%     icav = icav + 1;
%     switch lower(varargin{icav})
%       case {'plot','plotting'}
%         mp.flagPlot = true;  %--Value to use for turning on plots
%       otherwise
%         error('falco_wfsc_loop: Unknown keyword: %s\n', varargin{icav});
%           
%     end
% end

%% Save the config file    
mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.F4.corr.Rin),'_OWA',num2str(mp.F4.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

mp.folders.dummy = 1; %--Initialize the folders structure in case it doesn't already exist
if(isfield(mp.folders,'brief')==false)
    mainPath = pwd;
    mp.folders.m = mainPath;
    mp.folders.brief = [mainPath filesep 'data' filesep 'brief' filesep];      % Store minimal data to re-construct the data from the run: the config files and "out" structure go here
end

fn_config = [mp.runLabel,'_config.mat'];

cd(mp.folders.brief)
save(fn_config)
cd(mp.folders.m)
fprintf('Saved the config file: \t%s\n',fn_config)


%% Get configuration data from a function file
[mp,out] = falco_init_ws(fn_config, mp.folders.brief);

%% Jacobian storage
% G_mat_fname = sprintf('G_%s_%dDM_%dx%dx%dact_%dpix_%dpctBW_at%dnm_Nsbp%02d.mat',...  %--Name of the Jacobian file if it is saved
%     mp.coro, mp.num_dms,mp.dm1.Nact,...
%     mp.dm2.Nact,mp.dm9.NactTotal,length(mp.cor_ele),round(mp.fracBW*100), round(mp.lambda0*1e9), mp.Nsbp);

% cd(mp.folders.jac)
% if( exist(G_mat_fname, 'file') ~= 2) %--if 2, then the file by that name exists
%     flagCalcJac = true; %--Calculate the starting Jacobian again if the file does not exist
% end
% cd(mp.folders.m)

%% Initializations of Arrays for Data Storage 

%--Raw contrast (broadband)
InormHist = zeros(mp.Nitr,1); % Measured, mean raw contrast in scoring regino of dark hole.

% ImHist = single( zeros(mp.F4.Neta,mp.F4.Nxi,mp.Nitr+1) ); %--Full PSF after each correction step
%
% %--Store the DM surfaces (REQUIRES LOTS OF STORAGE)
% DM1S_array = single(zeros(mp.dm1.compact.Ndm,mp.dm1.compact.Ndm,mp.Nitr+1));
% DM2S_array = single(zeros(mp.dm2.compact.Ndm,mp.dm2.compact.Ndm,mp.Nitr+1));


%% Take initial broadband images
%EfieldCorrTrue = zeros(mp.F4.corr.Npix,mp.jac.Nmode,mp.Nitr+1); % (Simulation only) Vectorized true starlight E-field at each pixel and wavelength

if(mp.flagPlot); figure(101); imagesc(mp.P1.full.mask);axis image; colorbar; title('pupil');drawnow; end

if(mp.flagPlot && (length(mp.P4.full.mask)==length(mp.P1.full.mask))); figure(102); imagesc(mp.P4.full.mask);axis image; colorbar; title('Lyot stop');drawnow; end
if(mp.flagPlot && isfield(mp,'P3.full.mask')); figure(103); imagesc(padOrCropEven(mp.P1.full.mask,mp.P3.full.Narr).*mp.P3.full.mask);axis image; colorbar; drawnow; end

%% Take initial broadband image 

Im = falco_get_summed_image(mp);
% ImHist(:,:,1) = falco_get_summed_image(mp);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin the Correction Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Itr=1:mp.Nitr

    %--Start of new estimation+control iteration
    fprintf(['Iteration: ' num2str(Itr) '/' num2str(mp.Nitr) '\n' ]);

    %--Re-compute the starlight normalization factor for the compact and full models (to convert images to normalized intensity). No tip/tilt necessary.
    mp = falco_get_PSF_norm_factor(mp);
    
    %% Updated DM data
    %--Change the selected DMs if using the scheduled EFC controller
    switch mp.controller
        case{'plannedEFC'} 
            mp.dm_ind = mp.dm_ind_sched{Itr};
    end
    %--Report which DMs are used in this iteration
    fprintf('DMs to be used in this iteration = ['); for jj = 1:length(mp.dm_ind); fprintf(' %d',mp.dm_ind(jj)); end; fprintf(' ]\n');

    %--Fill in History of DM commands to Store
    if(isfield(mp,'dm1')); if(isfield(mp.dm1,'V'));  out.dm1.Vall(:,:,Itr) = mp.dm1.V;  end;  end
    if(isfield(mp,'dm2')); if(isfield(mp.dm2,'V'));  out.dm2.Vall(:,:,Itr) = mp.dm2.V;  end;  end
    if(isfield(mp,'dm8')); if(isfield(mp.dm8,'V'));  out.dm8.Vall(:,Itr) = mp.dm8.V(:);  end;  end
    if(isfield(mp,'dm9')); if(isfield(mp.dm9,'V'));  out.dm9.Vall(:,Itr) = mp.dm9.V(:);  end;  end
%     if(any(mp.dm_ind==1)); out.dm1.Vall(:,:,Itr) = mp.dm1.V; end
%     if(any(mp.dm_ind==2)); out.dm2.Vall(:,:,Itr) = mp.dm2.V; end 
%     if(any(mp.dm_ind==8)); out.dm8.Vall(:,Itr) = mp.dm8.V(:); end
%     if(any(mp.dm_ind==9)); out.dm9.Vall(:,Itr) = mp.dm9.V; end

    %--Compute the DM surfaces
    if(any(mp.dm_ind==1)); DM1surf =  falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  end
    if(any(mp.dm_ind==2)); DM2surf =  falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  end
    % if(any(mp.dm_ind==9)); DM9phase =  padOrCropEven(falco_dm_surf_from_cube(mp.dm9,mp.dm9.compact), mp.dm9.compact.NxiFPM); end
    switch mp.coro
        case{'EHLC'}
            mp.DM8surf = falco_gen_EHLC_FPM_surf_from_cube(mp.dm8,'compact'); %--Metal layer profile [m]
            mp.DM9surf = falco_gen_EHLC_FPM_surf_from_cube(mp.dm9,'compact'); %--Dielectric layer profile [m]
        case{'HLC','APHLC','SPHLC'}
            mp.DM8surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm8,'compact'); %--Metal layer profile [m]
            mp.DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9,'compact'); %--Dielectric layer profile [m]
    end

    %% Updated plot and reporting
    %--Calculate the core throughput (at higher resolution to be more accurate)
    [mp,thput] = falco_compute_thput(mp);
    mp.thput_vec(Itr) = thput;

    %--Compute the current contrast level
    InormHist(Itr) = mean(Im(mp.F4.corr.maskBool));

    %--Plot the updates to the DMs and PSF
    if(Itr==1); hProgress.master = 1; end %--dummy value to intialize the handle variable
    hProgress = falco_plot_progress(hProgress,mp,Itr,InormHist,Im,DM1surf,DM2surf);

    %% Updated selection of Zernike modes targeted by the controller
    %--Decide with Zernike modes to include in the Jacobian
    if(Itr==1)
        mp.jac.zerns0 = mp.jac.zerns;
       %mp.jac.maxZnoll0 = mp.jac.maxZnoll; 
    end
    if( (any(mp.dm_ind==8)==false) && (any(mp.dm_ind==9)==false) )
        mp.jac.zerns = 1;
    else
        mp.jac.zerns = mp.jac.zerns0;
    end
    fprintf('Zernike modes used in this Jacobian:\t'); fprintf('%d ',mp.jac.zerns); fprintf('\n');
    %--Re-compute the Jacobian weights
    mp = falco_config_jac_weights(mp); 
    


    %% Wavefront Estimation
    %----------------------      Wavefront Estimation    ----------------------
    switch lower(mp.estimator)
        case{'perfect'}
            EfieldVec  = falco_est_perfect_Efield_with_Zernikes(mp);
        case{'batch'}

    end

    %-----------------------------------------------------------------------------------------
    %% Wavefront Control

    cvar.Itr = Itr;
    cvar.EfieldVec = EfieldVec;
    cvar.InormHist = InormHist(Itr);
    
    %--Relinearize about the DMs only at the iteration numbers in mp.relinItrVec.
    if(any(mp.relinItrVec==Itr))
        cvar.flagRelin=true;
    else
        cvar.flagRelin=false;
    end
    
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
        %--Re-include all actuators in the basis set.
        if(any(mp.dm_ind==1)); mp.dm1.act_ele = 1:mp.dm1.NactTotal; end
        if(any(mp.dm_ind==2)); mp.dm2.act_ele = 1:mp.dm2.NactTotal; end
        if(any(mp.dm_ind==8)); mp.dm8.act_ele = 1:mp.dm8.NactTotal; end
        if(any(mp.dm_ind==9)); mp.dm9.act_ele = 1:mp.dm9.NactTotal; end
        %--Update the number of elements used per DM
        if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
        if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
        if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
        if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end
    end

    %% Compute the control Jacobians for each DM
    if( (Itr==1) || cvar.flagRelin )
        jacStruct =  model_Jacobian(mp); %--Get structure containing Jacobians
    end
    
    % %--Save or load a previous Jacobian (esp. useful for testbeds)
    %     if(Itr==1)
    %         cd(mp.folders.jac)
    % %             save(G_mat_fname,'G1','G2','G9','-v7.3');
    % %             save(G_mat_fname,'G9','-v7.3');
    %         cd(mp.folders.m)
    %     end
    % elseif( Itr==1 && flagCalcJac==0 )    
    %     cd(mp.folders.jac)
    %         %load(G_mat_fname)
    %     cd(mp.folders.m)    
    
    %% Cull actuators, but only if(cvar.flagCullAct && cvar.flagRelin)
    [mp,jacStruct] = falco_ctrl_cull(mp,cvar,jacStruct);

    % Add spatially-dependent weighting to the control Jacobians
    if(any(mp.dm_ind==1)); jacStruct.G1 = jacStruct.G1.*repmat(mp.WspatialVec,[1,mp.dm1.Nele,mp.jac.Nmode]); end
    if(any(mp.dm_ind==2)); jacStruct.G2 = jacStruct.G2.*repmat(mp.WspatialVec,[1,mp.dm2.Nele,mp.jac.Nmode]); end  
    if(any(mp.dm_ind==8)); jacStruct.G8 = jacStruct.G8.*repmat(mp.WspatialVec,[1,mp.dm8.Nele,mp.jac.Nmode]); end 
    if(any(mp.dm_ind==9)); jacStruct.G9 = jacStruct.G9.*repmat(mp.WspatialVec,[1,mp.dm9.Nele,mp.jac.Nmode]); end 

    %fprintf('Total Jacobian Calcuation Time: %.2f\n',toc);

    %--Compute the number of total actuators for all DMs used. 
    cvar.NeleAll = mp.dm1.Nele + mp.dm2.Nele + mp.dm3.Nele + mp.dm4.Nele + mp.dm5.Nele + mp.dm6.Nele + mp.dm7.Nele + mp.dm8.Nele + mp.dm9.Nele; %--Number of total actuators used 

    %% Zero out Jacobian response for railed actuators
    
%     if(any(mp.dm_ind==1))
%         dm1_act_ele_railed = find( (mp.dm1.V(mp.dm1.act_ele) < -mp.dm1.maxAbsV) | (mp.dm1.V(mp.dm1.act_ele) > mp.dm1.maxAbsV) ); 
%         jacStruct.G1(:,dm1_act_ele_railed,:) = 0*jacStruct.G1(:,dm1_act_ele_railed,:);
%     end
%     if(any(mp.dm_ind==2))
%         dm2_act_ele_railed = find( (mp.dm2.V(mp.dm2.act_ele) < -mp.dm2.maxAbsV) | (mp.dm2.V(mp.dm2.act_ele) > mp.dm2.maxAbsV) ); 
%         jacStruct.G2(:,dm2_act_ele_railed,:) = 0*jacStruct.G2(:,dm2_act_ele_railed,:);
%     end
%     if(any(mp.dm_ind==8))
%         dm8_act_ele_railed = find( (mp.dm8.V(mp.dm8.act_ele) < mp.dm8.Vmin) | (mp.dm8.V(mp.dm8.act_ele) > mp.dm8.Vmax) );
%         jacStruct.G8(:,dm8_act_ele_railed,:) = 0*jacStruct.G8(:,dm8_act_ele_railed,:);
%     end    
%     if(any(mp.dm_ind==9))
%         dm9_act_ele_railed = find( (mp.dm9.V(mp.dm9.act_ele) < mp.dm9.Vmin) | (mp.dm9.V(mp.dm9.act_ele) > mp.dm9.Vmax) );
%         jacStruct.G9(:,dm9_act_ele_railed,:) = 0*jacStruct.G9(:,dm9_act_ele_railed,:);
%     end
 
    %% Remove Jacobian response for railed actuators
%     
%     %dm1_subset=[]; dm2_subset=[]; dm3_subset=[]; dm4_subset=[]; dm5_subset=[]; dm6_subset=[]; dm7_subset=[]; dm8_subset=[]; dm9_subset=[]; 
%     if(any(mp.dm_ind==1))
%         dm1_subset = find( (mp.dm1.V(mp.dm1.act_ele) > -mp.dm1.maxAbsV) & (mp.dm1.V(mp.dm1.act_ele) < mp.dm1.maxAbsV) ); %--Relative, unrailed-actuator indices of the original subset
%         jacStruct.G1 = jacStruct.G1(:,dm1_subset,:);
%         mp.dm1.act_ele = intersect( mp.dm1.act_ele, find( (mp.dm1.V > -mp.dm1.maxAbsV) & (mp.dm1.V < mp.dm1.maxAbsV) ) ); %--Absolute indices out of the entire basis set
%         mp.dm1.Nele = length(mp.dm1.act_ele);
%     end
%     
%     if(any(mp.dm_ind==2))
%         dm2_subset = find( (mp.dm2.V(mp.dm2.act_ele) > -mp.dm2.maxAbsV) & (mp.dm2.V(mp.dm2.act_ele) < mp.dm2.maxAbsV) ); %--Relative, unrailed-actuator indices of the original subset
%         jacStruct.G2 = jacStruct.G2(:,dm2_subset,:);
%         mp.dm2.act_ele = intersect( mp.dm2.act_ele, find( (mp.dm2.V > -mp.dm2.maxAbsV) & (mp.dm2.V < mp.dm2.maxAbsV) ) ); %--Absolute indices out of the entire basis set
%         mp.dm2.Nele = length(mp.dm2.act_ele);
%     end
%     if(any(mp.dm_ind==8))
%         dm8_subset = find( (mp.dm8.V(mp.dm8.act_ele) > mp.dm8.Vmin) & (mp.dm8.V(mp.dm8.act_ele) < mp.dm8.Vmax) ); %--Relative, unrailed-actuator indices of the original subset
%         jacStruct.G8 = jacStruct.G8(:,dm8_subset,:);
%         mp.dm8.act_ele = intersect( mp.dm8.act_ele, find( (mp.dm8.V > mp.dm8.Vmin) & (mp.dm8.V < mp.dm8.Vmax) ) ); %--Absolute indices out of the entire basis set
%         mp.dm8.Nele = length(mp.dm8.act_ele);
%     end    
%     if(any(mp.dm_ind==9))
%         dm9_subset = find( (mp.dm9.V(mp.dm9.act_ele) > mp.dm9.Vmin) & (mp.dm9.V(mp.dm9.act_ele) < mp.dm9.Vmax) ); %--Relative, unrailed-actuator indices of the original subset
%         jacStruct.G9 = jacStruct.G9(:,dm9_subset,:);
%         mp.dm9.act_ele = intersect( mp.dm9.act_ele, find( (mp.dm9.V > mp.dm9.Vmin) & (mp.dm9.V < mp.dm9.Vmax) ) ); %--Absolute indices out of the entire basis set
%         mp.dm9.Nele = length(mp.dm9.act_ele);
%     end
%     %--Compute the number of total actuators for all DMs used. 
%     cvar.NeleAll = mp.dm1.Nele + mp.dm2.Nele + mp.dm3.Nele + mp.dm4.Nele + mp.dm5.Nele + mp.dm6.Nele + mp.dm7.Nele + mp.dm8.Nele + mp.dm9.Nele; %--Number of total actuators used 
    
% %     %--Compute indices in the total Jacobian to keep
% %     cvar.dm_subset = [...
% %         dm1_subset; ...
% %         mp.dm1.Nele+dm2_subset; ...
% %         mp.dm1.Nele+mp.dm2.Nele + dm8_subset; ...
% %         mp.dm1.Nele+mp.dm2.Nele+mp.dm8.Nele + dm9_subset ];
    
    %% Wavefront Control
    [mp,cvar] = falco_ctrl(mp,cvar,jacStruct);
    %--Save out regularization used.
    out.log10regHist(Itr) = cvar.log10regUsed; 
    

    

%     switch mp.controller
%         case{'plannedEFC','gridsearchEFC'}
%             %--Remove railed actuators from the basis set
%             if(any(mp.dm_ind==1)); mp.dm1.act_ele = intersect( mp.dm1.act_ele, find( (mp.dm1.V > -mp.dm1.maxAbsV) & (mp.dm1.V < mp.dm1.maxAbsV) ) ); end
%             if(any(mp.dm_ind==2)); mp.dm2.act_ele = intersect( mp.dm2.act_ele, find( (mp.dm2.V > -mp.dm2.maxAbsV) & (mp.dm2.V < mp.dm2.maxAbsV) ) ); end
%             if(any(mp.dm_ind==8)); mp.dm8.act_ele = intersect( mp.dm8.act_ele, find( (mp.dm8.V > mp.dm8.Vmin) & (mp.dm8.V < mp.dm8.Vmax) ) ); end   
%             if(any(mp.dm_ind==9)); mp.dm9.act_ele = intersect( mp.dm9.act_ele, find( (mp.dm9.V > mp.dm9.Vmin) & (mp.dm9.V < mp.dm9.Vmax) ) ); end   
%     end
%     %--Update the number of elements used per DM
%     if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
%     if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
%     if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
%     if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end
%-----------------------------------------------------------------------------------------
%% DM Stats

%--ID and OD of pupil in units of pupil diameters
% ID_pup = 0.303; % mp.P1.IDnorm
OD_pup = 1.0;

% DM1surf = DM1S_array(:,:,end);
% DM2surf = DM2S_array(:,:,end);
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
rms_ele = find(RS>=mp.P1.IDnorm/2 & RS<=OD_pup/2);

% out.dm1.Spv = max(DM1surf(:))-min(DM1surf(:));
% out.dm2.Spv = max(DM2surf(:))-min(DM2surf(:));
% fprintf('P-V surface of DM1 = %.1f nm\n', 1e9*out.dm1.Spv)
% fprintf('P-V surface of DM2 = %.1f nm\n', 1e9*out.dm2.Spv)
% 
% pvfsDM1 = max(mp.dm1.V(:))-min(mp.dm1.V(:));
% pvfsDM2 = max(mp.dm2.V(:))-min(mp.dm2.V(:));
% fprintf('P-V free stroke of DM1 = %.1f nm\n', pvfsDM1)
% fprintf('P-V free stroke of DM2 = %.1f nm\n', pvfsDM2)

%--Calculate and report updated P-V DM voltages.
if(any(mp.dm_ind==1))
    out.dm1.Vpv(Itr) = (max(max(mp.dm1.V))-min(min(mp.dm1.V)));
    Nrail1 = length(find( (mp.dm1.V <= -mp.dm1.maxAbsV) | (mp.dm1.V >= mp.dm1.maxAbsV) ));
    fprintf(' DM1 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm1.Vpv(Itr), Nrail1, mp.dm1.NactTotal, 100*Nrail1/mp.dm1.NactTotal); 
end
if(any(mp.dm_ind==2))
    out.dm2.Vpv(Itr) = (max(max(mp.dm2.V))-min(min(mp.dm2.V)));
    Nrail2 = length(find( (mp.dm2.V <= -mp.dm2.maxAbsV) | (mp.dm2.V >= mp.dm2.maxAbsV) ));
    fprintf(' DM2 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm2.Vpv(Itr), Nrail2, mp.dm2.NactTotal, 100*Nrail2/mp.dm2.NactTotal); 
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
    out.dm1.Spv = max(DM1surf(:))-min(DM1surf(:));
    out.dm1.Srms = rms(DM1surf(rms_ele));
    fprintf('RMS surface of DM1 = %.1f nm\n', 1e9*out.dm1.Srms)
end
if(any(mp.dm_ind==2))
    out.dm2.Spv = max(DM2surf(:))-min(DM2surf(:));
    out.dm2.Srms = rms(DM2surf(rms_ele));
    fprintf('RMS surface of DM2 = %.1f nm\n', 1e9*out.dm2.Srms)
end


% Take the next image to check the contrast level (in simulation only)
Im = falco_get_summed_image(mp);
% ImHist(:,:,Itr+1) = falco_get_summed_image(mp);

%--REPORTING NORMALIZED INTENSITY
if( (Itr==mp.Nitr) || (strcmpi(mp.controller,'conEFC') && (numel(mp.ctrl.muVec)==1)  ) ) 
    InormHist(Itr+1) = mean(Im(mp.F4.corr.maskBool));
%     ImBandAvg_current = ImHist(:,:,Itr+1);
%     InormHist(Itr+1) = mean(ImBandAvg_current(mp.F4.corr.maskBool));

    fprintf('Prev and New Measured Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
        InormHist(Itr), InormHist(Itr+1), InormHist(Itr)/InormHist(Itr+1) ); 
else
    fprintf('Prev and New Measured Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
        InormHist(Itr), cvar.cMin, InormHist(Itr)/cvar.cMin ); 
end
fprintf('\n\n');


% %--Save out DM commands after each iteration in case the run crashes part way through.
% fprintf('Saving DM commands for this iteration...')
% cd(mp.folders.ws_inprogress)
% 
%         if(any(mp.dm_ind==1)); DM1V = mp.dm1.V; else; DM1V = 0; end
%         if(any(mp.dm_ind==2)); DM2V = mp.dm2.V; else; DM2V = 0; end
%         if(any(mp.dm_ind==3)); DM3V = mp.dm3.V; else; DM3V = 0; end
%         if(any(mp.dm_ind==8)); DM8V = mp.dm8.V; else; DM8V = 0; end
%         if(any(mp.dm_ind==9)); DM9V = mp.dm9.V; else; DM9V = 0; end
%         Nitr = mp.Nitr;
%         thput_vec = mp.thput_vec;
% 
%         fnWS = sprintf('ws_%s_Iter%dof%d.mat',mp.runLabel,Itr,mp.Nitr);
%         save(fnWS,'Nitr','Itr','DM1V','DM2V','DM3V','DM8V','DM9V','InormHist','thput_vec')
% cd(mp.folders.m)
% fprintf('done.\n\n')


end %--END OF ESTIMATION + CONTROL LOOP
%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------


%% Update plot one last time
Itr = Itr + 1;

%--Compute the DM surfaces
if(any(mp.dm_ind==1)); DM1surf =  falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);  end
if(any(mp.dm_ind==2)); DM2surf =  falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);  end
%if(any(mp.dm_ind==9)); DM9phase =  padOrCropEven(falco_dm_surf_from_cube(mp.dm9,mp.dm9.compact), mp.dm9.compact.NxiFPM); end

%--Data to store
if(any(mp.dm_ind==1)); out.dm1.Vall(:,:,Itr) = mp.dm1.V; end % DM1S_array(:,:,Itr) = single(DM1surf); end
if(any(mp.dm_ind==2)); out.dm2.Vall(:,:,Itr) = mp.dm2.V; end %DM2S_array(:,:,Itr) = single(DM2surf); end
if(any(mp.dm_ind==8)); out.dm8.Vall(:,Itr) = mp.dm8.V;  end
if(any(mp.dm_ind==9)); out.dm9.Vall(:,Itr) = mp.dm9.V;  end

% modvar.x_offset = mp.thput_eval_x;
% modvar.y_offset = mp.thput_eval_y;
% modvar.sbpIndex = mp.si_ref; 
% modvar.whichSource = 'offaxis';

%--Calculate the core throughput (at higher resolution to be more accurate)
[mp,thput] = falco_compute_thput(mp);
mp.thput_vec(Itr) = thput;

hProgress = falco_plot_progress(hProgress,mp,Itr,InormHist,Im,DM1surf,DM2surf);

%% Save the final DM commands separately for faster reference
if(isfield(mp,'dm1')); if(isfield(mp.dm1,'V')); out.DM1V = mp.dm1.V; end; end
if(isfield(mp,'dm2')); if(isfield(mp.dm2,'V')); out.DM2V = mp.dm2.V; end; end
switch mp.coro
    case{'HLC','EHLC'}
        if(isfield(mp.dm8,'V')); out.DM8V = mp.dm8.V;  end
        if(isfield(mp.dm9,'V')); out.DM9V = mp.dm9.V;  end
end
% if(isfield(mp,'dm8')); if(isfield(mp.dm8,'V')); out.DM8V = mp.dm8.V; end; end
% if(isfield(mp,'dm8')); if(isfield(mp.dm9,'V')); out.DM9V = mp.dm9.V; end; end


%% Save out an abridged workspace

%--Variables to save out:
% contrast vs iter
% regularization history
%  DM1surf,DM1V, DM2surf,DM2V, DM8surf,DM9surf, fpm sampling, base pmgi thickness, base nickel thickness, dm_tilts, aoi, ...
% to reproduce your basic design results of NI, throughput, etc..

out.thput = mp.thput_vec;
out.Nitr = mp.Nitr;
out.InormHist = InormHist;

fnOut = [mp.runLabel,'_snippet.mat'];

disp(['\nSaving abridged workspace to file ' fnOut '...'])
cd(mp.folders.brief)
save(fnOut,'out');
fprintf('done.\n\n')
cd(mp.folders.m)


%% Save out the data from the workspace
clear cvar G* h* jacStruct; % Save a ton of space when storing the workspace

% Don't bother saving the large 2-D, floating point maps in the workspace (they take up too much space)
mp.P1.full.mask=1; mp.P1.compact.mask=1;
mp.P3.full.mask=1; mp.P3.compact.mask=1;
mp.P4.full.mask=1; mp.P4.compact.mask=1;
mp.F3.full.mask=1; mp.F3.compact.mask=1;

mp.P1.full.E = 1; mp.P1.compact.E=1; mp.Eplanet=1; 
mp.dm1.full.mask = 1; mp.dm1.compact.mask = 1; mp.dm2.full.mask = 1; mp.dm2.compact.mask = 1;
mp.complexTransFull = 1; mp.complexTransCompact = 1;

mp.dm1.compact.inf_datacube = 0;
mp.dm2.compact.inf_datacube = 0;
mp.dm8.compact.inf_datacube = 0;
mp.dm9.compact.inf_datacube = 0;
mp.dm8.inf_datacube = 0;
mp.dm9.inf_datacube = 0;

fnAll = [mp.runLabel,'_all.mat'];
% fnWS = ['ws_',mp.runLabel,'_',num2str(mp.Nitr),'its.mat'];

disp(['Saving workspace to file ' fnAll '...'])
cd(mp.folders.ws)
save(fnAll);
fprintf('done.\n\n')
cd(mp.folders.m)


end %--END OF main FUNCTION
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I placed the function falco_ctrl_cull.m here as a nested
% function in order to save RAM since the output of the Jacobian structure
% is large and I do not want it copied.
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
                %if(mp.flagPlot); figure(81); imagesc(log10(G1intNorm),[-6 0]); axis xy equal tight; colorbar; end
                clear G1intNorm
            end
                if(any(mp.dm_ind==2))
                G2intNorm = zeros(mp.dm2.Nact);
                G2intNorm(1:end) = sum( mean(abs(jacStruct.G2).^2,3),1);
                G2intNorm = G2intNorm/max(max(G2intNorm));
                mp.dm2.act_ele = find(G2intNorm>=10^(mp.logGmin));
                %if(mp.flagPlot); figure(82); imagesc(log10(G2intNorm),[-6 0]); axis xy equal tight; colorbar; end
                clear G2intNorm
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

            %--Update the number of elements used per DM
            if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
            if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
            if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
            if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end
            %mp.NelePerDMvec = [length(mp.dm1.Nele), length(mp.dm2.Nele), length(mp.dm3.Nele), length(mp.dm4.Nele), length(mp.dm5.Nele), length(mp.dm6.Nele), length(mp.dm7.Nele), length(mp.dm8.Nele), length(mp.dm9.Nele) ];

            if(any(mp.dm_ind==1)); fprintf('  DM1: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm1.Nele, mp.dm1.NactTotal,100*mp.dm1.Nele/mp.dm1.NactTotal); end
            if(any(mp.dm_ind==2)); fprintf('  DM2: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm2.Nele, mp.dm2.NactTotal,100*mp.dm2.Nele/mp.dm2.NactTotal); end
            if(any(mp.dm_ind==8)); fprintf('  DM8: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm8.Nele, mp.dm8.NactTotal,100*mp.dm8.Nele/mp.dm8.NactTotal); end
            if(any(mp.dm_ind==9)); fprintf('  DM9: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm9.Nele, mp.dm9.NactTotal,100*mp.dm9.Nele/mp.dm9.NactTotal); end
            
            %--Crop out unused actuators from the control Jacobian
            if(any(mp.dm_ind==1)); jacStruct.G1 = jacStruct.G1(:,mp.dm1.act_ele,:); end
            if(any(mp.dm_ind==2)); jacStruct.G2 = jacStruct.G2(:,mp.dm2.act_ele,:); end
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
% out.log10regHist(Itr) = cvar.log10regUsed;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Control Algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    fprintf('Computing linearized control matrices from the Jacobian...'); tic;

    %--Compute matrices for linear control with regular EFC
%     if(cvar.flagRelin==true);  cvar.GstarG_wsum  = zeros(cvar.NeleAll,cvar.NeleAll);  end %--Re-use cvar.GstarG_wsum since no re-linearization was done.
    cvar.GstarG_wsum  = zeros(cvar.NeleAll,cvar.NeleAll); 
    cvar.RealGstarEab_wsum = zeros(cvar.NeleAll, 1);

    for im=1:mp.jac.Nmode
    %     si = mp.jac.sbp_inds(im);
    %     zi = mp.jac.zern_inds(im);

        Gstack = [jacStruct.G1(:,:,im), jacStruct.G2(:,:,im), jacStruct.G3(:,:,im), jacStruct.G4(:,:,im), jacStruct.G5(:,:,im), jacStruct.G6(:,:,im), jacStruct.G7(:,:,im), jacStruct.G8(:,:,im), jacStruct.G9(:,:,im)];

        %--Square matrix part stays the same if no re-linearization has occurrred. 
        cvar.GstarG_wsum  = cvar.GstarG_wsum  + mp.jac.weights(im)*real(Gstack'*Gstack); 
        %if(cvar.flagRelin==true);  cvar.GstarG_wsum  = cvar.GstarG_wsum  + mp.jac.weights(im)*real(Gstack'*Gstack);  end

        %--The G^*E part changes each iteration because the E-field changes.
        Eweighted = mp.WspatialVec.*cvar.EfieldVec(:,im); %--Apply 2-D spatial weighting to E-field in dark hole pixels.
        cvar.RealGstarEab_wsum = cvar.RealGstarEab_wsum + mp.jac.weights(im)*real(Gstack'*Eweighted); %--Apply the Jacobian weights and add to the total.

    end
    clear GallCell Gstack Eweighted % save RAM

    %%%%%cvar.GstarG_wsum = cvar.GstarG_wsum(cvar.dm_subset,cvar.dm_subset);
    %%%%cvar.RealGstarEab_wsum = cvar.RealGstarEab_wsum(cvar.dm_subset);
    
    %--Make the regularization matrix. (Define only the diagonal here to save RAM.)
    cvar.EyeGstarGdiag = max(diag(cvar.GstarG_wsum ))*ones(cvar.NeleAll,1);
%     if(cvar.flagRelin==true);   cvar.EyeGstarGdiag = max(diag(cvar.GstarG_wsum ))*ones(cvar.NeleAll,1);  end %--Re-use cvar.GstarG_wsum since no re-linearization was done. %--No relative weighting among the DMs
    fprintf(' done. Time: %.3f\n',toc);

    %--Call the Controller Function
    fprintf('Control beginning ...\n'); tic
    switch mp.controller

        case{'plannedEFC'} %--EFC regularization is scheduled ahead of time
            [dDM,cvar] = falco_ctrl_planned_EFC(mp,cvar);

        case{'gridsearchEFC'}  %--Empirical grid search of EFC. Scaling factor for DM commands too.
            [dDM,cvar] = falco_ctrl_grid_search_EFC(mp,cvar);

        case{'conEFC'} %--Constrained EFC. The quadratic cost function is solved directly with CVX rather than by inverting.
            cvar.dummy = 1;
            [dDM,cvar] = falco_ctrl_EFC_constrained(mp,cvar);

    end
    fprintf(' done. Time: %.3f sec\n',toc);
    

    %% Updates to DM commands

    %--Update the DM commands by adding the new control signal
    if(any(mp.dm_ind==1))
        mp.dm1.dV = dDM.dDM1V;
        mp.dm1.V = mp.dm1.V + mp.dm1.dV; 
    end
    if(any(mp.dm_ind==2))
        mp.dm2.dV = dDM.dDM2V;
        mp.dm2.V = mp.dm2.V + mp.dm2.dV; 
    end
    if(any(mp.dm_ind==8))
        mp.dm8.dV = dDM.dDM8V(:);
        mp.dm8.V = mp.dm8.V + mp.dm8.dV(:);
    end
    if(any(mp.dm_ind==9))
        mp.dm9.dV = dDM.dDM9V;
        mp.dm9.V = mp.dm9.V + mp.dm9.dV;
    end


end %--END OF FUNCTION
