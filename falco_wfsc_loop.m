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

function [] = falco_wfsc_loop(fn_config,varargin)

% Set default values of input parameters
mp.flagPlot = false; % flag to plot PSF correction in real time

%--Enable different arguments values by using varargin
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'plot','plotting'}
        mp.flagPlot = true;  %--Value to use for turning on plots
      otherwise
        error('falco_wfsc_loop: Unknown keyword: %s\n', varargin{icav});
          
    end
end

%% Get configuration data from a function file

[mp,cp,ep,DM,folders] = falco_init_ws(fn_config,mp.flagPlot);

%%
% G_mat_fname = sprintf(
% cd(folders.jac)
% if( exist(G_mat_fname, 'file') ~= 2) %--if 2, then the file by that name exists
%     flagCalcJac = true; %--Calculate the starting Jacobian again if the file does not exist
% end
% cd(folders.m)

%% Initializations of Arrays for Data Storage 
ImBandAvg_array = single( zeros(mp.F4.full.Neta,mp.F4.full.Nxi,mp.Nitr+1) ); %--Full PSF after each correction step

%--Raw contrast (broadband)
contrast_bandavg = zeros(mp.Nitr,1); % Measured/raw contrast in ScoreMask
% contrast_left_bandavg = zeros(mp.Nitr+1,1);
% contrast_right_bandavg = zeros(mp.Nitr+1,1);

%--Store the DM surfaces (REQUIRES LOTS OF STORAGE)
if(any(DM.dm_ind==1)); DM1S_array = single(zeros(DM.dm1.compact.Ndm,DM.dm1.compact.Ndm,mp.Nitr+1)); else; DM1S_array = zeros(2,2,mp.Nitr+1); end
if(any(DM.dm_ind==2)); DM2S_array = single(zeros(DM.dm2.compact.Ndm,DM.dm2.compact.Ndm,mp.Nitr+1)); else; DM2S_array = zeros(2,2,mp.Nitr+1); end

%% Take initial broadband images
EfieldCorrTrue = zeros(length(mp.F4.compact.corr.inds),mp.Nttlam,mp.Nitr+1); % (Simulation only) Vectorized true starlight E-field at each pixel and wavelength

if(mp.flagPlot); figure(101); imagesc(mp.P1.full.mask);axis image; colorbar; title('pupil');drawnow; end

if(mp.flagPlot && (length(mp.P4.full.mask)==length(mp.P1.full.mask))); figure(102); imagesc(mp.P4.full.mask);axis image; colorbar; title('Lyot stop');drawnow; end;
if(mp.flagPlot && isfield(mp,'P3.full.mask')); figure(103); imagesc(padOrCropEven(mp.P1.full.mask,mp.P3.full.Narr).*mp.P3.full.mask);axis image; colorbar; drawnow; end;

%% Take initial broadband image 

[~, ImBandAvg_array(:,:,1)] = falco_est_perfect_Efield_full(mp,DM);
[EfieldCorrTrue(:,:,1), ~]  = falco_est_perfect_Efield_compact(mp,DM);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin the Correction Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Itr=1:mp.Nitr
fprintf(['Iteration: ' num2str(Itr) '/' num2str(mp.Nitr) '\n' ]);

%--Compute the DM surfaces
if(any(DM.dm_ind==1)); DM1surf =  falco_gen_dm_surf(DM.dm1, DM.dm1.compact.dx, DM.dm1.compact.Ndm);  end
if(any(DM.dm_ind==2)); DM2surf =  falco_gen_dm_surf(DM.dm2, DM.dm2.compact.dx, DM.dm2.compact.Ndm);  end

%--Data to store
if(any(DM.dm_ind==1)); DM.dm1.Vall(:,:,Itr) = DM.dm1.V; DM1S_array(:,:,Itr) = single(DM1surf); end
if(any(DM.dm_ind==2)); DM.dm2.Vall(:,:,Itr) = DM.dm2.V; DM2S_array(:,:,Itr) = single(DM2surf); end


%--Calculate the core throughput (at same resolution as in models)

modvar.x_offset = mp.thput_eval_x;
modvar.y_offset = mp.thput_eval_y;
modvar.sbpIndex = mp.si_ref; 
modvar.whichSource = 'offaxis';
ImTemp = falco_sim_image(mp, modvar, DM);
if(mp.flagPlot); figure(324); imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,ImTemp); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end

maskHM = 0*mp.FP4.compact.RHOS;
maskHM(ImTemp>=1/2*max(max(ImTemp))) = 1;
mp.maskHMcore = maskHM.*mp.maskCore;
% figure(325); imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,mp.maskCore); axis xy equal tight; drawnow;
mp.thput_vec(Itr) = sum(ImTemp(mp.maskHMcore==1))/mp.sumPupil*mp.F4.full.I00(mp.si_ref);
fprintf('>=Half-max core throughput within a %.2f lambda/D radius = %.2f%% \tat separation = (%.1f, %.1f) lambda/D.\n',mp.thput_radius,100*mp.thput_vec(Itr),mp.thput_eval_x,mp.thput_eval_y);
%------------------

%--Re-compute the starlight normalization factor for the compact and full models 
%  (to convert images to normalized intensity). No tip/tilt necessary.
mp = falco_get_PSF_norm_factor(mp, DM);

%------------------

% Relinearize about the DMs only at the iteration numbers in mp.relinItrVec.
if(any(mp.relinItrVec==Itr))
    flagRelin=true;
else
    flagRelin=false;
end

%--Calculate (or re-calculate) the control Jacobian for each DM.
%    Relinearize about the DMs only at the iteration numbers in mp.relinItrVec.
% if( (Itr==1 && flagCalcJac==1) || (flagRelin==1 && Itr>1) )
if( (Itr==1) || (flagRelin==true) )
    
    modvar.flagCalcJac = true; 
    modvar.wpsbpIndex = mp.wi_ref;
    modvar.whichSource = 'star'; 
    
    %--Re-initialize the Jacobian arrays to full size
    G1=zeros(1,1,mp.Nttlam); G2=zeros(1,1,mp.Nttlam); G3=zeros(1,1,mp.Nttlam); G4=zeros(1,1,mp.Nttlam); G5=zeros(1,1,mp.Nttlam); G6=zeros(1,1,mp.Nttlam); G7=zeros(1,1,mp.Nttlam); G8=zeros(1,1,mp.Nttlam); G9=zeros(1,1,mp.Nttlam); %--Initialize for bookkeeping in cells later 
    if(any(DM.dm_ind==1)); G1 = zeros(length(mp.F4.compact.corr.inds),DM.dm1.NactTotal,mp.Nttlam); end % control Jacobian for DM1
    if(any(DM.dm_ind==2)); G2 = zeros(length(mp.F4.compact.corr.inds),DM.dm2.NactTotal,mp.Nttlam); end % control Jacobian for DM2
    
%--Compute the number of total actuators for all DMs used. 
if(Itr==1)
    GallCell1 = {squeeze(G1(:,:,1)),squeeze(G2(:,:,1)),squeeze(G3(:,:,1)),squeeze(G4(:,:,1)),squeeze(G5(:,:,1)),squeeze(G6(:,:,1)),squeeze(G7(:,:,1)),squeeze(G8(:,:,1)),squeeze(G9(:,:,1))}; % Create the cell array. Placeholders for non-existent Jacobians to have consistent numbering
    NeleAll = 0;
    NeleVec = []; %--Vector of total number of used actuators for each used DM
    for ii=1:numel(DM.dm_ind)
        dm_index = DM.dm_ind(ii);
        NeleAll = NeleAll + size(GallCell1{dm_index},2);
        NeleVec = [NeleVec; size(GallCell1{dm_index},2) ];
    end
    clear GallCell1
end

%--Compute the control Jacobians for each DM
jacStruct =  model_Jacobian(mp, DM);
if(any(DM.dm_ind==1)); G1 = jacStruct.G1; end
if(any(DM.dm_ind==2)); G2 = jacStruct.G2; end
if(any(DM.dm_ind==9)); G9 = jacStruct.G9; end
clear jacStruct  %--Save RAM



%--MOVE TO A FUNCTION
%--Reduce the number of actuators used based on their relative strength in the first Jacobian
if(Itr==1)
    fprintf('Weeding out weak actuators from the control Jacobian...\n'); 
    if(any(DM.dm_ind==1))
        %--Crop out very weak-effect actuators
        G1intNorm = zeros(DM.dm1.Nact);
        G1intNorm(1:end) = sum( abs(G1(:,:,1)).^2, 1);
        G1intNorm = G1intNorm/max(max(G1intNorm));
        DM.dm1.act_ele = find(G1intNorm>=10^(mp.logGmin));
        %if(mp.flagPlot); figure(81); imagesc(log10(G1int),[-6 0]); axis xy equal tight; colorbar; end
    end
        if(any(DM.dm_ind==2))
        G2intNorm = zeros(DM.dm2.Nact);
        G2intNorm(1:end) = sum( abs(G2(:,:,1)).^2,1);
        G2intNorm = G2intNorm/max(max(G2intNorm));
        DM.dm2.act_ele = find(G2intNorm>=10^(mp.logGmin));
        %if(mp.flagPlot); figure(82); imagesc(log10(G2int),[-6 0]); axis xy equal tight; colorbar; end
        end


    %--Update the number of elements used per DM
    if(any(DM.dm_ind==1)); DM.dm1.Nele = length(DM.dm1.act_ele); end
    if(any(DM.dm_ind==2)); DM.dm2.Nele = length(DM.dm2.act_ele); end
    DM.NelePerDMvec = [length(DM.dm1.Nele), length(DM.dm2.Nele), length(DM.dm3.Nele), length(DM.dm4.Nele), length(DM.dm5.Nele), length(DM.dm6.Nele), length(DM.dm7.Nele), length(DM.dm8.Nele), length(DM.dm9.Nele) ];

    %if(tsi==1);
        if(any(DM.dm_ind==1)); fprintf('  DM1: %d/%d (%.2f%%) actuators kept for Jacobian\n', DM.dm1.Nele, DM.dm1.NactTotal,100*DM.dm1.Nele/DM.dm1.NactTotal); end
        if(any(DM.dm_ind==2)); fprintf('  DM2: %d/%d (%.2f%%) actuators kept for Jacobian\n', DM.dm2.Nele, DM.dm2.NactTotal,100*DM.dm2.Nele/DM.dm2.NactTotal); end
    %end
end    


%--Crop out unused actuators from the control Jacobian
if(any(DM.dm_ind==1)); G1 = G1(:,DM.dm1.act_ele,:); end
if(any(DM.dm_ind==2)); G2 = G2(:,DM.dm2.act_ele,:); end

% Add spatially-dependent weighting to the control Jacobians
if(any(DM.dm_ind==1)); G1 = G1.*repmat(mp.Wspatial_ele,[1,DM.dm1.Nele,mp.Nttlam]); end
if(any(DM.dm_ind==2)); G2 = G2.*repmat(mp.Wspatial_ele,[1,DM.dm2.Nele,mp.Nttlam]); end  

fprintf('Total Jacobian Calcuation Time: %.2f\n',toc);
    

%--Compute the number of total actuators for all DMs used. 
GallCell1 = {squeeze(G1(:,:,1)),squeeze(G2(:,:,1)),squeeze(G3(:,:,1)),squeeze(G4(:,:,1)),squeeze(G5(:,:,1)),squeeze(G6(:,:,1)),squeeze(G7(:,:,1)),squeeze(G8(:,:,1)),squeeze(G9(:,:,1))}; % Create the cell array. Placeholders for non-existent Jacobians to have consistent numbering
NeleAll = 0; %--Number of total actuators used 
NeleVec = []; %--Vector of total number of used actuators for each used DM
Nele12 = 0; %--Number of total actuators used (for DMs 1 and 2 only)
Nele12Vec = []; %--Vector of total number of used actuators for each used DM (of DMs 1 and 2)
for ii=1:numel(DM.dm_ind)
    dm_index = DM.dm_ind(ii);
    NeleAll = NeleAll + size(GallCell1{dm_index},2);
    NeleVec = [NeleVec; size(GallCell1{dm_index},2) ]; %--Number of elements per mode
    if(dm_index==1 || dm_index==2)
        Nele12 = Nele12 + size(GallCell1{dm_index},2);
        Nele12Vec = [Nele12Vec; size(GallCell1{dm_index},2) ]; %--Number of elements per mode
    end
end
clear GallCell1

%     if(Itr==1)
%         cd(folders.jac)
% %             save(G_mat_fname,'G1','G2','G9','-v7.3');
% %             save(G_mat_fname,'G9','-v7.3');
%         cd(folders.m)
%     end

% elseif( Itr==1 && flagCalcJac==0 )    
%     cd(folders.jac)
%         %load(G_mat_fname)
%     cd(folders.m)    
end

%--Compute the current contrast level
ImBandAvg_current = ImBandAvg_array(:,:,Itr);
contrast_bandavg(Itr) = mean(ImBandAvg_current(mp.F4.full.corr.inds));
% contrast_left_bandavg(Itr) = mean(ImBandAvg_current(mp.cor_ele_left));
% contrast_right_bandavg(Itr) = mean(ImBandAvg_current(mp.cor_ele_right));


%% Plot the updates to the DMs and PSF
if(Itr==1); hProgress.master = 1; end %--dummy value to start the handle
hProgress = falco_plot_progress(hProgress,mp,Itr,contrast_bandavg,ImBandAvg_array,DM1S_array,DM2S_array);

%% MOVE THIS TO THE CONTROLLER FUNCTION
%-----------------------------------------------------------------------------------------

% switch mp.controller
% 
%     case{'EFC'}

        %--Compute matrices for linear control with regular EFC
        GstarG_wsum = zeros(NeleAll,NeleAll);
        RealGstarEab_wsum = zeros(NeleAll, 1);
        
        for tsi=1:mp.Nttlam
            GallCell = {squeeze(G1(:,:,tsi)),squeeze(G2(:,:,tsi)),squeeze(G3(:,:,tsi)),squeeze(G4(:,:,tsi)),squeeze(G5(:,:,tsi)),squeeze(G6(:,:,tsi)),squeeze(G7(:,:,tsi)),squeeze(G8(:,:,tsi)),squeeze(G9(:,:,tsi))}; % Create the cell array. Placeholders for non-existent Jacobians to have consistent numbering
            Gstack = []; %--For all DMs
            for ii=1:numel(DM.dm_ind)
                dm_index = DM.dm_ind(ii);
                Gtemp = GallCell{dm_index};
                Gstack = [Gstack, Gtemp]; % stacked control Jacobian at this wavelength and T/T setting. Dimensions: [Npix x NeleAll]

            end
            GstarG_wsum = GstarG_wsum + mp.WttlamVec(tsi)*real(Gstack'*Gstack);
            RealGstarEab_wsum = RealGstarEab_wsum + mp.WttlamVec(tsi)*real(Gstack'*squeeze(EfieldCorrTrue(:,tsi,Itr)));  
        end
        clear GallCell Gtemp Gstack Gstack12 Gstack9% save RAM

        %--Make the diagonal regularization matrix. (Define only the diagonal here
        % to save RAM.)
        %--Decide on the relative DM weighting within the controller's regularization matrix.
        
        EyeGstarGdiag = max(diag(GstarG_wsum))*ones(NeleAll,1);
       
        clear GstarG_wsum12

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
cvar.Itr = Itr; 
cvar.contrast_bandavg = contrast_bandavg(Itr);
cvar.RealGstarEab_wsum = RealGstarEab_wsum;
cvar.EyeGstarGdiag = EyeGstarGdiag;
cvar.GstarG_wsum = GstarG_wsum;
cvar.NeleVec = NeleVec;

fprintf('Control beginning ...\n'); tic
switch mp.controller

    case{'EFC'}  %--Empirical grid search of EFC. Weighted controller for tip/tilt modes. Scaling factor for DM commands too.   #NEWFORTIPTILT

        [dDM,cp,cvar] = falco_ctrl_EFC(DM, cp, cvar, mp);
        %dDMvec = -cp.dmfacBest(Itr)*(diag(EyeGstarGdiag)/cp.muBest(Itr) + GstarG_wsum)\RealGstarEab_wsum;  
        
    case{'conEFC'} %--Constrained EFC. The quadratic cost function is solved directly with CVX rather than by inverting.
        cvar.dummy = 1;
        dDM = falco_ctrl_EFC_constrained(DM, cvar);
        
%         Efield = squeeze(EfieldCorrTrue(:,:,Itr));
%         cvar.dummy = 1;
%         GsCell = {G1,G2,G3,G4,G5,G6,G7,G8,G9};
%         dDM = falco_constrained_EFC(mp, DM, cp, cvar, Efield, GsCell);
%         clear GsCell;
        
    case{'RLP'} %--Robust linear program (NOT FINISHED YET)
        %dDMvec = falco_ctrl_robustLP(x, G, deltaG, P)

end
fprintf(' done. Time: %.3f\n',toc);


%--Update the DM commands with the new control signal
if(any(DM.dm_ind==1))
    DM.dm1.dV = dDM.dDM1V;
    DM.dm1.V = DM.dm1.V + DM.dm1.dV; 
end
if(any(DM.dm_ind==2))
    DM.dm2.dV = dDM.dDM2V;
    DM.dm2.V = DM.dm2.V + DM.dm2.dV; 
end


switch mp.controller
    case 'EFC'
        %--Remove railed actuators from the basis set
        if(any(DM.dm_ind==1)); DM.dm1.act_ele = intersect( DM.dm1.act_ele, find( (DM.dm1.V > -DM.dm1.maxAbsV) & (DM.dm1.V < DM.dm1.maxAbsV) ) ); end
        if(any(DM.dm_ind==2)); DM.dm2.act_ele = intersect( DM.dm2.act_ele, find( (DM.dm2.V > -DM.dm2.maxAbsV) & (DM.dm2.V < DM.dm2.maxAbsV) ) ); end
end
%--Update the number of elements used per DM
if(any(DM.dm_ind==1)); DM.dm1.Nele = length(DM.dm1.act_ele); end
if(any(DM.dm_ind==2)); DM.dm2.Nele = length(DM.dm2.act_ele); end

%--Calculate and report updated P-V DM commands.
if(any(DM.dm_ind==1))
    DM.dm1.Vpv(Itr) = (max(max(DM.dm1.V))-min(min(DM.dm1.V)));
    Nrail1 = length(find( (DM.dm1.V <= -DM.dm1.maxAbsV) | (DM.dm1.V >= DM.dm1.maxAbsV) ));
    fprintf(' DM1 P-V in volts: %.3f\t\t%d/%d (%.1f%%) railed actuators \n', DM.dm1.Vpv(Itr), Nrail1, DM.dm1.NactTotal, 100*Nrail1/DM.dm1.NactTotal); 
% %     DM.dm1.maxAbsV(Itr) = max(max(abs(DM.dm1.V)));
end
if(any(DM.dm_ind==2))
    DM.dm2.Vpv(Itr) = (max(max(DM.dm2.V))-min(min(DM.dm2.V)));
    Nrail2 = length(find( (DM.dm2.V <= -DM.dm2.maxAbsV) | (DM.dm2.V >= DM.dm2.maxAbsV) ));
    fprintf(' DM2 P-V in volts: %.3f\t\t%d/%d (%.1f%%) railed actuators \n', DM.dm2.Vpv(Itr), Nrail2, DM.dm2.NactTotal, 100*Nrail2/DM.dm2.NactTotal); 
% %     DM.dm2.maxAbsV(Itr) = max(max(abs(DM.dm2.V)));
end


% Take the next image to check the contrast level (in simulation only)
[~, ImBandAvg_array(:,:,Itr+1)] = falco_est_perfect_Efield_full(mp,DM);
[EfieldCorrTrue(:,:,Itr+1), ~]  = falco_est_perfect_Efield_compact(mp,DM);


if( (Itr==mp.Nitr) || strcmpi(mp.controller,'conEFC') )
    ImBandAvg_current = ImBandAvg_array(:,:,Itr+1);
    contrast_bandavg(Itr+1) = mean(ImBandAvg_current(mp.F4.full.corr.inds));

    fprintf('Prev and New Measured Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
        contrast_bandavg(Itr), contrast_bandavg(Itr+1), contrast_bandavg(Itr)/contrast_bandavg(Itr+1) ); 
else
    fprintf('Prev and New Measured Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
        contrast_bandavg(Itr), cvar.cMin, contrast_bandavg(Itr)/cvar.cMin ); 
end

fprintf('\n');

%--Save out DM commands in case the run crashes part way through.
fprintf('Saving DM commands for this iteration...')
cd(folders.ws_inprogress)

        if(any(DM.dm_ind==1)); DM1V = DM.dm1.V; else DM1V = 0; end
        if(any(DM.dm_ind==2)); DM2V = DM.dm2.V; else DM2V = 0; end
        if(any(DM.dm_ind==3)); DM3V = DM.dm3.V; else DM3V = 0; end
        if(any(DM.dm_ind==8)); DM8V = DM.dm8.V; else DM8V = 0; end
        if(any(DM.dm_ind==9)); DM9V = DM.dm9.V; else DM9V = 0; end
        Nitr = mp.Nitr;
        thput_vec = mp.thput_vec;

        fnWS = sprintf('ws_%s_Iter%dof%d.mat',mp.runLabel,Itr,mp.Nitr);
        save(fnWS,'Nitr','Itr','DM1V','DM2V','DM3V','DM8V','DM9V','contrast_bandavg','thput_vec')
cd(folders.m)
fprintf('done.\n')


end %--END OF ESTIMATION + CONTROL LOOP
%-----------------------------------------------------------------------------------------

%% Update plot one last time
Itr = Itr + 1;

%--Compute the DM surfaces
if(any(DM.dm_ind==1)); DM1surf =  falco_gen_dm_surf(DM.dm1, DM.dm1.compact.dx, DM.dm1.compact.Ndm);  end
if(any(DM.dm_ind==2)); DM2surf =  falco_gen_dm_surf(DM.dm2, DM.dm2.compact.dx, DM.dm2.compact.Ndm);  end

%--Data to store
if(any(DM.dm_ind==1)); DM.dm1.Vall(:,:,Itr) = DM.dm1.V; DM1S_array(:,:,Itr) = single(DM1surf); end
if(any(DM.dm_ind==2)); DM.dm2.Vall(:,:,Itr) = DM.dm2.V; DM2S_array(:,:,Itr) = single(DM2surf); end

modvar.x_offset = mp.thput_eval_x;
modvar.y_offset = mp.thput_eval_y;
modvar.sbpIndex = mp.si_ref; 
modvar.whichSource = 'offaxis';
ImTemp = falco_sim_image(mp, modvar, DM);
if(mp.flagPlot); figure(324); imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,ImTemp); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end

maskHM = 0*mp.FP4.compact.RHOS;
maskHM(ImTemp>=1/2*max(max(ImTemp))) = 1;
mp.maskHMcore = maskHM.*mp.maskCore;
% figure(325); imagesc(mp.F4.full.xisDL,mp.F4.full.etasDL,mp.maskCore); axis xy equal tight; drawnow;
mp.thput_vec(Itr) = sum(ImTemp(mp.maskHMcore==1))/mp.sumPupil*mp.F4.full.I00(mp.si_ref);
fprintf('>=Half-max core throughput within a %.2f lambda/D radius = %.2f%% \tat separation = (%.1f, %.1f) lambda/D.\n',mp.thput_radius,100*mp.thput_vec(Itr),mp.thput_eval_x,mp.thput_eval_y);

hProgress = falco_plot_progress(hProgress,mp,Itr,contrast_bandavg,ImBandAvg_array,DM1S_array,DM2S_array);

%% Save out the data from the workspace
clear cvar GstarG_wsum G1 G2 G9 Gs RealGstarEab* h*; % Save a ton of space when storing the workspace


% Don't bother saving the large 2-D, floating point maps in the workspace (they take up too much space)
mp.P1.full.mask=1; mp.P3.full.mask=1; mp.P3.compact.mask = 1; mp.P1.full.E = 1; mp.Eplanet=1; mp.P1.compact.E=1; mp.dm1.full.mask = 1; mp.dm1.compact.mask = 1; mp.dm2.full.mask = 1; mp.dm2.compact.mask = 1;
DM.dm1.compact.inf_datacube = 0;
DM.dm2.compact.inf_datacube = 0;
mp.P1.compact.mask=1; mp.P1.compact.E=1; mp.LScompact = 1;
clear DM_config_temp

cd(folders.ws)
fnWS = ['ws_',mp.runLabel,'_',num2str(mp.Nitr),'its.mat'];

disp(['Saving workspace to file ' fnWS '...'])
save(fnWS);
disp('done.')
cd(folders.m)



end %--END OF main FUNCTION
