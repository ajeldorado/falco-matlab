% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to perform wavefront sensing and control using an absolute
% Zernike wavefront sensor.  This function is intended to be parallel
% falco_wfsc_loop and follow a similar structure, except for the fact that
% it doesn't do EFC.

function [out, varargout] = falco_Zwfsc_loop(mp)

%% Sort out file paths and save the config file    

%--Add the slash or backslash to the FALCO path if it isn't there.
if( strcmp(mp.path.falco(end),'/')==false && strcmp(mp.path.falco(end),'\')==false )
    mp.path.falco = [mp.path.falco filesep];
end

mp.path.dummy = 1; %--Initialize the folders structure in case it doesn't already exist

%--Store minimal data to re-construct the data from the run: the config files and "out" structure after a trial go here
if(isfield(mp.path,'config')==false)
    mp.path.config = [mp.path.falco filesep 'data' filesep 'config' filesep];     
end

%--Entire final workspace from FALCO gets saved here.
if(isfield(mp.path,'ws')==false)
    mp.path.ws = [mp.path.falco filesep 'data' filesep 'ws' filesep];
end

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config.mat'];
save(fn_config)
fprintf('Saved the config file: \t%s\n',fn_config)

%% Get configuration data from a function file
if(~mp.flagSim);  bench = mp.bench;  end %--Save the testbed structure "mp.bench" into "bench" so it isn't overwritten by falco_init_ws
[mp,out] = falco_init_ws(fn_config);
if(~mp.flagSim);  mp.bench = bench;  end

%% Initializations of Arrays for Data Storage 

%--Raw contrast (broadband)

InormHist = zeros(mp.Nitr,1); % Measured, mean raw contrast in scoring region of dark hole.

%% Take initial broadband image 

Im = falco_get_summed_image(mp);

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
    [mp,thput] = falco_compute_thput(mp);
    if(mp.flagFiber)
        mp.thput_vec(Itr) = max(thput);
    else
        mp.thput_vec(Itr) = thput; %--record keeping
    end
    
    %--Compute the current contrast level
    InormHist(Itr) = mean(Im(mp.Fend.corr.maskBool));

    %--Plot the updates to the DMs and PSF
    if(Itr==1); hProgress.master = 1; end %--dummy value to intialize the handle variable
    if(isfield(mp,'testbed') )
        InormHist_tb.total = InormHist; 
        Im_tb.Im = Im;
        Im_tb.E = zeros(size(Im));
        if(Itr>1)
            InormHist_tb.mod(Itr-1) = mean(abs(EfieldVec(:)).^2);
            InormHist_tb.unmod(Itr-1) = mean(IincoVec(:));
            Im_tb.E(mp.Fend.corr.mask) = EfieldVec(:,ceil(mp.Nsbp/2));
        else
            InormHist_tb.mod = NaN;
            InormHist_tb.unmod = NaN;
        end
        hProgress = falco_plot_progress_gpct(hProgress,mp,Itr,InormHist_tb,Im_tb,DM1surf,DM2surf);
    else
        hProgress = falco_plot_progress(hProgress,mp,Itr,InormHist,Im,DM1surf,DM2surf);
    end

%% WFSC goes here
    
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
Im = falco_get_summed_image(mp);
fprintf('done. Time = %.1f s\n',toc);

%--REPORTING NORMALIZED INTENSITY
InormHist(Itr+1) = mean(Im(mp.Fend.corr.maskBool));
fprintf('Prev and New Measured Contrast (LR):\t\t\t %.2e\t->\t%.2e\t (%.2f x smaller)  \n',...
    InormHist(Itr), InormHist(Itr+1), InormHist(Itr)/InormHist(Itr+1) ); 

fprintf('\n\n');

%--Save out DM commands after each iteration in case the trial crashes part way through.
if(mp.flagSaveEachItr)
    fprintf('Saving DM commands for this iteration...')
    if(any(mp.dm_ind==1)); DM1V = mp.dm1.V; else; DM1V = 0; end
    if(any(mp.dm_ind==2)); DM2V = mp.dm2.V; else; DM2V = 0; end
    if(any(mp.dm_ind==3)); DM3V = mp.dm3.V; else; DM3V = 0; end
    if(any(mp.dm_ind==4)); DM4V = mp.dm4.V; else; DM4V = 0; end
    if(any(mp.dm_ind==5)); DM5V = mp.dm5.V; else; DM5V = 0; end
    if(any(mp.dm_ind==6)); DM6V = mp.dm6.V; else; DM6V = 0; end
    if(any(mp.dm_ind==7)); DM7V = mp.dm7.V; else; DM7V = 0; end
    if(any(mp.dm_ind==8)); DM8V = mp.dm8.V; else; DM8V = 0; end
    if(any(mp.dm_ind==9)); DM9V = mp.dm9.V; else; DM9V = 0; end
    Nitr = mp.Nitr;
    thput_vec = mp.thput_vec;
    fnWS = sprintf('%sws_%s_Iter%dof%d.mat',mp.path.ws_inprogress,mp.runLabel,Itr,mp.Nitr);
    save(fnWS,'Nitr','Itr','DM1V','DM2V','DM3V','DM4V','DM5V','DM6V','DM7V','DM8V','DM9V','InormHist','thput_vec')
    fprintf('done.\n\n')
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
if(any(mp.dm_ind==5)); out.dm5.Vall(:,:,Itr) = mp.dm5.V; end
if(any(mp.dm_ind==8)); out.dm8.Vall(:,Itr) = mp.dm8.V; end
if(any(mp.dm_ind==9)); out.dm9.Vall(:,Itr) = mp.dm9.V; end

%--Calculate the core throughput (at higher resolution to be more accurate)
[mp,thput] = falco_compute_thput(mp);
if(mp.flagFiber)
    mp.thput_vec(Itr) = max(thput);
else
    mp.thput_vec(Itr) = thput; %--record keeping
end

if(isfield(mp,'testbed'))
    InormHist_tb.total = InormHist; 
    InormHist_tb.mod(Itr-1) = mean(abs(EfieldVec(:)).^2);
    InormHist_tb.unmod(Itr-1) = mean(IincoVec(:));
    Im_tb.Im = Im;
    Im_tb.E = zeros(size(Im));
    Im_tb.E(mp.Fend.corr.mask) = EfieldVec(:,ceil(mp.Nsbp/2));
    hProgress = falco_plot_progress_gpct(hProgress,mp,Itr,InormHist_tb,Im_tb,DM1surf,DM2surf);
else
    hProgress = falco_plot_progress(hProgress,mp,Itr,InormHist,Im,DM1surf,DM2surf);
end
%% Optional output variable: mp
varargout{1} = mp;

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

fnOut = [mp.path.config mp.runLabel,'_snippet.mat'];

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

    mp.P1.full.E = 1; mp.P1.compact.E=1; mp.Eplanet=1; 
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