% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to compute the flux ratio noise (FRN) input table (also called a
% Krist table) for CGI analysis.
%
% REVISION HISTORY:
% - Created by A.J. Riggs on 2019-04-30.
% -------------------------------------------------------------------------

function tableKrist = falco_FRN_Krist_table(mp)
    

    %%--Misc setup    
    arcsecPerLamD = mp.lambda0/mp.yield.Dtel*180/pi*3600;

    %--Define radial sampling and offset values
    dr = 1/mp.Fend.res; % step size between pixels[lambda0/D]
    rs = mp.yield.R0:dr:mp.yield.R1; % [lambda0/D]
    Noff = length(rs);
    
    %--Initializations
    matKrist = zeros(Noff,8);
    thput_vec = zeros(Noff,1);
    PSF_peak_vec = zeros(Noff,1);
    area_vec = zeros(Noff,1);
    
    %--For parfor
    %--Off-axis
    vals_list_off = allcomb(rs,1:mp.Nsbp,1:mp.Nwpsbp).';   %--dimensions: [3 x Noff*mp.Nsbp*mp.Nwpsbp ]
    inds_list_off = allcomb(1:Noff,1:mp.Nsbp,1:mp.Nwpsbp).';   %--dimensions: [3 x Noff*mp.Nsbp*mp.Nwpsbp ]
    NvalsOff = max(size(vals_list_off,2));
    %--On-axis
    vals_list_on = allcomb(mp.full.pol_conds,1:mp.Nsbp,1:mp.Nwpsbp).'; %--dimensions: [3 x length(mp.full.pol_conds)*mp.Nsbp*mp.Nwpsbp ]
    inds_list_on = allcomb(1:length(mp.full.pol_conds),1:mp.Nsbp,1:mp.Nwpsbp).'; %--dimensions: [3 x length(mp.full.pol_conds)*mp.Nsbp*mp.Nwpsbp ]
    NvalsOn = max(size(vals_list_on,2));
    
    %% Nested function for use in Sections with OFF-axis PSF calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The function func_offset_wavelength.m exists at all so that parfor
    % can use a linear indexing scheme from 1 to Nvals. 
    % This is a nested function to try to reduce RAM overhead in MATLAB.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:
    % -mp = structure of model parameters
    % -vals_list = structure containing combinations vectors to loop over
    % -jj = index of the loop
    % -flagFull = whether to use the full model or not
    %
    % OUTPUTS:
    % -Iout = Jacobian for the specified combo of DM, wavelength, and Zernike mode.

    function Iout = func_offset_wavelength(mp,vals_list,jj)

        mp.full.polaxis = 10; %--Use the average polarization for off-axis sources because we don't care about the PSF wings
        
        offsetVal = vals_list(1,jj); %--tip/tilt offset
        si = vals_list(2,jj); %--sub-bandpass index
        wi = vals_list(3,jj); %--wavelength index within a sub-bandpass
        
        modvar.whichSource = 'offaxis';
        mp.full.source_x_offset = offsetVal;
        mp.full.source_y_offset = 0;
        modvar.x_offset = offsetVal; % mp.thput_eval_x;
        modvar.y_offset = 0; % mp.thput_eval_y;
        modvar.sbpIndex = si; 
        modvar.wpsbpIndex = wi;
        modvar.zernIndex = 1; %--Piston only
        E2D = model_full(mp, modvar,'Unnorm'); %--un-normalized to PSF peak. Instead normalized to energy at the input pupil
        
        Iout = (abs(E2D).^2)*mp.sbp_weights(si)*mp.full.lambda_weights(wi); %--Include sbp and wavelength weights
        % fprintf('Image generated for      offset = %.2f        sbp %d/%d        wvl %d/%d\n',offsetVal,si,mp.Nsbp,wi,mp.Nwpsbp);
 
    end %--END OF FUNCTION func_offset_wavelength.m    

    funcOffaxis = @(ii) func_offset_wavelength(mp,vals_list_off,ii); %--Make a function handle for parfor to use
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Nested function for use in Sections with ON-axis PSF calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The function func_offset_wavelength.m exists at all so that parfor
    % can use a linear indexing scheme from 1 to Nvals. 
    % This is a nested function to try to reduce RAM overhead in MATLAB.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:
    % -mp = structure of model parameters
    % -vals_list = structure containing combinations vectors to loop over
    % -jj = index of the loop
    %
    % OUTPUTS:
    % -Iout = 2-D normalized intensity image

    function Iout = func_pol_wavelength(mp,vals_list,jj)
        
        mp.full.polaxis = vals_list(1,jj); %--polarization state
        si = vals_list(2,jj); %--sub-bandpass index
        wi = vals_list(3,jj); %--wavelength index within a sub-bandpass
        
        Npol = length(mp.full.pol_conds);
        modvar.whichSource = 'offaxis';
        mp.full.source_x_offset = 0;
        mp.full.source_y_offset = 0;
        modvar.x_offset = 0;
        modvar.y_offset = 0;
        modvar.sbpIndex = si; 
        modvar.wpsbpIndex = wi; 
        modvar.zernIndex = 1;
        E2D = model_full(mp, modvar,'Unnorm'); %--un-normalized to PSF peak. Instead normalized to energy at the input pupil

        Iout = (abs(E2D).^2)*mp.sbp_weights(si)*mp.full.lambda_weights(wi)/Npol;
        % fprintf('Image generated for      polaxis = %.2f        sbp %d/%d        wvl %d/%d\n',mp.full.polaxis,si,mp.Nsbp,wi,mp.Nwpsbp);
 
    end %--END OF FUNCTION func_offset_wavelength.m    

    funcOnaxis = @(ii) func_pol_wavelength(mp,vals_list_on,ii); %--Make a function handle for parfor to use
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Off-axis Calculations: Columns 1,2,5,6,7
    %--Can use polaxis=10 since care only about the PSF peak and not the wings.
    %--Fastest way to compute: Make a single parfor loop over both
    %  wavelength and PSF offset and save all images into a datacube. Then run a loop over radius to compute the various metrics. 
    
    

    %--Parallel/distributed computing
    tic; fprintf('Computing columns 1,2,5,6,7 of the Krist table... ')
    if(mp.flagParfor) 
        parfor ni=1:NvalsOff;  IoffaxisArray{ni} = feval(funcOffaxis, ni);  end    
    else
        for ni=NvalsOff:-1:1;  IoffaxisArray{ni} = feval(funcOffaxis, ni);  end
    end
    fprintf('done. Time = %.2f s\n',toc)

    %--Sum up sub-bandpasses and wavelengths at each radial offset
    IoffaxisCube = zeros(mp.Fend.Neta,mp.Fend.Nxi,Noff);
    for ni=1:NvalsOff
        ioff = inds_list_off(1,ni); %--index of offset
        IoffaxisCube(:,:,ioff) = IoffaxisCube(:,:,ioff) + IoffaxisArray{ni};
    end
    clear IoffaxisArray

    %--Compute table values at each radial offset
    for ioff=1:Noff
        xi_offset = rs(ioff);
        eta_offset = 0;
        Icam = IoffaxisCube(:,:,ioff);
        if(mp.full.flagPROPER==false);  Icam = Icam/mp.sumPupil;  end  %--PROPER already normalizes the energy
        if(mp.flagPlot); figure(324); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,Icam); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end

        %--Peak pixel value
        PSF_peak_vec(ioff) = max(Icam(:));

        %--Absolute energy within half-max isophote(s)
        maskHM = 0*mp.Fend.RHOS;
        maskHM(Icam>=1/2*max(max(Icam))) = 1;
        if(mp.flagPlot);  figure(325); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,maskHM); axis xy equal tight; drawnow;  end
        thput_vec(ioff) = sum(Icam(maskHM==1)); %--half-max throughput
        fprintf('Core throughput within the half-max isophote(s) = %.2f%% \tat separation = (%.1f, %.1f) lambda0/D.\n',100*thput_vec(ioff),xi_offset,eta_offset);

        area_vec(ioff) = (arcsecPerLamD/mp.Fend.res)^2*length(Icam(Icam>=1/2*max(max(Icam)))); %--area of the photometric aperture used for throughput [arcsec^2]
    end
    clear IoffaxisCube
    
    matKrist(:,1) = rs.';                % [lambda0/D]
    matKrist(:,2) = arcsecPerLamD*rs.';  % [arcseconds]
    % Column 3 is intensity
    % Column 4 is contrast
    matKrist(:,5) = thput_vec;
    matKrist(:,6) = PSF_peak_vec;
    matKrist(:,7) = area_vec;
    % Column 8 is the Lyot stop transmission

    % %--Shows that the conversion from intensity to contrast is pretty much the
    % %same whether you use the throughput curve or peak pixel curve
    % figure(198); plot(rs,1./(thput_vec/max(thput_vec)),rs,1./(PSF_peak_vec/max(PSF_peak_vec))); 
    
%% On-axis calculations (intensity,contrast)
    %--Have to use all polarization states since care about PSF wings.
    
    
    %--Parallel/distributed computing
    tic; fprintf('Computing columns 3 and 4 of the Krist table... ')
    if(mp.flagParfor) 
        parfor ii=1:NvalsOn;  IonaxisArray{ii} = feval(funcOnaxis, ii);  end    
    else
        for ii=NvalsOn:-1:1;  IonaxisArray{ii} = feval(funcOnaxis, ii);  end
    end
    fprintf('done. Time = %.2f s\n',toc)
    
    %--Sum up sub-bandpasses
    Icam = 0;
    for ii=1:NvalsOn
        Icam = Icam + IonaxisArray{ii};
    end
    clear IonaxisArray
    if(mp.full.flagPROPER==false); Icam = Icam/mp.sumPupil;  end %--PROPER already normalizes the energy
    if(mp.flagPlot); figure(324); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,log10(Icam)); axis xy equal tight; title('On-axis Intensity','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end             

    %--Average the intensity for each annular sector or annulus.
    int_vec = zeros(Noff,1);
    for ioff=1:Noff        
        min_r = rs(ioff) - dr/2;
        max_r = rs(ioff) + dr/2;
        
        %--Compute the software mask for the scoring region
        maskScore.pixresFP = mp.Fend.res;
        maskScore.rhoInner = min_r; %--lambda0/D
        maskScore.rhoOuter = max_r; %--lambda0/D
        maskScore.angDeg = mp.Fend.score.ang; %--degrees
        maskScore.centering = mp.centering;
        maskScore.FOV = mp.Fend.FOV;
        maskScore.whichSide = mp.Fend.sides; %--which (sides) of the dark hole have open
        if(isfield(mp.Fend,'shape'));  maskScore.shape = mp.Fend.shape;  end
        [maskPartial,xis,etas] = falco_gen_SW_mask(maskScore);
        
        %--Compute the average intensity over the selected region
        int_vec(ioff) = sum(sum(maskPartial.*Icam))/sum(sum(maskPartial));
    end
    c_vec = int_vec./(PSF_peak_vec); %--Convert intensity to contrast

    matKrist(:,3) = int_vec; % [intensity (not normalized)]
    matKrist(:,4) = c_vec; % [raw contrast]
%     figure(195); semilogy(rs,int_vec);  drawnow;
%     figure(197); semilogy(rs,c_vec);  drawnow;


%% Lyot Plane Transmission (Generate column 8 of the Krist table)
    %--Can use polaxis=10 since care only about the total energy.
    %--Fastest way to compute: Make a single parfor loop over both
    %wavelength and PSF offset and save all images into a datacube. Then run a loop over radius to compute the various metrics. 
    
    
    %--Create a new mp as mpTemp to change some values for this column's calculation.
    mp.full.use_field_stop = 0; %--Make sure the field stop is not used
    mpTemp = mp;
    mpTemp.Fend.res = 1; %--Can use low sampling when just summin the total energy in the  plane.
    FOV = 50; %--Use a large enough FOV [lambda0/D]
    mpTemp.full.output_dim = ceil_even(1 + mpTemp.Fend.res*(2*FOV)); %  dimensions of output in pixels (overrides output_dim0)
    mpTemp.full.final_sampling_lam0 = 1/mpTemp.Fend.res;	%   final sampling in lambda0/D
    funcOffaxis2 = @(ii) func_offset_wavelength(mpTemp,vals_list_off,ii); %--Make a function handle for parfor to use
    
    %--Parallel/distributed computing
    tic; fprintf('Computing column 8 of the Krist table (Lyot transmission)... ')
    if(mp.flagParfor) 
        parfor ni=1:NvalsOff;  IoffaxisArray{ni} = feval(funcOffaxis2, ni);  end    
    else
        for ni=NvalsOff:-1:1;  IoffaxisArray{ni} = feval(funcOffaxis2, ni);  end
    end
    fprintf('done. Time = %.2f s\n',toc)

    %--Sum up sub-bandpasses
    IoffaxisCube = zeros(mpTemp.full.output_dim,mpTemp.full.output_dim,Noff);
    for ni=1:NvalsOff
        ioff = inds_list_off(1,ni); %--index of the radial offset
        IoffaxisCube(:,:,ioff) = IoffaxisCube(:,:,ioff) + IoffaxisArray{ni};
    end
    clear IoffaxisArray
           
    %--Get the total intensity leaving the Lyot plane (= at the final focal plane assuming no field stop).
    occ_trans_vec = zeros(Noff,1);
    for ioff=1:Noff 
        Icam = IoffaxisCube(:,:,ioff);
        if(mp.full.flagPROPER==false);  Icam = Icam/mp.sumPupil;  end%--PROPER already normalizes the energy
%         if(mp.flagPlot); figure(324); imagesc(Icam); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end
%         if(mp.flagPlot); figure(325); imagesc(log10(Icam)); axis xy equal tight; title('Off-axis PSF for Throughput Calculation','Fontsize',20); set(gca,'Fontsize',20); colorbar; drawnow;  end
        occ_trans_vec(ioff) = sum(sum(Icam));
    end
    matKrist(:,8) = occ_trans_vec; %--[Lyot stop transmission = total focal plane energy]
    
    %% Make into a table for printing as a CSV file

    tableKrist = table(matKrist(:,1),matKrist(:,2),matKrist(:,3),matKrist(:,4),matKrist(:,5),matKrist(:,6),matKrist(:,7),matKrist(:,8));
    tableKrist.Properties.VariableNames = {'r_lamOverD', 'r_arcsec', 'I', 'contrast', 'core_thruput', 'PSF_peak', 'area_sq_arcsec', 'trans_Lyot'};


end %--END OF FUNCTION