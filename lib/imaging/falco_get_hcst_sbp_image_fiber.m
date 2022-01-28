% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to get an image in the specified sub-bandpass from the Caltech
% HCST testbed. This function will need to be replaced in order to run on a
% different testbed. Note that the number of pixels per lambda*F# is
% predetermined. 
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - si = index of sub-bandpass for which to take the image
%
% OUTPUTS
% - normI: Normalized intensity in the sub-bandpass
%          	(i.e. approximate raw contrast but normalized 
%           by a photometry measurement at a single offset)
%
% REVISION HISTORY
% - Modified 2020-02-11 by C. Coker to support single fiber imaging
% - Created on 2019-02-22 by G. Ruane 

function normI = falco_get_hcst_sbp_image_fiber(mp,si)

    bench = mp.bench;
    sbp_width = bench.info.sbp_width(si); %--Width of each sub-bandpass on testbed (meters)
    sbp_texp  = bench.info.sbp_texp(si);% Exposure time for each sub-bandpass (seconds)
    SMFInt0   = bench.info.SMFInt0s(si);% counts per second 
    
    % ND filter calibration
    if numel(mp.NDfilter_cal)>1
        NDfilter_cal = interp1(mp.NDfilter_wave_cal,mp.NDfilter_cal,mp.sbp_centers(si)*1e9,'linear','extrap');
    else
        NDfilter_cal=mp.NDfilter_cal;
    end

    %----- Send commands to the DM -----
    disp('Sending current DM voltages to testbed') 

    
    %--Quantization of DM actuation steps based on least significant bit of the
    % DAC (digital-analog converter). In height, so called HminStep_tb
    % If HminStep_tb (minimum step in H) is defined, then quantize the DM voltages
    % This is for simulating the LSB effect (Dan's experiment) 
    if(isfield(mp.dm1,'HminStep_tb') && ~any(isnan(mp.dm1.HminStep_tb(:))))

        % If desired method is not defined, set it to the default. 
        if(~isfield(mp.dm1,'HminStepMethod')); mp.dm1.HminStepMethod = 'round'; end

        % Discretize/Quantize the DM voltages (creates dm.Vquantized)
        mp.dm1 = falco_discretize_dm_surf(mp.dm1, mp.dm1.HminStepMethod, 'tb');
        mp.dm1.V = mp.dm1.Vquantized;

    end
    
    
    map = (mp.dm1.V'); % There's a transpose between Matlab and BMC indexing
%     map = flipud(map); % There's a transpose between Matlab and BMC indexing
%     
    % Send the commands to the DM. 
    % Notes: bench.DM.flatvec contains the commands to flatten the DM. 
    %        mp.dm1.V is added to the flat commands inside
    %        hcst_DM_apply2Dmap. 
    %        FALCO allows for searches through multiplicative gains which
    %        are handled outside of this function. 
    cmds = hcst_DM_apply2Dmap(bench,map,1);% Returns actual DM commands 
    
    %Change exp time
    if ~mp.est.flag_performingEst
        if sbp_texp~=bench.andor.tint;hcst_andor_setExposureTime(bench,sbp_texp);end
    end
    
    if(isfield(bench.info,'source') && strcmpi(bench.info.source,'nkt'))
        %----- Get image from the testbed -----
        disp(['Getting image from testbed in band ',num2str(si)])
    
        % Set wavelength
        lam0 = mp.sbp_centers(si);
        lam1 = lam0 - sbp_width/2;
        lam2 = lam0 + sbp_width/2;
        tb_NKT_setWvlRange(bench,lam1*1e9,lam2*1e9);
        pause(0.5)
    else
        disp('Getting image from testbed (using laser source)')
    end
    
    if si==1; resPos = hcst_FEUzaber_move(bench,bench.FEUzaber.posIn);end
    
    if ~mp.flagUseCamera4EFCSMF
        pause(1)
        FemtoV = hcst_readFemtoOutput_adaptive_inV(bench,bench.Femto.averageNumReads); % read voltage from Femto
        normI = FemtoV/SMFInt0;
        %     FemtoV = hcst_readFemtoOutput_inV(bench,bench.Femto.averageNumReads); % read voltage from Femto
        if normI<0
            while normI<0
                disp('Negative reading out of the Femto')
                averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,1e-9);
                FemtoV = hcst_readFemtoOutput_adaptive_inV(bench,averageNumReads);
                normI = FemtoV/SMFInt0;
            end
        elseif(normI<1.3e-8 && bench.Femto.averageNumReads<1000)
            averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,normI);
            FemtoV = hcst_readFemtoOutput_adaptive_inV(bench,averageNumReads);
            normI = FemtoV/SMFInt0;
            while normI<0
                disp('Negative reading out of the Femto')
                averageNumReads =  hcst_fiu_computeNumReadsNeeded(bench,1e-9);
                FemtoV = hcst_readFemtoOutput_adaptive_inV(bench,averageNumReads);
                normI = FemtoV/SMFInt0;
            end
        end
    else

        hcst_andor_setSubwindow(bench,bench.andor.FEURow,...
            bench.andor.FEUCol,128,false);
        dark4EFCSMF = hcst_andor_loadDark(bench,[bench.info.path2darks,'dark_tint',num2str(bench.andor.tint,2),'_coadds1.fits']);
        im = hcst_andor_getImage(bench)-dark4EFCSMF;
        Vsmf = hcst_fiu_aperturePhotometryOnAndor(bench,im,true);%max(im(:));
        normI = Vsmf/(SMFInt0/mp.peakPSFtint(si)*bench.andor.tint*NDfilter_cal);

%         % plot image to see how it looks lik
        figure(111)
%         subplot(2,2,1);
        imagesc((im))
        axis image
        colorbar
        title(['Wvl: ',num2str(si),'/',num2str(mp.Nsbp)])
%         text(8,25,['# negative: ',num2str(num_neg)],'Color','w','FontSize',16);
% 
%         subplot(2,2,3); 
%         if si==1;ma_im_arr = zeros(mp.Nsbp,1)*nan;end
%         ma_im_arr(si)=max(im(:));
%         plot(bench.info.sbp_width,ma_im_arr)
%         xlabel(['wavelength'])
%         ylabel(['max counts image'])
%         
%         subplot(2,2,4); 
%         plot(bench.info.sbp_width,bench.info.sbp_texp)
%         xlabel(['wavelength'])
%         ylabel(['tint [sec]'])
%         drawnow
    end

    
%     % Load the dark with the correct tint. It must exist in the dark
%     % library. 
%     dark = hcst_andor_loadDark(bench,[bench.info.path2darks,'dark_tint',num2str(bench.andor.tint,2),'_coadds1.fits']);
%     
%     % Scale the PSF photometry by the current integration time
%     peakPSF = mp.peakPSF/mp.peakPSFtint*bench.andor.tint*mp.NDfilter_cal; 
%     
%     % Get normalized intensity (dark subtracted and normalized by peakPSF)
%     normI = (hcst_andor_getImage(bench)-dark)/peakPSF; 
%     disp(['***************',num2str(max(map(:)))])
end 
