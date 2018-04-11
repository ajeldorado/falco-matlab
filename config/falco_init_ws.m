% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to provide input parameters as structures
% Setup parameters for a mock SPLC at the HCIT.
%
% Created by A.J. Riggs on 2017-10-31.

function [mp,cp,ep,DM,folders] = falco_init_ws(fn_config, flagPlot)

%% File Paths
%(Defined such that the entire folder can be downloaded anywhere and run without a problem)
mainPath = pwd;
addpath(genpath(mainPath));

folders.m = mainPath;
folders.maps = [mainPath filesep 'maps' filesep];      % Maps go here
% folders.init = [mainPath '/init'];  % Store initialization maps and files here
folders.jac = [mainPath filesep 'data' filesep 'jac' filesep];    % Store the control Jacobians here
folders.images = [mainPath filesep 'data' filesep 'images' filesep];  % Store all full, reduced images here
folders.dm = [mainPath filesep 'data' filesep 'DMmaps' filesep];      % Store DM command maps here
folders.ws = [mainPath filesep 'data' filesep 'ws' filesep];      % Store final workspace data here
folders.ws_inprogress = [mainPath filesep 'data' filesep 'ws_inprogress' filesep];      % Store in progress workspace data here
folders.DoNotPush = [mainPath filesep 'DoNotPush' filesep]; % For extraneous or very large files not to push with git. (E.g., generated mask representations) 
folders.DoesNotSync = [mainPath filesep 'DoesNotSync' filesep]; % For extraneous or very large files not to sync with git. (E.g., prototype functions and scripts) 

%% Read inputs as structures from a .mat config file
load(fn_config);
mp.flagPlot = flagPlot;

disp(['DM 1-2 Fresnel number = ',num2str((mp.P2.D/2)^2/(mp.d_dm1_dm2*mp.lambda0))]);

%% DM
DM.num_dms = length(DM.dm_ind);%2; % 1 or 2, number of DMs to use.

%% Bandwidth and Wavelength Specs
mp.si_ref = ceil(mp.Nsbp/2);
mp.wi_ref = ceil(mp.Nwpsbp/2);
mp.fracBWsbp = mp.fracBW/mp.Nsbp;
fprintf(' Using %d discrete wavelength(s) in each of %d sub-bandpasses over a %.1f%% total bandpass \n\n', mp.Nwpsbp, mp.Nsbp,100*mp.fracBW) ;

%% Tip/Tilt and Spatial Weighting of the Control Jacobian  #NEWFORTIPTILT
mp.ti_ref = ceil(mp.Ntt/2);
mas2lam0D = 1/(mp.lambda0/mp.P1.D*180/pi*3600*1000); %--Conversion factor: milliarcseconds (mas) to lambda0/D
%--Define the (x,y) values for each tip/tilt offset in units of lambda0/D
if(mp.Ntt == 5)
    mp.ttx = mas2lam0D*mp.TToffset*[0,cosd(0),cosd(90),cosd(180),cosd(270)];
    mp.tty = mas2lam0D*mp.TToffset*[0,sind(0),sind(90),sind(180),sind(270)];
elseif(mp.Ntt == 4)
    mp.ttx = mas2lam0D*mp.TToffset*[0,cosd(0),cosd(120),cosd(240)];
    mp.tty = mas2lam0D*mp.TToffset*[0,sind(0),sind(120),sind(240)];
elseif(mp.Ntt == 1)
    mp.ttx = 0;
    mp.tty = 0;
else
    disp('######## ERROR: Number of tip-tilt modes not specified properly. ##########')
    return    
end

mp.Wttlam = zeros(mp.Nsbp,mp.Ntt); %--Initialize weighting matrix of each tip/tilt-wavelength mode for the controller
if(isinf(mp.NlamForTT))
    mp.Wttlam = ones(mp.Nsbp,mp.Ntt); %--Full usage and equal weighting for all T/T's and sub-bandpasses.
elseif(mp.NlamForTT==3)
    %--Set tip/tilt offsets at the middle and both end sub-bandpasses.
    mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
    mp.Wttlam(1,:) = ones(1,mp.Ntt);
    mp.Wttlam(mp.si_ref,:) = ones(1,mp.Ntt);
    mp.Wttlam(end,:) = ones(1,mp.Ntt);
elseif(mp.NlamForTT==2)
    %--Set tip/tilt offsets at only both end sub-bandpasses.
    mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
    mp.Wttlam(1,:) = ones(1,mp.Ntt);
    mp.Wttlam(end,:) = ones(1,mp.Ntt);
elseif(mp.NlamForTT==1)
    %--Set tip/tilt offsets at only the middle sub-bandpass.
    mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
    mp.Wttlam(mp.si_ref,:) = ones(1,mp.Ntt);
elseif(mp.NlamForTT==0)
    %--Set tip/tilt offsets at no sub-bandpasses.
    mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
end

mp.Wsum = sum(sum(mp.Wttlam)); %--Sum of all the control Jacobian weights

mp.Wttlam_ele = find(mp.Wttlam>0); %--Indices of the non-zero control Jacobian modes in the weighting matrix
mp.WttlamVec = mp.Wttlam(mp.Wttlam_ele); %--Vector of control Jacobian mode weights
mp.Nttlam = length(mp.Wttlam_ele); %--Number of modes in the control Jacobian

%--Get the wavelength indices for the nonzero values in the weight matrix. 
temp = (1:mp.Nsbp).';
tempMat = repmat(temp,[1,mp.Ntt]);
mp.Wttlam_si = tempMat(mp.Wttlam_ele);

%--Get the tip/tilt indices for the nonzero values in the weight matrix. 
temp = 1:mp.Ntt;
tempMat = repmat(temp,[mp.Nsbp,1]);
mp.Wttlam_ti = tempMat(mp.Wttlam_ele);

%--Chromatic weighting


%% Coronagraphic Masks

%% Input pupil plane resolution, masks, and coordinates
%--Resolution at input pupil and DM1 and DM2
mp.P2.full.dx = mp.P2.D/mp.P1.full.Nbeam; 
mp.P2.compact.dx = mp.P2.D/mp.P1.compact.Nbeam;
% mp.Npup1factor = mp.P1.full.Nbeam/mp.P1.compact.Nbeam; %--Scaling factor between full model and compact model of the input pupil array
% % % mp.Npup1factor = mp.P1.full.Narr/mp.P1.compact.Narr; %--Scaling factor between full model and compact model of the input pupil array
% % % mp.P2.compact.dx = mp.Npup1factor*mp.P2.full.dx;


%--Generate/Load Input Pupil
switch mp.whichPupil
    case{'Simple'}
        
        inputs.Nbeam = mp.P1.full.Nbeam;     % number of points across the pupil diameter
        inputs.OD = mp.P1.ODnorm;
        inputs.ID = mp.P1.IDnorm;
        inputs.num_strut = mp.P1.num_strut;
        inputs.strut_angs = mp.P1.strut_angs;%Angles of the struts 
        inputs.strut_width = mp.P1.strut_width;% spider width (fraction of the pupil diameter)
        inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));% 

        mp.P1.full.mask = falco_gen_pupil_Simple( inputs );
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam;     % number of points across usable pupil   
        inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));%     % number of points across usable pupil 
        mp.P1.compact.mask = falco_gen_pupil_Simple( inputs );
        
    case{'WFIRST20180103'}
        %--Generate high-res input pupil for the 'full' model
        mp.P1.full.mask = falco_gen_pupil_WFIRST_20180103(mp.P1.full.Nbeam, mp.centering);
        
        %--Generate low-res input pupil for the 'compact' model
        mp.P1.compact.mask = falco_gen_pupil_WFIRST_20180103(mp.P1.compact.Nbeam, mp.centering);
        
    case{'WFIRST_onaxis'}
        inputs.strut_width = mp.pup_strut_width; %--0.0261  is nominal from 2014 on-axis (in pupil diameters)
        %--Generate input pupil
        inputs.Nbeam = mp.P1.full.Nbeam;     % number of points across usable pupil  
        inputs.Dbeam = mp.P2.D;
        %inputs.Narray = mp.P1.full.Narr;  % number of points across output array
        inputs.centering = mp.centering;
        inputs.clock_deg = 0; % clocking angle of the pupil (in degrees)
        inputs.magfacD = 1; %magnification factor of the pupil diameter
        inputs.xshift = 0;  % translation in x of pupil (in diameters)
        inputs.yshift = 0;  % translation in y of pupil (in diameters)
        %--Generate high-res input pupil for the 'full' model
        mp.P1.full.mask = falco_gen_pupil_WFIRSTcycle6_mag_rot_trans(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam;     % number of points across usable pupil    
        inputs.Dbeam = mp.P2.D;
        %inputs.Narray = mp.P1.compact.Narr;  % number of points across output array
        mp.P1.compact.mask = falco_gen_pupil_WFIRSTcycle6_mag_rot_trans(inputs);
        
    case{'LUVOIRA5'}
        inputs.centering = mp.centering;
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
        
    case{'LUVOIRA0'}
        
        
    case 'LUVOIR_B_offaxis'
        input.Nbeam = mp.P1.full.Nbeam/0.925; % number of points across the pupil diameter
        input.gapWidth = 6e-3/7.989*mp.P1.full.Nbeam; % samples
        input.numRings = 4;% Number of rings in hexagonally segmented mirror 
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.ID = 0; % central obscuration radius 
        input.OD = 1; % pupil outer diameter, can be < 1
        input.num_strut = 0;% Number of struts 
        input.strut_angs = [];%Angles of the struts (deg)
        input.strut_width = []; % Width of the struts (fraction of pupil diam.)

        missingSegments = ones(1,hexSegMirror_numSegments( input.numRings ));
        for index = 0:5
            missingSegments(38+index*4) = 0;
        end
        input.missingSegments = missingSegments;

        mp.P1.full.mask = falco_gen_pupil_customHex( input );
        
        
        input.Nbeam = mp.P1.compact.Nbeam/0.925; % number of points across the pupil diameter
        input.gapWidth = 6e-3/7.989*mp.P1.compact.Nbeam; % samples
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        mp.P1.compact.mask = falco_gen_pupil_customHex( input );
        
    case 'HabEx_B_offaxis'
        input.Nbeam = mp.P1.full.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.aperture_num = mp.P1.aperture_num;
        input.gap_size = mp.P1.gap_size;
        
        mp.P1.full.mask = falco_gen_pupil_HabEx_B(input);
        
        input.Nbeam = mp.P1.compact.Nbeam; % number of points across the pupil diameter
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        mp.P1.compact.mask = falco_gen_pupil_HabEx_B( input );
        
end
mp.P1.full.Narr = length(mp.P1.full.mask);  %--Total number of pixels across array containing the pupil in the full model. Add 2 pixels to Nbeam when the beam is pixel-centered.
mp.P1.compact.Narr = length(mp.P1.compact.mask);  %--Number of pixels across the array containing the input pupil in the compact model
mp.sumPupil = sum(sum(abs(mp.P1.full.mask).^2));
% if(mp.useGPU)
%     mp.P1.full.mask = gpuArray(mp.P1.full.mask);
% end

%--NORMALIZED (in pupil diameter) coordinate grids in the input pupil for making the tip/tilted input wavefront within the compact and full models
if(strcmpi(mp.centering,'interpixel') )
    mp.P2.full.xsDL = (- (mp.P1.full.Narr-1)/2:(mp.P1.full.Narr-1)/2)*mp.P2.full.dx/mp.P2.D;
    mp.P2.compact.xsDL = (-(mp.P1.compact.Narr-1)/2:(mp.P1.compact.Narr-1)/2)*mp.P2.compact.dx/mp.P2.D;
else
    mp.P2.full.xsDL = ( -mp.P1.full.Narr/2:(mp.P1.full.Narr/2 -1) )*mp.P2.full.dx/mp.P2.D;
    mp.P2.compact.xsDL = ( -mp.P1.compact.Narr/2:(mp.P1.compact.Narr/2 - 1) )*mp.P2.compact.dx/mp.P2.D;
end
[mp.P2.full.XsDL,mp.P2.full.YsDL] = meshgrid(mp.P2.full.xsDL);
[mp.P2.compact.XsDL,mp.P2.compact.YsDL] = meshgrid(mp.P2.compact.xsDL);

%--Interpolate the lower resolution input pupil for the compact model. 
%  Have to use interpolate instead of imresize() if the pupil is pixel-centered.
% % % mp.P1.compact.mask = interp2(mp.P2.full.XsDL,mp.P2.full.YsDL,mp.P1.full.mask,mp.P2.compact.XsDL,mp.P2.compact.YsDL,'linear',0);

%% DM aperture stops
mp.dm1.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm1.Dstop,mp.centering);
mp.dm1.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm1.Dstop,mp.centering);

mp.dm2.full.mask = falco_gen_DM_stop(mp.P2.full.dx,mp.dm2.Dstop,mp.centering);
mp.dm2.compact.mask = falco_gen_DM_stop(mp.P2.compact.dx,mp.dm2.Dstop,mp.centering);


%% Apodizer Masks for full and compact models

switch mp.coro
    case{'SPLC'}
        mp.P3.full.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.full.dx,mp.P2.D,mp.centering);
        mp.P3.full.Narr= length(mp.P3.full.mask);
        %mp.P3.compact.Nbeam = mp.P1.compact.Narr;  %--Number of pixels across the array containing the SP pupil in the compact model  
        
        mp.P3.full.dx = mp.P2.full.dx;
        mp.P3.compact.dx = mp.P2.compact.dx;

        %--Generate lower-resolution SP for the compact model
        mp.P3.compact.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.compact.dx,mp.P2.D,mp.centering);
        mp.P3.compact.Narr = length(mp.P3.compact.mask);
        
        %--Shaped Pupil Plane Coordinates (meters)
        if( strcmpi(mp.centering,'interpixel') )
            % mp.P3.full.xs = ( -(mp.P3.full.Narr-1)/2:(mp.P3.full.Narr-1)/2 ).'*mp.P3.full.dx;
            %mp.P3.compact.xs = ( -(mp.P3.compact.Narr-1)/2:(mp.P3.compact.Narr-1)/2 ).'*mp.P3.compact.dx;
        else
            % mp.P3.full.xs = (-mp.P3.full.Narr/2:(mp.P3.full.Narr/2-1)).'*mp.P3.full.dx;
            %mp.P3.compact.xs = (-mp.P3.compact.Narr/2:(mp.P3.compact.Narr/2-1)).'*mp.P3.compact.dx;
        end
        
        
%         %--Downsample the SP
%         %--Use the same values as for the regular input pupil. 
%         SPpad = padOrCropEven(mp.P3.full.mask,mp.P1.full.Narr); %--Need to pad the SP to get the grid sizes to match the input pupil
%         SPcompact = interp2(mp.P2.full.XsDL,mp.P2.full.YsDL,SPpad,mp.P2.compact.XsDL,mp.P2.compact.YsDL,'cubic',0);
%         
%         %--Crop down the low-resolution SP to get rid of extra zero padding. Speeds up the compact model.
%         SPsum = sum(SPcompact(:));
%         SPdiff = 0; counter = 2;
%         while( abs(SPdiff) <= 1e-7)
%             mp.P3.compact.Narr = length(SPcompact)-counter; %--Number of points across the cropped-down Lyot stop
%             SPdiff = SPsum - sum(sum( padOrCropEven(SPcompact, mp.P3.compact.Narr-2) )); %--Subtract an extra 2 to negate the extra step that overshoots.
%             counter = counter + 2;
%         end
%         mp.P3.compact.mask = padOrCropEven(SPcompact,mp.P3.compact.Narr); %--The cropped-down Lyot stop for the compact model       

    case{'Vortex','vortex','VC','AVC'}
        
        if(nnz(strcmp(mp.whichPupil,{'LUVOIR_B_offaxis','HabEx_B_offaxis'}))>0 && mp.flagApod)
            % Full aperture stop 
            mp.P3.full.Narr = 2^(nextpow2(mp.P1.full.Nbeam));
%             mp.P3.full.dx = mp.P2.full.dx;

            inputs.Nbeam = mp.P1.full.Nbeam;     % number of points across incoming beam 
            inputs.Npad = mp.P3.full.Narr;% 
            inputs.OD = mp.P3.ODnorm;
            inputs.ID = mp.P3.IDnorm;
            inputs.num_strut = 0;
            inputs.strut_angs = [];%Angles of the struts 
            inputs.strut_width = 0;% spider width (fraction of the pupil diameter)

            mp.P3.full.mask= falco_gen_pupil_Simple( inputs );

            % Compact aperture stop 
            inputs.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
            inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));% 
            
            mp.P3.compact.Narr = 2^(nextpow2(mp.P1.compact.Nbeam));
%            mp.P3.compact.dx = mp.P1.compact.dx;
            mp.P3.compact.mask = falco_gen_pupil_Simple( inputs );
        else
            disp('Using vortex without apodizer or aperture stop.')
        end
        
end

%% Lyot plane resolution, coordinates, and cropped-down mask for compact model
%--Resolution at Lyot Plane
mp.P4.full.dx = mp.P4.D/mp.P4.full.Nbeam;
mp.P4.compact.dx = mp.P4.D/mp.P4.compact.Nbeam; %mp.NlyotFactor*mp.P4.full.dx;

switch mp.whichPupil
    case{'Simple'}
        inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam 
        inputs.Npad = 2^(nextpow2(mp.P4.full.Nbeam));% 
        inputs.OD = mp.P4.ODnorm;
        inputs.ID = mp.P4.IDnorm;
        inputs.num_strut = mp.P4.num_strut;
        inputs.strut_angs = mp.P4.strut_angs;%Angles of the struts 
        inputs.strut_width = mp.P4.strut_width;% spider width (fraction of the pupil diameter)

        mp.P4.full.mask = falco_gen_pupil_Simple( inputs );
        
        inputs.Nbeam = mp.P4.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
        inputs.Npad = 2^(nextpow2(mp.P4.compact.Nbeam));% 
        
        mp.P4.compact.mask = falco_gen_pupil_Simple( inputs );
        
        
    case{'WFIRST_onaxis','WFIRST20180103'}
        %--Define Lyot stop generator function inputs for the 'full' optical model
        inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam  
        inputs.Dbeam = mp.P4.D; %--diameter of the beam at the mask (meters)
        %inputs.Narray = mp.Nlyot;   % number of points across output array
        inputs.ID = mp.P4.IDnorm;
        inputs.OD = mp.P4.ODnorm;
        inputs.strut_width = mp.LS_strut_width;
        inputs.centering = mp.centering;

        %--Make or read in Lyot stop (LS) for the 'full' model
        mp.P4.full.mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,'ROT180');

        %--Make or read in Lyot stop (LS) for the 'compact' model
        inputs.Nbeam = mp.P4.compact.Nbeam;     % number of points across incoming beam           
        %inputs.Narray = mp.NlyotCompact;   % number of points across output array
        mp.P4.compact.mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,'ROT180');

        
	case{'LUVOIRA5','LUVOIRA0'}
        
        %--Define Lyot stop generator function inputs for the 'full' optical model
        inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam  
        inputs.Dbeam = mp.P4.D; %--diameter of the beam at the mask (meters)
        %inputs.Narray = mp.Nlyot;   % number of points across output array
        inputs.ID = mp.P4.IDnorm;
        inputs.OD = mp.P4.ODnorm;
        inputs.strut_width = 0;
        inputs.centering = mp.centering;

        %--Make or read in Lyot stop (LS) for the 'full' model
        mp.P4.full.mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,'ROT180');

        %--Make or read in Lyot stop (LS) for the 'compact' model
        inputs.Nbeam = mp.P4.compact.Nbeam;     % number of points across incoming beam           
        %inputs.Narray = mp.NlyotCompact;   % number of points across output array
        mp.P4.compact.mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,'ROT180');
        
    case {'LUVOIR_B_offaxis','HabEx_B_offaxis'}
        inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam 
        inputs.Npad = 2^(nextpow2(mp.P4.full.Nbeam));% 
        inputs.OD = mp.P4.ODnorm;
        inputs.ID = mp.P4.IDnorm;
        inputs.num_strut = 0;
        inputs.strut_angs = [];%Angles of the struts 
        inputs.strut_width = 0;% spider width (fraction of the pupil diameter)

        mp.P4.full.mask = falco_gen_pupil_Simple( inputs );
        
        inputs.Nbeam = mp.P4.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
        inputs.Npad = 2^(nextpow2(mp.P4.compact.Nbeam));% 
        
        mp.P4.compact.mask = falco_gen_pupil_Simple( inputs );
       
        
end
if(mp.useGPU)
    mp.P4.full.mask = gpuArray(mp.P4.full.mask);
end

%--NORMALIZED Lyot plane coordinates (over the ENTIRE beam diameter, not cropped down to LS opening yet)
% if(strcmpi(mp.centering,'interpixel') )
%     xsLS = (-(mp.Nlyot-1)/2:(mp.Nlyot-1)/2)*mp.P4.full.dx/mp.P4.D;
%     xsLScompact = (-(mp.NlyotCompact-1)/2:(mp.NlyotCompact-1)/2)*mp.P4.compact.dx/mp.P4.D;
% else
%     xsLS = (-mp.Nlyot/2:(mp.Nlyot/2-1))*mp.P4.full.dx/mp.P4.D;
%     xsLScompact = (-mp.NlyotCompact/2:(mp.NlyotCompact/2-1))*mp.P4.compact.dx/mp.P4.D;
% end
% [XSls,YSls] = meshgrid(xsLS);
% [XSlsComp,YSlsComp] = meshgrid(xsLScompact);
% LScompact = interp2(XSls,YSls,mp.P4.full.mask,XSlsComp,YSlsComp,'linear',0);
%figure(202); imagesc(mp.P4.full.mask); axis xy equal tight;
%figure(203); imagesc(mp.P4.full.mask-padOrCropEven(mp.P4.full.croppedMask,mp.Nlyot)); axis xy equal tight; colorbar;




%% Crop down the Lyot stop(s) to get rid of extra zero padding for the full model
switch mp.coro
    case{'Vortex','vortex','AVC','VC'}
        mp.P4.full.Narr = length(mp.P4.full.mask);
        mp.P4.full.croppedMask = mp.P4.full.mask;
        mp.P4.compact.Narr = length(mp.P4.compact.mask);
        mp.P4.compact.croppedMask = mp.P4.compact.mask;
    otherwise
    
        % --Crop down the high-resolution Lyot stop to get rid of extra zero padding
        LSsum = sum(mp.P4.full.mask(:));
        LSdiff = 0; counter = 2;
        while( abs(LSdiff) <= 1e-7)
            mp.P4.full.Narr = length(mp.P4.full.mask)-counter;
            LSdiff = LSsum - sum(sum( padOrCropEven(mp.P4.full.mask, mp.P4.full.Narr-2) )); %--Subtract an extra 2 to negate the extra step that overshoots.
            counter = counter + 2;
        end
        mp.P4.full.croppedMask = padOrCropEven(mp.P4.full.mask,mp.P4.full.Narr); %--The cropped-down Lyot stop for the full model.


        % --Crop down the low-resolution Lyot stop to get rid of extra zero padding. Speeds up the compact model.
        LSsum = sum(mp.P4.compact.mask(:));
        LSdiff = 0; counter = 2;
        while( abs(LSdiff) <= 1e-7)
            mp.P4.compact.Narr = length(mp.P4.compact.mask)-counter; %--Number of points across the cropped-down Lyot stop
            LSdiff = LSsum - sum(sum( padOrCropEven(mp.P4.compact.mask, mp.P4.compact.Narr-2) )); %--Subtract an extra 2 to negate the extra step that overshoots.
            counter = counter + 2;
        end
        mp.P4.compact.croppedMask = padOrCropEven(mp.P4.compact.mask,mp.P4.compact.Narr); %--The cropped-down Lyot stop for the compact model       




end


%--(METERS) Lyot plane coordinates (over the cropped down to Lyot stop mask) for MFTs in the compact model from the FPM to the LS.
if(strcmpi(mp.centering,'interpixel') )
    mp.P4.compact.xs = (-(mp.P4.compact.Narr-1)/2:(mp.P4.compact.Narr-1)/2)*mp.P4.compact.dx;
else
    mp.P4.compact.xs = (-mp.P4.compact.Narr/2:(mp.P4.compact.Narr/2-1))*mp.P4.compact.dx;
end
mp.ysLScompactCrop = mp.P4.compact.xs.';




%% Generate/Load FPM

switch mp.coro
    case{'Vortex','vortex','AVC','VC'}
        %--Vortex FPM is generated as needed
    otherwise
        %--Make or read in focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.IWA = mp.F3.Rin;% inner radius of the focal plane mask, in lambda0/D
        FPMgenInputs.OWAmask = mp.F3.Rout; % outer radius of the focal plane mask, in lambda0/D
        %FPMgenInputs.flagOdd = false; % flag to specify odd or even-sized array
        FPMgenInputs.centering = mp.centering;
        FPMgenInputs.ang = mp.F3.ang; % angular opening on each side of the focal plane mask, in degrees
        FPMgenInputs.magx = 1; % magnification factor along the x-axis
        FPMgenInputs.magy = 1; % magnification factor along the y-axis
end


switch mp.coro
    
    case {'LC','DMLC','APLC'} %--Occulting spot FPM
        
        
        %--Make or read in focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = mp.FPMampFac; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
        mp.F3.full.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        % figure(204); imagesc(mp.F3.full.mask.amp); axis xy equal tight;

        mp.F3.full.Nxi = size(mp.F3.full.mask.amp,2);
        mp.F3.full.Neta= size(mp.F3.full.mask.amp,1);   
        
        %--Number of points across the FPM in the compact model
        if(mp.F3.Rout==inf)
            switch mp.centering
            case 'pixel'
                mp.F3.compact.Nxi = ceil_even((2*(mp.F3.Rin*mp.F3.compact.res + 1/2)));
            case 'interpixel'
                mp.F3.compact.Nxi = ceil_even((2*mp.F3.Rin*mp.F3.compact.res));
            end
        else
            switch mp.centering
                case 'pixel'
                    mp.F3.compact.Nxi = ceil_even((2*(mp.F3.Rout*mp.F3.compact.res + 1/2)));
                case 'interpixel'
                    mp.F3.compact.Nxi = ceil_even((2*mp.F3.Rout*mp.F3.compact.res));
            end
        end
        mp.F3.compact.Neta = mp.F3.compact.Nxi;
        
        %--Coordinates for the FPMs in the full and compact models
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.Nxi,2)==1  )
            mp.F3.full.xisDL  = (-(mp.F3.full.Nxi -1)/2:(mp.F3.full.Nxi -1)/2)/mp.F3.full.res;
            mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
            
            mp.F3.compact.xisDL  = (-(mp.F3.compact.Nxi -1)/2:(mp.F3.compact.Nxi -1)/2)/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2)/mp.F3.compact.res;
        else
            mp.F3.full.xisDL  = (-mp.F3.full.Nxi/2:(mp.F3.full.Nxi/2-1))/mp.F3.full.res;
            mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
            
            mp.F3.compact.xisDL  = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1))/mp.F3.compact.res;
        end
        
        
        %--Make or read in focal plane mask (FPM) amplitude for the compact model
        FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
        mp.F3.compact.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        
%         %--Downsample the FPM for the compact model --> DO NOT USE: Does not handle the true width of the mask correctly
%         [XIS,ETAS] = meshgrid(mp.F3.full.xisDL,mp.F3.full.etasDL);
%         [XIScompact,ETAScompact] = meshgrid(mp.F3.compact.xisDL,mp.F3.compact.etasDL);
%         if(mp.F3.Rout==inf); extrapval = 1; else extrapval = 0; end
%         mp.F3.compact.mask.amp = interp2(XIS,ETAS,mp.F3.full.mask.amp,XIScompact,ETAScompact,'cubic',extrapval);
%         % figure(205); imagesc(mp.F3.compact.mask.amp); axis xy equal tight;
         
        
    case{'SPLC'}
        
        %--Generate the focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        FPMgenInputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
        FPMgenInputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
        FPMgenInputs.FPMampFac = mp.FPMampFac; % amplitude transmission of inner FPM spot
        FPMgenInputs.centering = mp.centering;
        mp.F3.full.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        %figure(204); imagesc(mp.F3.full.mask.amp); axis xy equal tight; drawnow;

        mp.F3.full.Nxi = size(mp.F3.full.mask.amp,2);
        mp.F3.full.Neta= size(mp.F3.full.mask.amp,1);   
        
        %--Number of points across the FPM in the compact model
        switch mp.centering
            case 'pixel'
                mp.F3.compact.Nxi = ceil_even(2*(mp.F3.Rout*mp.F3.compact.res + 0.5));
            case 'interpixel'
                mp.F3.compact.Nxi = ceil_even(2*mp.F3.Rout*mp.F3.compact.res);
        end
        mp.F3.compact.Neta = mp.F3.compact.Nxi;
        
        %--Coordinates for the FPMs in the full and compact models
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.Nxi,2)==1  )
            mp.F3.full.xisDL  = (-(mp.F3.full.Nxi -1)/2:(mp.F3.full.Nxi -1)/2)/mp.F3.full.res;
            mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
            
            mp.F3.compact.xisDL  = (-(mp.F3.compact.Nxi -1)/2:(mp.F3.compact.Nxi -1)/2)/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2)/mp.F3.compact.res;
        else
            mp.F3.full.xisDL  = (-mp.F3.full.Nxi/2:(mp.F3.full.Nxi/2-1))/mp.F3.full.res;
            mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
            
            mp.F3.compact.xisDL  = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))/mp.F3.compact.res;
            mp.F3.compact.etasDL = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1))/mp.F3.compact.res;
        end
        
        %--Generate the FPM amplitude for the compact model
        FPMgenInputs.pixresFPM = mp.F3.compact.res; %--pixels per lambda_c/D
        mp.F3.compact.mask.amp = falco_gen_annular_FPM(FPMgenInputs);
        %figure(205); imagesc(mp.F3.compact.mask.amp); axis xy equal tight; drawnow;
        
%         %--DOWNSAMPLING MIGHT BE CAUSING BIG MODEL DISCREPANCIES:
%         %--Downsample the FPM for the compact model
%         [XIS,ETAS] = meshgrid(mp.F3.full.xisDL,mp.F3.full.etasDL);
%         [XIScompact,ETAScompact] = meshgrid(mp.F3.compact.xisDL,mp.F3.compact.etasDL);
%         mp.F3.compact.mask.amp = interp2(XIS,ETAS,mp.F3.full.mask.amp,XIScompact,ETAScompact,'cubic',0);
%         %figure(205); imagesc(mp.F3.compact.mask.amp); axis xy equal tight;
        


end

%% FPM coordinates

switch mp.coro

    case{'Vortex','vortex','AVC','VC'}
        %--Nothing needed to run the vortex model
        
    otherwise
        mp.F3.full.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.full.res;
        mp.F3.full.deta = mp.F3.full.dxi;
        mp.F3.compact.dxi = (mp.fl*mp.lambda0/mp.P2.D)/mp.F3.compact.res;
        mp.F3.compact.deta = mp.F3.compact.dxi;
        %--Compute coordinates in plane of FPM in the compact model (in meters)
        if(strcmpi(mp.centering,'interpixel') || mod(mp.F3.compact.Nxi,2)==1  )
            %mp.xisFPMcompact = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)*mp.F3.compact.dxi;
            %mp.etasFPMcompact = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2).'*mp.F3.compact.deta;
            
            mp.F3.compact.xis  = (-(mp.F3.compact.Nxi -1)/2:(mp.F3.compact.Nxi -1)/2)*mp.F3.compact.dxi;
            mp.F3.compact.etas = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2).'*mp.F3.compact.deta;
        else
            %mp.xisFPMcompact = (-mp.F3.compact.Nxi/2:(mp.F3.compact.Nxi/2-1))*mp.F3.compact.dxi;
            %mp.etasFPMcompact = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1)).'*mp.F3.compact.deta;
            
            mp.F3.compact.xis  = (-mp.F3.compact.Nxi/2: (mp.F3.compact.Nxi/2 -1))*mp.F3.compact.dxi;
            mp.F3.compact.etas = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1)).'*mp.F3.compact.deta;
        end
end

%% Sampling/Resolution and Scoring/Correction Masks for 2nd Focal Plane

mp.F4.compact.dxi = mp.fl*mp.lambda0/mp.P4.D/mp.F4.compact.res;
mp.F4.compact.deta = mp.F4.compact.dxi;       

mp.F4.full.dxi = mp.fl*mp.lambda0/mp.P4.D/mp.F4.full.res;
mp.F4.full.deta = mp.F4.full.dxi;    

% mp.P1.full.Narr = length(mp.P1.full.mask); %--Number of points across at entrance pupil
% mp.Ndm1 = mp.P1.full.Narr; %--May want low-res versions as well
% mp.Ndm2 = mp.P1.full.Narr; %--May want low-res versions as well

%--Make the scoring and correction masks for the camera's focal plane
% maskFP2.centering = mp.centering; % pixel centering
% maskFP2.whichSide = mp.F4.sides; %--options: 'left', 'right', or 'both'
% maskFP2.pixres = mp.F4.compact.res; % pixels per lambda0/D at the camera
% maskFP2.FOV = mp.F4.FOV; % minimum desired field of view (along xi axis) in lambda0/D

%--Software Mask for Correction 
%--Set Inputs
maskCorr.pixresFP = mp.F4.compact.res;
maskCorr.rhoInner = mp.F4.corr.Rin; %--lambda0/D
maskCorr.rhoOuter = mp.F4.corr.Rout ; %--lambda0/D
maskCorr.angDeg = mp.F4.corr.ang; %--degrees
maskCorr.centering = mp.centering;
maskCorr.FOV = mp.F4.FOV;
maskCorr.whichSide = mp.F4.sides; %--which (sides) of the dark hole have open
%--Compact Model: Generate Software Mask for Correction 
[mp.F4.compact.corr.mask, mp.F4.compact.xisDL, mp.F4.compact.etasDL] = falco_gen_SW_mask(maskCorr); 
mp.F4.compact.corr.settings = maskCorr; %--Store values for future reference
%--Full Model: Generate Software Mask for Correction 
maskCorr.pixresFP = mp.F4.full.res;
[mp.F4.full.corr.mask, mp.F4.full.xisDL, mp.F4.full.etasDL] = falco_gen_SW_mask(maskCorr); 
mp.F4.full.corr.settings = maskCorr; %--Store values for future reference



%--Software Mask for Scoring Contrast 
%--Set Inputs
maskScore.pixresFP = mp.F4.compact.res;
maskScore.rhoInner = mp.F4.score.Rin; %--lambda0/D
maskScore.rhoOuter = mp.F4.score.Rout ; %--lambda0/D
maskScore.angDeg = mp.F4.score.ang; %--degrees
maskScore.centering = mp.centering;
maskScore.FOV = mp.F4.FOV;
maskScore.whichSide = mp.F4.sides; %--which (sides) of the dark hole have open
%--Compact Model: Generate Software Mask for Scoring Contrast 
[mp.F4.compact.score.mask,~,~] = falco_gen_SW_mask(maskScore); 
mp.F4.compact.score.settings = maskScore; %--Store values for future reference
%--Full Model: Generate Software Mask for Scoring Contrast 
maskScore.pixresFP = mp.F4.full.res;
[mp.F4.full.score.mask,~,~] = falco_gen_SW_mask(maskScore); 
mp.F4.full.score.settings = maskScore; %--Store values for future reference

%--One-sided scoring and correction masks
% maskCMr = maskCorr; maskCMr.whichSide = 'right'; [mp.MaskCorrRight,~,~] = falco_gen_SW_mask(maskCMr);
% maskCMl = maskCorr; maskCMl.whichSide = 'left';  [mp.MaskCorrLeft,~,~]  = falco_gen_SW_mask(maskCMl);
% maskSMr = maskScore; maskSMr.whichSide = 'right'; [mp.MaskScoreRight,~,~] = falco_gen_SW_mask(maskSMr);
% maskSMl = maskScore; maskSMl.whichSide = 'left';  [mp.MaskScoreLeft,~,~]  = falco_gen_SW_mask(maskSMl);

%--Indices of dark hole pixels
mp.F4.compact.corr.inds = find(mp.F4.compact.corr.mask~=0);
mp.F4.compact.score.inds = find(mp.F4.compact.score.mask~=0);

mp.F4.full.score.inds = find(mp.F4.full.score.mask~=0);
mp.F4.full.corr.inds = find(mp.F4.full.corr.mask~=0);

% mp.cor_ele_right = find(mp.MaskCorrRight~=0);
% mp.cor_ele_left = find(mp.MaskCorrLeft~=0);
% mp.score_ele_right = find(mp.MaskScoreRight~=0);
% mp.score_ele_left = find(mp.MaskScoreLeft~=0);

mp.F4.compact.Nxi  = size(mp.F4.compact.score.mask,2);
mp.F4.compact.Neta = size(mp.F4.compact.score.mask,1); 

mp.F4.full.Nxi  = size(mp.F4.full.score.mask,2);
mp.F4.full.Neta = size(mp.F4.full.score.mask,1); 


%% --Spatial weighting of pixel intensity (compac model only since for control)
[XISLAMD,ETASLAMD] = meshgrid(mp.F4.compact.xisDL, mp.F4.compact.etasDL);
RHOS = sqrt(XISLAMD.^2+ETASLAMD.^2);

mp.Wspatial = mp.F4.compact.corr.mask;
if( isfield(mp,'WspatialDef') ) %--If there are spatial weights for the Jacobian
    if(isempty(mp.WspatialDef)==false)
        mp.WspatialDef = mp.WspatialDef;
    
        for kk=1:size(mp.WspatialDef,1)
            Wannulus = 1+(sqrt(mp.WspatialDef(kk,3))-1)*( (RHOS>=mp.WspatialDef(kk,1)) & (RHOS<mp.WspatialDef(kk,2)) );
            %figure(70);imagesc(mp.F4.compact.xisDL,mp.F4.compact.etasDL,Wannulus); axis xy equal tight; colorbar;
            mp.Wspatial = mp.Wspatial.*Wannulus;
            %figure(71);imagesc(mp.F4.compact.xisDL,mp.F4.compact.etasDL,mp.Wspatial); axis xy equal tight; colorbar; pause(1);
        end
    else
        %--nothing
    end
    
else
 %--nothing
end
mp.Wspatial_ele = mp.Wspatial(mp.F4.compact.corr.inds);


%% Bandpass parameters                                                            
if mp.Nsbp > 1
%     Set so that the edge of the extremum sub-bandpasses match the edge of the entire bandpass.
    sbp_center_vec = linspace(mp.lambda0*(1-(mp.fracBW-mp.fracBWsbp)/2),...
                          mp.lambda0*(1+(mp.fracBW-mp.fracBWsbp)/2),...
                          mp.Nsbp); % vector of central wavelengths for each sub-bandpass
else
    sbp_center_vec = [mp.lambda0];
end
mp.sbp_center_vec = sbp_center_vec;

if mp.Nwpsbp > 1 && mp.Nsbp==1 % Fill the whole band if just one sub-bandpass
    mp.lamFac_vec = linspace(1-mp.fracBWsbp/2,...
                          1+mp.fracBWsbp/2,...
                          mp.Nwpsbp); % vector of lambda factors about the center wavelength of each sub-bandpass
                      
elseif mp.Nwpsbp > 1 
    % Evenly spaced so that in full bandpass no wavelength is weighted more heavily.
    mp.lamFac_vec = linspace( 1-(mp.fracBWsbp/2)*(1-1/mp.Nwpsbp),...
                           1+(mp.fracBWsbp/2)*(1-1/mp.Nwpsbp),...
                          mp.Nwpsbp);
else
    mp.lamFac_vec = 1;
end
mp.lam_array = mp.lambda0*mp.lamFac_vec;
% fprintf('Wavelengths cover %3.1f%% for the %3.1f%% sub-bandpass (%3.1f%% coverage). \n', (max(lamFac_vec)-min(lamFac_vec))*100 , mp.fracBWsbp*100,(max(lamFac_vec)-min(lamFac_vec))/mp.fracBWsbp*100);

%% Deformable Mirror (DM) Parameters
DM.dm1.NactTotal=0; DM.dm2.NactTotal=0; DM.dm3.NactTotal=0; DM.dm4.NactTotal=0; DM.dm5.NactTotal=0; DM.dm6.NactTotal=0; DM.dm7.NactTotal=0; DM.dm8.NactTotal=0; DM.dm9.NactTotal=0; %--Initialize for bookkeeping later.

%--Centering of DM surfaces on array
DM.dm1.centering = mp.centering;
DM.dm2.centering = mp.centering;


%--Create influence function datacubes for each DM
if( any(DM.dm_ind==1) ) %if(isfield(DM.dm1,'inf_datacube')==0 && any(DM.dm_ind==1) )
    DM.dm1.compact = DM.dm1;
%     DM.dm1.dx = mp.P2.full.dx;
%     DM.dm1.compact.dx = mp.P2.compact.dx;
    DM.dm1 = falco_gen_dm_poke_cube(DM.dm1, mp, mp.P2.full.dx,'NOCUBE');
    DM.dm1.compact = falco_gen_dm_poke_cube(DM.dm1.compact, mp, mp.P2.compact.dx);
%     DM.dm1 = falco_gen_dm_poke_cube_PROPER(DM.dm1,mp,'NOCUBE');
%     DM.dm1.compact = falco_gen_dm_poke_cube_PROPER(DM.dm1.compact,mp);
end
if( any(DM.dm_ind==2) ) %if(isfield(DM.dm2,'inf_datacube')==0 && any(DM.dm_ind==2) )
    DM.dm2.compact = DM.dm2;
    DM.dm2.dx = mp.P2.full.dx;
    DM.dm2.compact.dx = mp.P2.compact.dx;
    
    DM.dm2 = falco_gen_dm_poke_cube(DM.dm2, mp, mp.P2.full.dx, 'NOCUBE');
    DM.dm2.compact = falco_gen_dm_poke_cube(DM.dm2.compact, mp, mp.P2.compact.dx);
%     DM.dm2 = falco_gen_dm_poke_cube_PROPER(DM.dm2,mp,'NOCUBE');
%     DM.dm2.compact = falco_gen_dm_poke_cube_PROPER(DM.dm2.compact,mp);
end



% % DM.dm1.Ndm = min([ceil_even(mp.P1.full.Narr*(DM.dm1.dm_spacing*DM.dm1.Nact)/(mp.Ddm1)), DM.dm1.NdmPad]); %--Number of points across to crop the surface to for storage
% % DM.dm2.Ndm = min([ceil_even(mp.P1.full.Narr*(DM.dm2.dm_spacing*DM.dm2.Nact)/(mp.Ddm1)), DM.dm2.NdmPad]); %--Number of points across to crop the surface to for storage
% % DM.dm1.compact.Ndm = min([ceil_even(mp.P1.compact.Narr*(DM.dm1.compact.dm_spacing*DM.dm1.compact.Nact)/(mp.Ddm1)), DM.dm1.compact.NdmPad]); %--Number of points across to crop the surface to for storage
% % DM.dm2.compact.Ndm = min([ceil_even(mp.P1.compact.Narr*(DM.dm2.compact.dm_spacing*DM.dm1.compact.Nact)/(mp.Ddm1)), DM.dm2.compact.NdmPad]); %--Number of points across to crop the surface to for storage

%--Initial DM voltages
DM.dm1.V = zeros(DM.dm1.Nact);
DM.dm2.V = zeros(DM.dm2.Nact);

%--Peak-to-Valley DM voltages
DM.dm1.Vpv = zeros(mp.Nitr,1);
DM.dm2.Vpv = zeros(mp.Nitr,1);

%--First delta DM settings are zero (for covariance calculation in Kalman filters or robust controllers)
DM.dm1.dV = zeros(DM.dm1.Nact,DM.dm1.Nact);  % delta voltage on DM1;
DM.dm2.dV = zeros(DM.dm2.Nact,DM.dm2.Nact);  % delta voltage on DM2;

%--Store the DM commands 
DM.dm1.Vall = zeros(DM.dm1.Nact,DM.dm1.Nact,mp.Nitr+1);
DM.dm2.Vall = zeros(DM.dm2.Nact,DM.dm2.Nact,mp.Nitr+1);

%--Initialize the number of actuators used.
DM.dm1.Nele=[]; DM.dm2.Nele=[];  DM.dm3.Nele=[];  DM.dm4.Nele=[];  DM.dm5.Nele=[];  DM.dm6.Nele=[];  DM.dm7.Nele=[];  DM.dm8.Nele=[];  DM.dm9.Nele=[]; %--Initialize for Jacobian calculations later. 
if(any(DM.dm_ind==1)); DM.dm1.Nele = length(DM.dm1.act_ele); end
if(any(DM.dm_ind==2)); DM.dm2.Nele = length(DM.dm2.act_ele); end
if(any(DM.dm_ind==9)); DM.dm9.Nele = length(DM.dm9.act_ele); end
DM.NelePerDMvec = [length(DM.dm1.Nele), length(DM.dm2.Nele), length(DM.dm3.Nele), length(DM.dm4.Nele), length(DM.dm5.Nele), length(DM.dm6.Nele), length(DM.dm7.Nele), length(DM.dm8.Nele), length(DM.dm9.Nele) ];
DM.NactTotals = [DM.dm1.NactTotal, DM.dm2.NactTotal, DM.dm3.NactTotal, DM.dm4.NactTotal, DM.dm5.NactTotal, DM.dm6.NactTotal, DM.dm7.NactTotal, DM.dm8.NactTotal, DM.dm9.NactTotal]; 
DM.NactTotal = sum(DM.NactTotals);

%% Array Sizes for Angular Spectrum Propagation with FFTs

%--Compact Model: Set nominal DM plane array sizes as a power of 2 for angular spectrum propagation with FFTs
if( any(DM.dm_ind==1) && any(DM.dm_ind==2) )
    NdmPad = 2.^ceil(1 + log2(max([DM.dm1.compact.NdmPad,DM.dm2.compact.NdmPad]))); 
elseif(  any(DM.dm_ind==1) )
    NdmPad = 2.^ceil(1 + log2(DM.dm1.compact.NdmPad));
elseif(  any(DM.dm_ind==2) )
    NdmPad = 2.^ceil(1 + log2(DM.dm2.compact.NdmPad));
end
while( (NdmPad < min(mp.lam_array)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) || (NdmPad < min(mp.lam_array)*abs(mp.d_P2_dm1)/mp.P2.compact.dx^2) ) %--Double the zero-padding until the angular spectrum sampling requirement is not violated
    NdmPad = 2*NdmPad; 
end
DM.compact.NdmPad = NdmPad;

%--Full Model: Set nominal DM plane array sizes as a power of 2 for angular spectrum propagation with FFTs
if( any(DM.dm_ind==1) && any(DM.dm_ind==2) )
    NdmPad = 2.^ceil(1 + log2(max([DM.dm1.NdmPad,DM.dm1.NdmPad]))); 
elseif(  any(DM.dm_ind==1) )
    NdmPad = 2.^ceil(1 + log2(DM.dm1.NdmPad));
elseif(  any(DM.dm_ind==2) )
    NdmPad = 2.^ceil(1 + log2(DM.dm2.NdmPad));
end
while( (NdmPad < min(mp.lam_array)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) || (NdmPad < min(mp.lam_array)*abs(mp.d_P2_dm1)/mp.P2.full.dx^2) ) %--Double the zero-padding until the angular spectrum sampling requirement is not violated
    NdmPad = 2*NdmPad; 
end
DM.full.NdmPad = NdmPad;

% %--For propagation from pupil P2 to DM1.
% mp.P2.compact.Nfft =  2.^ceil(1 + log2(mp.P1.compact.Nbeam));
% while( (mp.P2.compact.Nfft < min(mp.lam_array)*abs(mp.d_dm1_dm2)/mp.P2.compact.dx^2) ||  (mp.P2.compact.Nfft < min(mp.lam_array)*abs(mp.d_P2_dm1)/mp.P2.compact.dx^2) ); 
%     mp.P2.compact.Nfft = 2*mp.P2.compact.Nfft; %--Double the zero-padding until the angular spectrum sampling requirement is not violated
% end 
% 
% mp.P2.full.Nfft =  2.^ceil(1 + log2(mp.P1.full.Nbeam));
% while( (mp.P2.full.Nfft < min(mp.lam_array)*abs(mp.d_dm1_dm2)/mp.P2.full.dx^2) ||  (mp.P2.full.Nfft < min(mp.lam_array)*abs(mp.d_P2_dm1)/mp.P2.full.dx^2) ); 
%     mp.P2.full.Nfft = 2*mp.P2.full.Nfft; %--Double the zero-padding until the angular spectrum sampling requirement is not violated
% end 


%% % % Initial Electric Fields for Star and Exoplanet
% Starlight. Can add some propagation here to create an aberrate wavefront
%   starting from a primary mirror.
mp.P1.full.E  = ones(mp.P1.full.Narr,mp.P1.full.Narr,mp.Nwpsbp,mp.Nsbp); % Input E-field at entrance pupil
% mp.P1.full.E = mp.Estar; %--Re-define in the main code if pupil phase flattening is done.
mp.Eplanet = mp.P1.full.E; %--Initialize the input E-field for the planet at the entrance pupil. Will apply the phase ramp later
mp.P1.compact.E = ones(mp.P1.compact.Narr,mp.P1.compact.Narr,mp.Nwpsbp,mp.Nsbp);

%% Off-axis, incoherent point source (exoplanet)
mp.c_planet = 1;%1e-14;%4e-10;%3e-10;%1e-8; % contrast of exoplanet
mp.x_planet = 6;%4; % x position of exoplanet in lambda0/D
mp.y_planet = 0;%1/2; % 7 position of exoplanet in lambda0/D

%% Field Stop

mp.F4.compact.mask = ones(mp.F4.compact.Neta,mp.F4.compact.Nxi);
% figure; imagesc(Lam0Dxi,Lam0Deta,FPMstop); axis xy equal tight; colormap gray; colorbar; title('Field Stop');



%% Initialize Storage Arrays for Controller
switch mp.controller
    case{'EFC'} % 'EFC' = empirical grid search over both overall scaling coefficient and Lagrange multiplier
        % Take images for different Lagrange multiplier values and overall command gains and pick the value pair that gives the best contrast
        cp.contrast_array_mus_meas = zeros(length(cp.muVec),mp.Nitr);
        cp.contrast_array_mus_linMod = zeros(length(cp.muVec),mp.Nitr);
        cp.PtoVdu1V = zeros(length(cp.muVec),mp.Nitr);
        cp.PtoVdu2V = zeros(length(cp.muVec),mp.Nitr);
        cp.muBest = zeros(mp.Nitr,1);
        cp.dmfacBest = zeros(mp.Nitr,1);
    case{'conEFC'}
        cvx_startup
end

%% Contrast to Normalized Intensity Map Calculation 

%%
fprintf('Beginning Trial %d of Series %d.\n',mp.TrialNum,mp.SeriesNum);

%% Get the starlight normalization factor for the compact and full models (to convert images to normalized intensity)
mp = falco_get_PSF_norm_factor(mp, DM);

% %% Get normalization factors for the models
% mp = falco_get_PSF_norm_factor(mp, DM);


%% Check that the normalization of the coronagraphic PSF is correct

modvar.ttIndex = 1;
modvar.flagCalcJac = 0; 
modvar.sbpIndex = mp.si_ref;
modvar.wpsbpIndex = mp.wi_ref;
modvar.whichSource = 'star';     

E0c = model_compact(mp, DM, modvar);
I0c = abs(E0c).^2;
if(mp.flagPlot)
    figure(501); imagesc(log10(I0c)); axis xy equal tight; colorbar;
    title('(Compact Model: Normalization Check Using Starting PSF)'); 
    drawnow;
end
E0f = model_full(mp, DM, modvar);
I0f = abs(E0f).^2;
if(mp.flagPlot)
    figure(502); imagesc(log10(I0f)); axis xy equal tight; colorbar;
    title('(Full Model: Normalization Check Using Starting PSF)'); drawnow;
end

%% Coordinates and scoring mask for core throughput calculation (full model)

[XIS,ETAS] = meshgrid(mp.F4.full.xisDL - mp.x_planet, mp.F4.full.etasDL - mp.y_planet);
mp.FP4.compact.RHOS = sqrt(XIS.^2 + ETAS.^2);
mp.maskHMcore = 0*mp.FP4.compact.RHOS;
mp.maskCore  = 0*mp.FP4.compact.RHOS;
mp.maskCore(mp.FP4.compact.RHOS<=mp.thput_radius) = 1;

mp.thput_vec = zeros(mp.Nitr,1);


end %--END OF FUNCTION
