% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_gen_chosen_pupil(mp)
%
% Function to generate the apodizer representation based on configuration settings.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-05-29 by A.J. Riggs.




function mp = falco_config_gen_chosen_pupil(mp)

%% Input pupil plane resolution, masks, and coordinates
%--Resolution at input pupil and DM1 and DM2
mp.P2.full.dx = mp.P2.D/mp.P1.full.Nbeam; 
mp.P2.compact.dx = mp.P2.D/mp.P1.compact.Nbeam;
% mp.Npup1factor = mp.P1.full.Nbeam/mp.P1.compact.Nbeam; %--Scaling factor between full model and compact model of the input pupil array
% % % mp.Npup1factor = mp.P1.full.Narr/mp.P1.compact.Narr; %--Scaling factor between full model and compact model of the input pupil array
% % % mp.P2.compact.dx = mp.Npup1factor*mp.P2.full.dx;


%--Generate/Load Input Pupil
switch mp.whichPupil
    case{'Simple','SimplePROPER'}
        
        if(strcmpi(mp.whichPupil,'simpleproper'))
            inputs.flagPROPER = true;
        end
        
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

    case{'WFIRST180718'}
        %--Generate high-res input pupil for the 'full' model
        mp.P1.full.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P1.full.Nbeam, mp.centering);
        
        %--Generate low-res input pupil for the 'compact' model
        mp.P1.compact.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P1.compact.Nbeam, mp.centering);        
        
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
        if(isfield(mp.P1,'strut_width'))
            inputs.strut_width = mp.P1.strut_width;% spider width (fraction of the pupil diameter)
        end
        
        if(isfield(mp.P1,'gap_width_m'))
            inputs.gap_width_m = mp.P1.gap_width_m;% spider width (fraction of the pupil diameter)
        end
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
        
        
%         %--TESTING TESTING TESTING: Trying downsampling pupil for full model
%         switch mp.centering
%             case 'interpixel'
%                 inputs.Nbeam = 5000;
%                 mask_temp = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
%                 mp.P1.full.mask = imresize(mask_temp, mp.P1.full.Nbeam/inputs.Nbeam, 'cubic');
%             case 'pixel'
%                 Nstep2 = 2*mp.P1.full.Nbeam;
%                 Nstep1 = 5000;%floor(5000/Nstep2)*Nstep2;
%                 inputs.centering = 'interpixel'; %--temporarily switch
%                 inputs.Nbeam = Nstep1;
%                 mask_temp = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
%                 mask_temp = imresize(mask_temp,[Nstep2,Nstep2],'bicubic');
%                 mask_temp = padOrCropEven(mask_temp, 2*(mp.P1.full.Nbeam+1) );
%                 mask_temp = imresize(mask_temp,1/2,'bicubic');
%                 mp.P1.full.mask = zeros(mp.P1.full.Nbeam+2);
%                 mp.P1.full.mask(2:end,2:end) = mask_temp;   
%         end
        
        
        
%         %--TESTING TESTING TESTING: Trying downsampling pupil for compact model
%         switch mp.centering
%             case 'interpixel'
%                 mp.P1.compact.mask = imresize(mp.P1.full.mask, mp.P1.compact.Nbeam/mp.P1.full.Nbeam, 'cubic');
%             case 'pixel'
%                 Nstep2 = 2*mp.P1.compact.Nbeam;
%                 Nstep1 = 5000;%floor(5000/Nstep2)*Nstep2;
%                 inputs.centering = 'interpixel'; %--temporarily switch
%                 inputs.Nbeam = Nstep1;
%                 mask_temp = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
%                 mask_temp = imresize(mask_temp,[Nstep2,Nstep2],'bicubic');
%                 mask_temp = padOrCropEven(mask_temp, 2*(mp.P1.compact.Nbeam+1) );
%                 mask_temp = imresize(mask_temp,1/2,'bicubic');
%                 mp.P1.compact.mask = zeros(mp.P1.compact.Nbeam+2);
%                 mp.P1.compact.mask(2:end,2:end) = mask_temp;   
%         end
        
        clear inputs
        
%         %--DEBUGGING: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         mp.P1.full.mask = ones(size(mp.P1.full.mask));
%         mp.P1.compact.mask = ones(size(mp.P1.compact.mask));
        
    case{'LUVOIRA0'}
        inputs.centering = mp.centering;
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_0(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_0(inputs);
        
%         %--DEBUGGING: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         mp.P1.full.mask = ones(size(mp.P1.full.mask));
%         mp.P1.compact.mask = ones(size(mp.P1.compact.mask));
        
        
    case 'LUVOIR_B_offaxis'
        input.Nbeam = mp.P1.full.Nbeam/0.925; % number of points across the pupil diameter
        input.gapWidth = mp.P1.gapWidth*mp.P1.full.Nbeam; % samples
        %input.gapWidth = 6e-3/7.989*mp.P1.full.Nbeam; % samples
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
        %input.gapWidth = 6e-3/7.989*mp.P1.compact.Nbeam; % samples
        input.gapWidth = mp.P1.gapWidth*mp.P1.compact.Nbeam; % samples
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        mp.P1.compact.mask = falco_gen_pupil_customHex( input );
        
    case 'HabEx_B_offaxis'
        input.Nbeam = mp.P1.full.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.aperture_num = mp.P1.aperture_num;
        input.gap_size = mp.P1.gap_size;
        input.fivemmRadii = mp.P1.fivemmRadii;
        
        mp.P1.full.mask = falco_gen_pupil_HabEx_B(input);
        
        input.Nbeam = mp.P1.compact.Nbeam; % number of points across the pupil diameter
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        mp.P1.compact.mask = falco_gen_pupil_HabEx_B( input );
        
        
end
mp.P1.full.Narr = length(mp.P1.full.mask);  %--Total number of pixels across array containing the pupil in the full model. Add 2 pixels to Nbeam when the beam is pixel-centered.
mp.P1.compact.Narr = length(mp.P1.compact.mask);  %--Number of pixels across the array containing the input pupil in the compact model
mp.sumPupil = sum(sum(abs(mp.P1.full.mask).^2));

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


end %--END OF FUNCTION



