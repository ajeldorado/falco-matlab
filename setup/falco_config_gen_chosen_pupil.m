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

%--Generate/Load Input Pupil
switch upper(mp.whichPupil)
    case{'SIMPLE','SIMPLEPROPER'}
        
        if(strcmpi(mp.whichPupil,'simpleproper'))
            inputs.flagPROPER = true;
        end
        
        inputs.Nbeam = mp.P1.full.Nbeam;     % number of points across the pupil diameter
        inputs.OD = mp.P1.ODnorm;
        inputs.ID = mp.P1.IDnorm;
        inputs.Nstrut = mp.P1.Nstrut;
        inputs.angStrut = mp.P1.angStrut;%Angles of the struts 
        inputs.wStrut = mp.P1.wStrut;% spider width (fraction of the pupil diameter)
        inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));% 
        inputs.stretch = mp.P1.stretch; 
        
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
        
    case{'WFIRST_ONAXIS'}
        inputs.wStrut = mp.pup_wStrut; %--0.0261  is nominal from 2014 on-axis (in pupil diameters)
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
        if(isfield(mp.P1,'wStrut'))
            inputs.wStrut = mp.P1.wStrut;% spider width (fraction of the pupil diameter)
        end
        
        if(isfield(mp.P1,'wGap_m'))
            inputs.wGap_m = mp.P1.wGap_m;% spider width (fraction of the pupil diameter)
        end
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs);
        clear inputs

    case{'LUVOIRA0'}
        inputs.centering = mp.centering;
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_0(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_0(inputs);
        
    case{'LUVOIR_B_OFFAXIS'}
        input.Nbeam = mp.P1.full.Nbeam/0.925; % number of points across the pupil diameter
        input.wGap = mp.P1.wGap*mp.P1.full.Nbeam; % samples
        %input.wGap = 6e-3/7.989*mp.P1.full.Nbeam; % samples
        input.numRings = 4;% Number of rings in hexagonally segmented mirror 
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.ID = 0; % central obscuration radius 
        input.OD = 1; % pupil outer diameter, can be < 1
        input.Nstrut = 0;% Number of struts 
        input.angStrut = [];%Angles of the struts (deg)
        input.wStrut = []; % Width of the struts (fraction of pupil diam.)

        if(isfield(mp.P1,'pistons'))
            input.pistons = mp.P1.pistons;%Tilts on segment in vertical direction (waves/apDia)
        end
        if(isfield(mp.P1,'tiltxs'))
            input.tiltxs = mp.P1.tiltxs;%Tilts on segments (waves/apDia)
            input.tiltys = mp.P1.tiltys;%Tilts on segments (waves/apDia)
        end

        missingSegments = ones(1,hexSegMirror_numSegments( input.numRings ));
        for index = 0:5
            missingSegments(38+index*4) = 0;
        end
        input.missingSegments = missingSegments;

        mp.P1.full.mask = falco_gen_pupil_customHex( input );
        
        input.Nbeam = mp.P1.compact.Nbeam/0.925; % number of points across the pupil diameter
        %input.wGap = 6e-3/7.989*mp.P1.compact.Nbeam; % samples
        input.wGap = mp.P1.wGap*mp.P1.compact.Nbeam; % samples
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        mp.P1.compact.mask = falco_gen_pupil_customHex( input );
        
        if(isfield(mp.P1,'pistons') || isfield(mp.P1,'tiltxs') || isfield(mp.P1,'tiltys'))
            
            % Compact model has one field per sub-bandpass 
            for sbpIndex = 1:mp.Nsbp
                lambda = mp.sbp_centers(sbpIndex);
                phz = angle(mp.P1.compact.mask)*mp.lambda0/lambda;
                mp.P1.compact.E(:,:,sbpIndex) = exp(1i*phz);
            end
            
            % Full model can have more than one field per sub-bandpass
            for sbpIndex = 1:mp.Nsbp
                for wpsbpIndex = 1:mp.Nwpsbp
                    lambda = mp.sbp_centers(sbpIndex)*mp.full.sbp_facs(wpsbpIndex);
                    phz = angle(mp.P1.full.mask)*mp.lambda0/lambda;
                    mp.P1.full.E(:,:,wpsbpIndex,sbpIndex) = exp(1i*phz);
                end
            end
            
            % Masks are amplitude-only, the E cube holds the phase
            % information 
            mp.P1.full.mask = abs(mp.P1.full.mask);
            mp.P1.compact.mask = abs(mp.P1.compact.mask);
            
            if(mp.flagPlot)
                figure;imagesc(angle(mp.P1.full.E(:,:,1,1))/2/pi);axis image;colorbar;title('Phase of telescope aperture (waves)');
                %figure;imagesc(angle(mp.P1.compact.E(:,:,1))/2/pi);axis image;colorbar;title('Phase of telescope aperture (waves)');
            end
        end
    
    case 'DST_LUVOIRB'
        
        mp.P1.full.mask    = falco_gen_pupil_dst_LUVOIR_B(mp.P1.full.Nbeam   ,2^(nextpow2(mp.P1.full.Nbeam   )));
        mp.P1.compact.mask = falco_gen_pupil_dst_LUVOIR_B(mp.P1.compact.Nbeam,2^(nextpow2(mp.P1.compact.Nbeam)));
    
    case 'HABEX_B_OFFAXIS'
        input.Nbeam = mp.P1.full.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.aperture_num = mp.P1.aperture_num;
        input.gap_size = mp.P1.gap_size;
        input.fivemmRadii = mp.P1.fivemmRadii;
        
        mp.P1.full.mask = falco_gen_pupil_HabEx_B(input);
        
        input.Nbeam = mp.P1.compact.Nbeam; % number of points across the pupil diameter
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        mp.P1.compact.mask = falco_gen_pupil_HabEx_B( input );
        
    case 'ISAT'
        input.numRaftRings = 3;
        
        input.Nbeam = mp.P1.full.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.raftDia = input.Nbeam/19*3;% samples
        input.wGap = mp.P1.wGap*input.Nbeam; % samples
        input.raftGap = mp.P1.raftGap*input.Nbeam; % samples
        mp.P1.full.mask    = falco_gen_pupil_iSAT( input );
        
        input.Nbeam = mp.P1.compact.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        input.raftDia = input.Nbeam/19*3;%samples
        input.wGap = mp.P1.wGap*input.Nbeam; % samples
        input.raftGap = mp.P1.raftGap*input.Nbeam; % samples
        mp.P1.compact.mask = falco_gen_pupil_iSAT( input );
    
end
mp.P1.full.Narr = length(mp.P1.full.mask);  %--Total number of pixels across array containing the pupil in the full model. Add 2 pixels to Nbeam when the beam is pixel-centered.
mp.P1.compact.Narr = length(mp.P1.compact.mask);  %--Number of pixels across the array containing the input pupil in the compact model
mp.sumPupil = sum(sum(abs(mp.P1.compact.mask).^2)); %--Throughput is computed with the compact model

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


end %--END OF FUNCTION



