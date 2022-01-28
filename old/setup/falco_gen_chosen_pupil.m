% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the chosen entrance pupil aperture.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_gen_chosen_pupil(mp)

%% Input pupil plane resolution, masks, and coordinates


%--Generate/Load Input Pupil
switch upper(mp.whichPupil)
    case{'SIMPLE'}
        
        inputs.Nbeam = mp.P1.full.Nbeam; % number of points across the pupil diameter
        inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        inputs.OD = mp.P1.ODnorm;
        if(isfield(mp.P1, 'IDnorm')); inputs.ID = mp.P1.IDnorm; end %else; inputs.ID = 0; end
        if(isfield(mp.P1, 'angStrut')); inputs.angStrut = mp.P1.angStrut; end
        if(isfield(mp.P1, 'wStrut')); inputs.wStrut = mp.P1.wStrut; end % spider width (fraction of the pupil diameter)
        if(isfield(mp.P1, 'stretch')); inputs.stretch = mp.P1.stretch; end % else; inputs.stretch = 1; end
        if(isfield(inputs,'centering')); inputs.centering = mp.P1.centering; else; inputs.centering = mp.centering; end
        if(isfield(inputs, 'clocking')); inputs.clocking = mp.P1.clocking; end % clocking [degrees]
        if(isfield(inputs, 'xShear')); inputs.xShear = mp.P1.xShear; end % x-shear [pupil diameters]
        if(isfield(inputs, 'yShear')); inputs.yShear = mp.P1.yShear;  end % y-shear [pupil diameters]
        if(isfield(mp.P1, 'flagHG')); inputs.flagHG = mp.P1.flagHG; end % whether to use hyper-Gaussians to generate the pupil instead. (Old behavior)

        mp.P1.full.mask = falco_gen_pupil_Simple(inputs);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; % number of points across usable pupil   
        inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam)); % number of points across usable pupil 
        mp.P1.compact.mask = falco_gen_pupil_Simple(inputs);

    case{'WFIRST20200513', 'ROMAN20200513', 'ROMAN'}
        if(mp.full.flagGenPupil);  mp.P1.full.mask = falco_gen_pupil_Roman_CGI_20200513(mp.P1.full.Nbeam, mp.centering);  end
        if(mp.compact.flagGenPupil);  mp.P1.compact.mask = falco_gen_pupil_Roman_CGI_20200513(mp.P1.compact.Nbeam, mp.centering);  end
    
    case{'WFIRST20191009'}
        if(mp.full.flagGenPupil);  mp.P1.full.mask = falco_gen_pupil_WFIRST_CGI_20191009(mp.P1.full.Nbeam, mp.centering);  end
        if(mp.compact.flagGenPupil);  mp.P1.compact.mask = falco_gen_pupil_WFIRST_CGI_20191009(mp.P1.compact.Nbeam, mp.centering);  end
    
    case{'WFIRST180718'}
        if(mp.full.flagGenPupil);  mp.P1.full.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P1.full.Nbeam, mp.centering);  end
        if(mp.compact.flagGenPupil);  mp.P1.compact.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P1.compact.Nbeam, mp.centering);  end

    case{'WFIRST2016', 'WFIRST2016ONAXIS'}
        if(isfield(mp.P1, 'wStrut')); changes.wStrut = mp.P1.wStrut; end % strut width (fraction of the pupil diameter)
        if(mp.full.flagGenPupil);  mp.P1.full.mask = falco_gen_pupil_WFIRST_2016_onaxis(mp.P1.full.Nbeam, mp.centering, changes);  end
        if(mp.compact.flagGenPupil);  mp.P1.compact.mask = falco_gen_pupil_WFIRST_2016_onaxis(mp.P1.compact.Nbeam, mp.centering, changes);  end
            
    case{'KECK'}
        inputs.centering = mp.centering; 
        
        %--Clocking angle, if defined [degrees]
        if(isfield(mp.P1,'clock_deg'))
            angRot = mp.P1.clock_deg;
        else
            angRot = 0;
        end
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        mp.P1.full.mask = falco_gen_pupil_Keck(inputs,'ROTATION',angRot);
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        mp.P1.compact.mask = falco_gen_pupil_Keck(inputs,'ROTATION',angRot); 
        
    case{'LUVOIRAFINAL'}
        if(mp.compact.flagGenPupil || mp.full.flagGenPupil)
            inputs.centering = mp.centering;
            if(isfield(mp.P1,'wGap_m')) %--Option to overwrite the default
                inputs.wGap_m = mp.P1.wGap_m;
            end
        end
        
        if(mp.full.flagGenPupil)
            %--Generate high-res input pupil for the 'full' model
            inputs.Nbeam = mp.P1.full.Nbeam; 
            mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_final(inputs);
        end
        
        if(mp.compact.flagGenPupil)
            %--Generate low-res input pupil for the 'compact' model
            inputs.Nbeam = mp.P1.compact.Nbeam; 
            mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_final(inputs);
        end
        
    case{'LUVOIRA5'}
        inputs.centering = mp.centering;
        if(isfield(mp.P1,'wStrut'))  %--Option to overwrite the default
            inputs.wStrut = mp.P1.wStrut; % spider width (fraction of the pupil diameter)
        end
        if(isfield(mp.P1,'wGap_m'))  %--Option to overwrite the default
            inputs.wGap_m = mp.P1.wGap_m; % spider width (fraction of the pupil diameter)
        end
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        if(mp.full.flagGenPupil); mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs); end
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        if(mp.compact.flagGenPupil); mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_5_mag_trans(inputs); end
        clear inputs

    case{'LUVOIRA0'}
        inputs.centering = mp.centering;
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam; 
        if(mp.full.flagGenPupil); mp.P1.full.mask = falco_gen_pupil_LUVOIR_A_0(inputs); end
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam; 
        if(mp.compact.flagGenPupil); mp.P1.compact.mask = falco_gen_pupil_LUVOIR_A_0(inputs); end
        
    case{'LUVOIR_B', 'LUVOIRB'} % simple version of LUVOIR B using PROPER
        inputs.centering = mp.centering;
        
        %--Generate high-res input pupil for the 'full' model
        inputs.Nbeam = mp.P1.full.Nbeam;
        Narray = 2^(nextpow2(inputs.Nbeam));
        if(mp.full.flagGenPupil)
            mp.P1.full.mask = pad_crop(falco_gen_pupil_LUVOIR_B(inputs), Narray);
        end
        
        %--Generate low-res input pupil for the 'compact' model
        inputs.Nbeam = mp.P1.compact.Nbeam;
        Narray = 2^(nextpow2(inputs.Nbeam));
        if(mp.compact.flagGenPupil)
            mp.P1.compact.mask = pad_crop(falco_gen_pupil_LUVOIR_B(inputs), Narray);
        end
        
    case{'LUVOIR_B_OFFAXIS'} % more complicated version with more features
        input.Nbeam = mp.P1.full.Nbeam/0.925; % number of points across the pupil diameter
        input.wGap = mp.P1.wGap*mp.P1.full.Nbeam; % samples
        input.numRings = 4; % Number of rings in hexagonally segmented mirror 
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.ID = 0; % central obscuration radius 
        input.OD = 1; % pupil outer diameter, can be < 1
        input.Nstrut = 0; % Number of struts 
        input.angStrut = []; % Angles of the struts (deg)
        input.wStrut = []; % Width of the struts (fraction of pupil diam.)

        if(isfield(mp.P1,'pistons'))
            input.pistons = mp.P1.pistons; % Tilts on segment in vertical direction (waves/apDia)
        end
        if(isfield(mp.P1,'tiltxs'))
            input.tiltxs = mp.P1.tiltxs; % Tilts on segments (waves/apDia)
        end
        if(isfield(mp.P1,'tiltys'))
            input.tiltys = mp.P1.tiltys; % Tilts on segments (waves/apDia)
        end
        if(isfield(mp.P1,'loworder_struct'))
            input.loworder_struct = mp.P1.loworder_struct; % Segment-level low order aberration structure
        end

        missingSegments = ones(1,hexSegMirror_numSegments(input.numRings));
        for index = 0:5
            missingSegments(38+index*4) = 0;
        end
        input.missingSegments = missingSegments;

        if(mp.full.flagGenPupil); mp.P1.full.mask = falco_gen_pupil_customHex(input); end
        
        input.Nbeam = mp.P1.compact.Nbeam/0.925; % number of points across the pupil diameter
        input.wGap = mp.P1.wGap*mp.P1.compact.Nbeam; % samples
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        if(mp.compact.flagGenPupil); mp.P1.compact.mask = falco_gen_pupil_customHex(input); end
        
        if(isfield(mp.P1,'pistons') || isfield(mp.P1,'tiltxs') || isfield(mp.P1,'tiltys') || isfield(mp.P1,'loworder_struct'))
            
            % Compact model has one field per sub-bandpass 
            for sbpIndex = 1:mp.Nsbp
                lambda = mp.sbp_centers(sbpIndex);
                phz = angle(mp.P1.compact.mask)*mp.lambda0/lambda;
                mp.P1.compact.E(:,:,sbpIndex) = exp(1i*phz);
            end
            
            % Full model can have more than one field per sub-bandpass
            for sbpIndex = 1:mp.Nsbp
                for wpsbpIndex = 1:mp.Nwpsbp
                    lambda = mp.full.lambdasMat(sbpIndex,wpsbpIndex);
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
            end
        end
    
    case 'DST_LUVOIRB'
        
        if(mp.full.flagGenPupil); mp.P1.full.mask = falco_gen_pupil_dst_LUVOIR_B(mp.P1.full.Nbeam, 2^(nextpow2(mp.P1.full.Nbeam))); end
        if(mp.compact.flagGenPupil); mp.P1.compact.mask = falco_gen_pupil_dst_LUVOIR_B(mp.P1.compact.Nbeam, 2^(nextpow2(mp.P1.compact.Nbeam))); end
    
    case 'HABEX_B_OFFAXIS'
        input.Nbeam = mp.P1.full.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.aperture_num = mp.P1.aperture_num;
        input.gap_size = mp.P1.gap_size;
        input.fivemmRadii = mp.P1.fivemmRadii;
        
        if(mp.full.flagGenPupil); mp.P1.full.mask = falco_gen_pupil_HabEx_B(input); end
        
        input.Nbeam = mp.P1.compact.Nbeam; % number of points across the pupil diameter
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        if(mp.compact.flagGenPupil); mp.P1.compact.mask = falco_gen_pupil_HabEx_B(input); end
        
    case 'ISAT'
        input.numRaftRings = 3;
        
        input.Nbeam = mp.P1.full.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
        input.raftDia = input.Nbeam/19*3; % samples
        input.wGap = mp.P1.wGap*input.Nbeam; % samples
        input.raftGap = mp.P1.raftGap*input.Nbeam; % samples
        if(isfield(mp.P1,'raftOffsets'))
            input.raftOffsets = mp.P1.raftOffsets*input.Nbeam;% samples 
        end
        if(mp.full.flagGenPupil); mp.P1.full.mask = falco_gen_pupil_iSAT(input); end
        
        input.Nbeam = mp.P1.compact.Nbeam;
        input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
        input.raftDia = input.Nbeam/19*3; % samples
        input.wGap = mp.P1.wGap*input.Nbeam; % samples
        input.raftGap = mp.P1.raftGap*input.Nbeam; % samples
        if(isfield(mp.P1,'raftOffsets'))
            input.raftOffsets = mp.P1.raftOffsets*input.Nbeam;% samples 
        end
        if(mp.compact.flagGenPupil); mp.P1.compact.mask = falco_gen_pupil_iSAT(input); end
end

end %--END OF FUNCTION