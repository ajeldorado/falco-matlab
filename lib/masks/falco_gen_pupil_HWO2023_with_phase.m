% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the HWO2023 telescope pupil and phase aberrations
% with FALCO's hexSegMirror functions.
%
function mp = falco_gen_pupil_HWO2023_with_phase(mp)

    input.Nbeam = mp.P1.full.Nbeam; % number of points across the pupil diameter
    input.wGap = mp.P1.wGap*mp.P1.full.Nbeam; % samples
    input.numRings = 2; % Number of rings in hexagonally segmented mirror 
    input.Npad = 2^(nextpow2(mp.P1.full.Nbeam));
    input.ID = 0; % central obscuration radius 
    input.OD = 1.1; % pupil outer diameter, can be < 1
    input.Nstrut = 0; % Number of struts 
    input.angStrut = []; % Angles of the struts (deg)
    input.wStrut = []; % Width of the struts (fraction of pupil diam.)
    input.hg_expon = 1000; 
    
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
% 
%     missingSegments = ones(1,hexSegMirror_numSegments(input.numRings));
%     for index = 0:5
%         missingSegments(38+index*4) = 0;
%     end
%     input.missingSegments = missingSegments;

    mp.P1.full.mask = falco_gen_pupil_customHex(input);

    input.Nbeam = mp.P1.compact.Nbeam; % number of points across the pupil diameter
    input.wGap = mp.P1.wGap*mp.P1.compact.Nbeam; % samples
    input.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));
    mp.P1.compact.mask = falco_gen_pupil_customHex(input);

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

        % Masks are amplitude-only, the E cube holds the phase information 
        mp.P1.full.mask = abs(mp.P1.full.mask);
        mp.P1.compact.mask = abs(mp.P1.compact.mask);

        if(mp.flagPlot)
            figure;imagesc(angle(mp.P1.full.E(:,:,1,1))/2/pi);axis image;colorbar;title('Phase of telescope aperture (waves)');
        end
    end

end %--END OF FUNCTION
