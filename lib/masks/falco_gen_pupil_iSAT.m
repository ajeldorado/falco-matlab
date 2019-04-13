% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Inputs: 
% input.Nbeam - Number of samples across the circumscribed pupil diameter
% input.Npad - Number of samples in padded array 
% input.numRaftRings - Number of rings of 7 segment rafts 
% input.raftDia - Diameter of the rafts in samples 
% input.wGap - Gap width within a raft in samples 
% input.raftGap - Gap width bewteen rafts in samples 
% input.raftOffsets - 2 column array of offsets for each raft (number of rafts)x2
% 
% NOTE: segment piston, tip, tilt is not working yet!

function PUPIL = falco_gen_pupil_iSAT( input )
%falco_gen_pupil_iSAT Generates an iSAT pupil

    f2f = input.raftDia/3;% flat-to-flat diameter of a single segment 
    numRaftRings = input.numRaftRings;
   
    % make the central raft of 7 segments
    input_sub = input;
    input_sub.apDia = input.raftDia;
    input_sub.numRings = 1;
    if(isfield(input,'raftOffsets'))
        input_sub.offset = input.raftOffsets(1,:);
    end
    if(isfield(input,'pistons'))
        PUPIL = hexSegMirror_getField( input_sub );
    else
        PUPIL = hexSegMirror_getSupport( input_sub );
    end
    
    % Make the rings of rafts
    count = 1;    
    for ringNum = 1:numRaftRings
    
        raftoffset_norm = ringNum*[-0.5,3*sqrt(3)/2];% vector to center of first raft in ring
        raftSep = ringNum*input.raftGap + sqrt(sum((raftoffset_norm*f2f).^2));% distance from origin to center of first raft in ring 
        theta0 = atan2(raftoffset_norm(1),raftoffset_norm(2));% angle of vector from origin to center of first raft in ring

        for direction_index = 0:5

            step_dir = pi/3*direction_index;
            steprow = raftSep*sin(step_dir+theta0);
            stepcol = raftSep*cos(step_dir+theta0);

            stepnum = 1;
            while(stepnum<=ringNum)

                if(isfield(input,'raftOffsets'))
                    raftPosition = [steprow+input.raftOffsets(1+count,1),stepcol+input.raftOffsets(1+count,2)];
                else
                    raftPosition = [steprow,stepcol];
                end

                ang2next = step_dir+theta0+4*pi/6;
                vec2next = raftSep/ringNum*[sin(ang2next),cos(ang2next)];
                input_sub.offset = raftPosition + (stepnum-1)*vec2next;
                
                % add the raft
                if(isfield(input_sub,'pistons'))
                    PUPIL = PUPIL + hexSegMirror_getField( input_sub );
                else
                    PUPIL = PUPIL + hexSegMirror_getSupport( input_sub );
                end
                count = count + 1;
                stepnum = stepnum + 1;
            end

        end
    end

end