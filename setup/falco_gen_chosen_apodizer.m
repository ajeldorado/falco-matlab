% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to generate the apodizer representation based on configuration settings.
%
% 
% REVISION HISTORY:
% ----------------
% Modified on 2019-06-19 by A.J. Riggs to be greatly simplified. Is now
% independent of the coronagraph type.
% Modified on 2019-06-13 by Garreth Ruane to update how the vortex apodizer
% is called.
% Created on 2018-05-29 by A.J. Riggs.

function mp = falco_gen_chosen_apodizer(mp)


if(mp.flagApod)
    
    switch lower(mp.apodType) %--mp.apodType is used only when generating certain types of analytical apodizers
        case{'simple'} %--A simple, circular aperture stop

            %--Inputs common to both the compact and full models
            inputs.OD = mp.P3.ODnorm;
            inputs.ID = mp.P3.IDnorm;
            inputs.Nstrut = 0;
            inputs.angStrut = []; %Angles of the struts 
            inputs.wStrut = 0; % spider width (fraction of the pupil diameter)

            %--Full model only
            inputs.Nbeam = mp.P1.full.Nbeam; % number of points across incoming beam 
            inputs.Npad = 2^(nextpow2(mp.P1.full.Nbeam)); 
            if(mp.full.flagGenApod)
                mp.P3.full.mask= falco_gen_pupil_Simple( inputs );
            else
                disp('*** Simple aperture stop to be loaded instead of generated for full model. ***')
            end

            % Compact model only 
            inputs.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
            inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam)); 
            if(mp.compact.flagGenApod)
                mp.P3.compact.mask= falco_gen_pupil_Simple( inputs );
            else
                disp('*** Simple aperture stop to be loaded instead of generated for compact model. ***')
            end

            
        case{'ring'} %--Concentric ring apodizer
            %--Full model
            if(mp.full.flagGenApod)
                mp.P3.full.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.full.dx,mp.P2.D,mp.centering); %--Generate binary-amplitude ring apodizer for the full model
            else
                disp('*** Concentric ring apodizer loaded instead of generated for full model. ***')
            end
            %--Compact model
            if(mp.compact.flagGenApod)
                mp.P3.compact.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.compact.dx,mp.P2.D,mp.centering); %--Generate binary-amplitude ring apodizer for the compact model
            else
                disp('*** Concentric ring apodizer loaded instead of generated for compact model. ***')
            end
            if(mp.flagPlot);  figure(504); imagesc(padOrCropEven(mp.P3.full.mask,length(mp.P1.full.mask)) + mp.P1.full.mask); axis xy equal tight; colorbar; drawnow;  end

            
        case{'grayscale','traditional','classical'} %--A grayscale aperture generated with a Gerchberg-Saxton type algorithm

            if(~strcmpi(mp.centering,'pixel'));  error('Use pixel centering for APLC');  end
            if(mp.P1.full.Nbeam ~= mp.P1.compact.Nbeam);  error('Tradional apodizer generation currently requires Nbeam for the full and compact (in order to use the same apodizer).');  end

            %--Generate the grayscale apodizer
            if(mp.full.flagGenApod && mp.compact.flagGenApod) 
                mp.P3.full.mask = falco_gen_tradApodizer(mp.P1.full.mask,mp.P1.full.Nbeam,mp.F3.Rin,(1+mp.fracBW/2)*mp.F3.Rout,mp.useGPU);
                mp.P3.full.Narr = length(mp.P3.full.mask);
                mp.P3.compact.mask = mp.P3.full.mask;
                mp.P3.compact.Narr = length(mp.P3.compact.mask);
            else
                disp('*** Grayscale apodizer loaded instead of generated for full and compact models. ***')
            end
            if(mp.flagPlot);  figure(504); imagesc(padOrCropEven(mp.P3.full.mask,length(mp.P1.full.mask)).*mp.P1.full.mask); axis xy equal tight; colorbar; drawnow;  end
            
            
        otherwise
            disp('No apodizer type specified for generation.')
    end

    mp.P3.full.dummy = 1; 
    if(isfield(mp.P3.full,'mask'))%    ==false || isfield(mp.P3.compact,'mask')==false)
        mp.P3.full.Narr = length(mp.P3.full.mask);
    else
        fprintf('*** If not generated or loaded in a PROPER model, the apodizer must be loaded \n    in the main script or config file into the variable mp.P3.full.mask ***\n')
    end
    
    mp.P3.compact.dummy = 1;
    if(isfield(mp.P3.compact,'mask'))%    ==false || isfield(mp.P3.compact,'mask')==false)
        mp.P3.compact.Narr = length(mp.P3.compact.mask);
    else
        fprintf('*** If not generated, the apodizer must be loaded in the main script or config \n    file into the variable mp.P3.compact.mask ***\n')
    end
    
end

% %--Set the pixel width [meters]
mp.P3.full.dx = mp.P2.full.dx; 
mp.P3.compact.dx = mp.P2.compact.dx;        
        
end %--END OF FUNCTION