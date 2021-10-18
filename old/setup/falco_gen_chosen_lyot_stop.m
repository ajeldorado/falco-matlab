% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate the chosen lyot stop.
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_gen_chosen_lyot_stop(mp)

%% Lyot plane resolution, coordinates, and cropped-down mask for compact model

switch upper(mp.whichPupil)
    case{'SIMPLE','DST_LUVOIRB','ISAT'}
        
        if(mp.compact.flagGenLS || mp.full.flagGenLS)            
            inputs.OD = mp.P4.ODnorm;
            if(isfield(mp.P4, 'IDnorm')); inputs.ID = mp.P4.IDnorm; end %else; inputs.ID = 0; end
            if(isfield(mp.P4, 'angStrut')); inputs.angStrut = mp.P4.angStrut; end
            if(isfield(mp.P4, 'wStrut')); inputs.wStrut = mp.P4.wStrut; end % spider width (fraction of the pupil diameter)
            if(isfield(mp.P4, 'stretch')); inputs.stretch = mp.P4.stretch; end % else; inputs.stretch = 1; end
            if(isfield(inputs,'centering')); inputs.centering = mp.P4.centering; else; inputs.centering = mp.centering; end
            if(isfield(inputs, 'clocking')); inputs.clocking = mp.P4.clocking; end % clocking [degrees]
            if(isfield(inputs, 'xShear')); inputs.xShear = mp.P4.xShear; end % x-shear [pupil diameters]
            if(isfield(inputs, 'yShear')); inputs.yShear = mp.P4.yShear;  end % y-shear [pupil diameters]
            if(isfield(mp.P4, 'flagHG')); inputs.flagHG = mp.P4.flagHG; end % whether to use hyper-Gaussians to generate the pupil instead. (Old behavior)
        end

        if(mp.full.flagGenLS)
            inputs.Nbeam = mp.P4.full.Nbeam; % number of points across incoming beam 
            inputs.Npad = 2^(nextpow2(mp.P4.full.Nbeam));
            mp.P4.full.mask = falco_gen_pupil_Simple(inputs); 
        end

        if(mp.compact.flagGenLS)
            inputs.Nbeam = mp.P4.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
            inputs.Npad = 2^(nextpow2(mp.P4.compact.Nbeam));
            mp.P4.compact.mask = falco_gen_pupil_Simple(inputs); 
        end
 
    case{'ROMAN', 'ROMAN20200513', 'WFIRST20200513','WFIRST20191009', 'WFIRST180718', 'WFIRST2016', 'WFIRST2016ONAXIS'}

        %--Define Lyot stop generator function inputs for the 'full' optical model
        if(mp.compact.flagGenLS || mp.full.flagGenLS)
            changes.flagLyot = true;
            changes.ID = mp.P4.IDnorm;
            changes.OD = mp.P4.ODnorm;
            changes.wStrut = mp.P4.wStrut;
            changes.flagRot180 = true;
        end
        
        switch upper(mp.whichPupil)
            case{'WFIRST20200513', 'ROMAN20200513', 'ROMAN'}
                if mp.P4.flagSymm == false
                    if(mp.full.flagGenLS); mp.P4.full.mask = falco_gen_pupil_Roman_CGI_20200513(mp.P4.full.Nbeam, mp.centering, changes); end
                    if(mp.compact.flagGenLS); mp.P4.compact.mask = falco_gen_pupil_Roman_CGI_20200513(mp.P4.compact.Nbeam, mp.centering, changes); end
                else
                    rocFillet = 0.02;
                    upsampleFactor = 100;
                    if(mp.full.flagGenLS); mp.P4.full.mask = propcustom_relay(falco_gen_Roman_CGI_lyot_stop_symm_fillet(mp.P4.full.Nbeam, changes.ID, changes.OD, changes.wStrut, rocFillet, upsampleFactor, mp.centering), 1, mp.centering); end
                    if(mp.compact.flagGenLS); mp.P4.compact.mask = propcustom_relay(falco_gen_Roman_CGI_lyot_stop_symm_fillet(mp.P4.compact.Nbeam, changes.ID, changes.OD, changes.wStrut, rocFillet, upsampleFactor, mp.centering), 1, mp.centering); end
                end
            
            case{'WFIRST20191009'}
                if(mp.full.flagGenLS); mp.P4.full.mask = falco_gen_pupil_WFIRST_CGI_20191009(mp.P4.full.Nbeam, mp.centering, changes); end
                if(mp.compact.flagGenLS); mp.P4.compact.mask = falco_gen_pupil_WFIRST_CGI_20191009(mp.P4.compact.Nbeam, mp.centering, changes); end
                
            case{'WFIRST180718'}
                if(mp.full.flagGenLS); mp.P4.full.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P4.full.Nbeam,mp.centering,changes); end
                if(mp.compact.flagGenLS); mp.P4.compact.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P4.compact.Nbeam,mp.centering,changes); end
            
            case{'WFIRST2016', 'WFIRST2016ONAXIS'}
                if(mp.full.flagGenLS); mp.P4.full.mask = falco_gen_pupil_WFIRST_2016_onaxis(mp.P4.full.Nbeam,mp.centering,changes); end
                if(mp.compact.flagGenLS); mp.P4.compact.mask = falco_gen_pupil_WFIRST_2016_onaxis(mp.P4.compact.Nbeam,mp.centering,changes); end
        end
        
        if(isfield(mp,'LSshape'))
            switch lower(mp.LSshape)
                case 'bowtie'
                    %--Define Lyot stop generator function inputs in a structure
                    inputs.Dbeam = mp.P4.D; % meters;
                    inputs.ID = mp.P4.IDnorm; % (pupil diameters)
                    inputs.OD = mp.P4.ODnorm; % (pupil diameters)
                    inputs.ang = mp.P4.ang; % (degrees)
                    inputs.centering = mp.centering; % 'interpixel' or 'pixel'
                    if(isfield(mp.P4, 'clocking')); inputs.clocking = mp.P4.clocking; end
                    if(isfield(mp.P4, 'Rfillet')); inputs.Rfillet = mp.P4.Rfillet; end

                    if(mp.full.flagGenLS)
                        inputs.Nbeam = mp.P4.full.Nbeam; 
                        mp.P4.full.mask = falco_gen_bowtie_LS(inputs);
                    end
                    
                    %--Make bowtie Lyot stop (LS) for the 'compact' model
                    if(mp.compact.flagGenLS)
                        inputs.Nbeam = mp.P4.compact.Nbeam; 
                        mp.P4.compact.mask = falco_gen_bowtie_LS(inputs);  
                    end
            end
        end
        
    case{'LUVOIRAFINAL'}
        if(mp.compact.flagGenLS || mp.full.flagGenLS)
            %--Define Lyot stop generator function inputs
            inputs.flagLyot = true;
            inputs.ID = mp.P4.IDnorm;
            inputs.OD = mp.P4.ODnorm;
            inputs.wStrut = mp.P4.wStrut;
            inputs.centering = mp.centering;
            inputs.flagRot180 = true;
        end
        
        if(mp.full.flagGenLS)
            %--Make or read in Lyot stop (LS) for the 'full' model
            inputs.Nbeam = mp.P4.full.Nbeam; % number of points across incoming beam  
            mp.P4.full.mask = falco_gen_pupil_LUVOIR_A_final(inputs);
        end
        
        if(mp.compact.flagGenLS)
            %--Make or read in Lyot stop (LS) for the 'compact' model
            inputs.Nbeam = mp.P4.compact.Nbeam;     % number of points across incoming beam           
            mp.P4.compact.mask = falco_gen_pupil_LUVOIR_A_final(inputs);
        end
        
	case{'LUVOIRA5'}
        
        %--Define Lyot stop generator function inputs for the 'full' optical model
        inputs.Nbeam = mp.P4.full.Nbeam; % number of points across incoming beam  
        inputs.Dbeam = mp.P1.D;
        inputs.ID = mp.P4.IDnorm;
        inputs.OD = mp.P4.ODnorm;
        inputs.wStrut = mp.P4.wStrut;
        inputs.centering = mp.centering;
        %--Make or read in Lyot stop (LS) for the 'full' model
        if(mp.full.flagGenLS); mp.P4.full.mask = falco_gen_pupil_LUVOIR_A_5_Lyot_struts(inputs,'ROT180'); end
        
        %--Make or read in Lyot stop (LS) for the 'compact' model
        inputs.Nbeam = mp.P4.compact.Nbeam;     % number of points across incoming beam           
        if(mp.compact.flagGenLS); mp.P4.compact.mask = falco_gen_pupil_LUVOIR_A_5_Lyot_struts(inputs,'ROT180'); end
        
    case {'LUVOIR_B_OFFAXIS','HABEX_B_OFFAXIS', 'LUVOIR_B', 'LUVOIRB'}
        %--Full model
        inputs.Nbeam = mp.P4.full.Nbeam; % number of points across incoming beam 
        inputs.Npad = 2^(nextpow2(mp.P4.full.Nbeam));
        inputs.OD = mp.P4.ODnorm;
        inputs.ID = mp.P4.IDnorm;
        inputs.Nstrut = 0;
        inputs.angStrut = []; % Angles of the struts 
        inputs.wStrut = 0; % spider width (fraction of the pupil diameter)

        if(mp.full.flagGenLS)
            mp.P4.full.mask = falco_gen_pupil_Simple(inputs);
        
            pad_pct = mp.P4.padFacPct;
            if(pad_pct>0) %--Also apply an eroded/padded version of the segment gaps

                pupil0 = mp.P1.full.mask;
                Nbeam = inputs.Nbeam;
                Npad = inputs.Npad;

                xsD = (-Npad/2:(Npad/2-1))/Nbeam; %--coordinates, normalized to the pupil diameter
                [XS,YS] = meshgrid(xsD);
                RS = sqrt(XS.^2 + YS.^2);

                pupil1 = 1-pupil0;

                spot = zeros(Npad);
                spot(RS <= pad_pct/100) = 1;

                pupil4 = ifftshift(ifft2(fft2(fftshift(pupil1)).*fft2(fftshift(spot))));
                pupil4 = abs(pupil4);
                pupil4 = pupil4/max(pupil4(:));

                pupil5 = 1-pupil4;

                thresh = 0.99;
                pupil5(pupil5<thresh) = 0;
                pupil5(pupil5>=thresh) = 1;

                mp.P4.full.mask = mp.P4.full.mask.*pupil5;            
            end
        end
        
        %--Compact model
        inputs.Nbeam = mp.P4.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
        inputs.Npad = 2^(nextpow2(mp.P4.compact.Nbeam));
        
        if(mp.compact.flagGenLS) 
            mp.P4.compact.mask = falco_gen_pupil_Simple(inputs);
        
            if(pad_pct>0) %--Also apply an eroded/padded version of the segment gaps
                pupil0 = mp.P1.compact.mask;
                Nbeam = inputs.Nbeam;
                Npad = inputs.Npad;

                xsD = (-Npad/2:(Npad/2-1))/Nbeam; %--coordinates, normalized to the pupil diameter
                [XS,YS] = meshgrid(xsD);
                RS = sqrt(XS.^2 + YS.^2);

                pupil1 = 1-pupil0;

                spot = zeros(Npad);
                spot(RS <= pad_pct/100) = 1;

                pupil4 = ifftshift(ifft2(fft2(fftshift(pupil1)).*fft2(fftshift(spot))));
                pupil4 = abs(pupil4);
                pupil4 = pupil4/max(pupil4(:));

                pupil5 = 1-pupil4;

                thresh = 0.99;
                pupil5(pupil5<thresh) = 0;
                pupil5(pupil5>=thresh) = 1;

                mp.P4.compact.mask = mp.P4.compact.mask.*pupil5;
            end
        end
end

end %--END OF FUNCTION
