% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_gen_chosen_LS(mp)
%
% Function to generate the Lyot stop representation based on configuration settings.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-05-29 by A.J. Riggs.




function mp = falco_config_gen_chosen_LS(mp)

%% Lyot plane resolution, coordinates, and cropped-down mask for compact model

%--Resolution at Lyot Plane
mp.P4.full.dx = mp.P4.D/mp.P4.full.Nbeam;
mp.P4.compact.dx = mp.P4.D/mp.P4.compact.Nbeam; %mp.NlyotFactor*mp.P4.full.dx;

switch mp.whichPupil
    case{'Simple','SimplePROPER'}
        
        if(strcmpi(mp.whichPupil,'simpleproper'))
            inputs.flagPROPER = true;
        end
        
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
 
        
    case{'WFIRST180718'}

        %--Define Lyot stop generator function inputs for the 'full' optical model
        changes.DCOBS = mp.P4.IDnorm;
        changes.ODpup = mp.P4.ODnorm;
        changes.wStrut = mp.P4.wStrut;
        changes.flagRot180 = true;

        %--Make or read in Lyot stop (LS) for the 'full' model
        mp.P4.full.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P4.full.Nbeam,mp.centering,changes);
        
        %--Make or read in Lyot stop (LS) for the 'compact' model
        mp.P4.compact.mask = falco_gen_pupil_WFIRST_CGI_180718(mp.P4.compact.Nbeam,mp.centering,changes);

        
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
        
        if(isfield(mp,'SPname'))
            switch mp.SPname
                case '20170714'
                    %--Define Lyot stop generator function inputs in a structure
                    inputs.Dbeam = mp.P4.D; % meters;
                    inputs.Nbeam = mp.P4.full.Nbeam; 
                    inputs.ID = mp.P4.IDnorm; % (pupil diameters)
                    inputs.OD = mp.P4.ODnorm; % (pupil diameters)
                    inputs.ang = mp.P4.ang; % (degrees)
                    inputs.centering = mp.centering; % 'interpixel' or 'pixel'


                    %--Make bowtie Lyot stop (LS) for the 'full' model
                    mp.P4.full.mask = falco_gen_bowtie_LS(inputs);

                    %--Make bowtie Lyot stop (LS) for the 'compact' model
                    inputs.Nbeam = mp.P4.compact.Nbeam; 
                    mp.P4.compact.mask = falco_gen_bowtie_LS(inputs);               
            end
        end
        
	case{'LUVOIRA5','LUVOIRA0'}
        
        
        %--Define Lyot stop generator function inputs for the 'full' optical model
        inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam  
        inputs.Dbeam = mp.P1.D;
        inputs.ID = mp.P4.IDnorm;
        inputs.OD = mp.P4.ODnorm;
        inputs.strut_width = mp.LS_strut_width;
        inputs.centering = mp.centering;
        %--Make or read in Lyot stop (LS) for the 'full' model
        mp.P4.full.mask = falco_gen_pupil_LUVOIR_A_5_Lyot_struts(inputs,'ROT180');
        
        %--Make or read in Lyot stop (LS) for the 'compact' model
        inputs.Nbeam = mp.P4.compact.Nbeam;     % number of points across incoming beam           
        mp.P4.compact.mask = falco_gen_pupil_LUVOIR_A_5_Lyot_struts(inputs,'ROT180');

%         %--Define Lyot stop generator function inputs for the 'full' optical model
%         inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam  
%         inputs.Dbeam = mp.P4.D; %--diameter of the beam at the mask (meters)
%         %inputs.Narray = mp.Nlyot;   % number of points across output array
%         inputs.ID = mp.P4.IDnorm;
%         inputs.OD = mp.P4.ODnorm;
%         %inputs.strut_width = mp.LS_strut_width;
%         inputs.strut_width = 0;
%         inputs.centering = mp.centering;
%         %--Make or read in Lyot stop (LS) for the 'full' model
%         mp.P4.full.mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,'ROT180');
% 
%         %--Make or read in Lyot stop (LS) for the 'compact' model
%         inputs.Nbeam = mp.P4.compact.Nbeam;     % number of points across incoming beam           
%         %inputs.Narray = mp.NlyotCompact;   % number of points across output array
%         mp.P4.compact.mask = falco_gen_pupil_WFIRSTcycle6_LS(inputs,'ROT180');
        
    case {'LUVOIR_B_offaxis','HabEx_B_offaxis'}
        %--Full model
        inputs.Nbeam = mp.P4.full.Nbeam;     % number of points across incoming beam 
        inputs.Npad = 2^(nextpow2(mp.P4.full.Nbeam));% 
        inputs.OD = mp.P4.ODnorm;
        inputs.ID = mp.P4.IDnorm;
        inputs.num_strut = 0;
        inputs.strut_angs = [];%Angles of the struts 
        inputs.strut_width = 0;% spider width (fraction of the pupil diameter)

        mp.P4.full.mask = falco_gen_pupil_Simple( inputs );
        
        pad_pct = mp.P4.padFacPct;
        if(pad_pct>0) %--Also apply an eroded/padded version of the segment gaps

            pupil0 = mp.P1.full.mask;
            Nbeam = inputs.Nbeam;
            Npad = inputs.Npad;

            xsD = ( -Npad/2:(Npad/2-1) )/Nbeam; %--coordinates, normalized to the pupil diameter
            [XS,YS] = meshgrid(xsD);
            RS = sqrt(XS.^2 + YS.^2);
        
            pupil1 = 1-pupil0;
            %figure(7); imagesc(xsD,xsD,pupil1); axis xy equal tight; colormap gray;

            spot = zeros(Npad);
            spot(RS <= pad_pct/100) = 1;
            % figure(8); imagesc(xsD,xsD,spot); axis xy equal tight; colormap gray;

            pupil4 = ifftshift(ifft2(  fft2(fftshift(pupil1)).*fft2(fftshift(spot))   ));
            pupil4 = abs(pupil4);
            pupil4 = pupil4/max(pupil4(:));
            % figure(9); imagesc(xsD,xsD,pupil4); axis xy equal tight; colormap gray;

            pupil5 = 1-pupil4;
            %figure(10); imagesc(xsD,xsD,pupil5); axis xy equal tight; colormap gray;

            thresh = 0.99;
            pupil5(pupil5<thresh) = 0;
            pupil5(pupil5>=thresh) = 1;
            %figure(12); imagesc(xsD,xsD,pupil5); axis xy equal tight; colormap gray;

            mp.P4.full.mask = mp.P4.full.mask.*pupil5;            
        end
        %figure(13); imagesc(xsD,xsD,mp.P4.full.mask); axis xy equal tight; colormap gray;

        
        %--Compact model
        inputs.Nbeam = mp.P4.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
        inputs.Npad = 2^(nextpow2(mp.P4.compact.Nbeam));% 
        
        mp.P4.compact.mask = falco_gen_pupil_Simple( inputs );
        
        if(pad_pct>0) %--Also apply an eroded/padded version of the segment gaps
            pupil0 = mp.P1.compact.mask;
            Nbeam = inputs.Nbeam;
            Npad = inputs.Npad;

            xsD = ( -Npad/2:(Npad/2-1) )/Nbeam; %--coordinates, normalized to the pupil diameter
            [XS,YS] = meshgrid(xsD);
            RS = sqrt(XS.^2 + YS.^2);

            pupil1 = 1-pupil0;
            %figure(7); imagesc(xsD,xsD,pupil1); axis xy equal tight; colormap gray;

            spot = zeros(Npad);
            spot(RS <= pad_pct/100) = 1;
            % figure(8); imagesc(xsD,xsD,spot); axis xy equal tight; colormap gray;

            pupil4 = ifftshift(ifft2(  fft2(fftshift(pupil1)).*fft2(fftshift(spot))   ));
            pupil4 = abs(pupil4);
            pupil4 = pupil4/max(pupil4(:));
            % figure(9); imagesc(xsD,xsD,pupil4); axis xy equal tight; colormap gray;

            pupil5 = 1-pupil4;
            %figure(10); imagesc(xsD,xsD,pupil5); axis xy equal tight; colormap gray;  drawnow;

            thresh = 0.99;
            pupil5(pupil5<thresh) = 0;
            pupil5(pupil5>=thresh) = 1;
            %figure(12); imagesc(xsD,xsD,pupil5); axis xy equal tight; colormap gray;

            mp.P4.compact.mask = mp.P4.compact.mask.*pupil5;
        end
        %figure(14); imagesc(xsD,xsD,mp.P4.compact.mask); axis xy equal tight; colormap gray;  drawnow;      
    
       
        
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
mp.P4.compact.ys = mp.P4.compact.xs.';



end %--END OF FUNCTION



