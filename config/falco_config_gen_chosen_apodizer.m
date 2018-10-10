% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_gen_chosen_apodizer(mp)
%
% Function to generate the apodizer representation based on configuration settings.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-05-29 by A.J. Riggs.




function mp = falco_config_gen_chosen_apodizer(mp)

switch mp.coro
    
	case{'LC','HLC','EHLC'}
         disp('Using Lyot coronagraph without apodizer or aperture stop.')
        
    case{'SPLC','SPHLC'}
        

        
        switch mp.SPname
            
            case '20170714'
                load('WFIRST_CGI_apod_SPC_20170714.mat','SP'); %--Contains the 1001x1001 matrix named "SP"
                temp = zeros(ceil_even(length(SP)));
                temp(2:end,2:end) = SP;
                mp.P3.full.mask = temp;
                
                %--NOTE: ADD A PART HERE TO DOWNSAMPLE THE APODIZER
                mp.P3.compact.mask = temp;
                
            otherwise
        
                mp.P3.full.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.full.dx,mp.P2.D,mp.centering);
                
                %--Generate lower-resolution SP for the compact model
                mp.P3.compact.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.compact.dx,mp.P2.D,mp.centering);

                %--Shaped Pupil Plane Coordinates (meters)
                if( strcmpi(mp.centering,'interpixel') )
                    % mp.P3.full.xs = ( -(mp.P3.full.Narr-1)/2:(mp.P3.full.Narr-1)/2 ).'*mp.P3.full.dx;
                    %mp.P3.compact.xs = ( -(mp.P3.compact.Narr-1)/2:(mp.P3.compact.Narr-1)/2 ).'*mp.P3.compact.dx;
                else
                    % mp.P3.full.xs = (-mp.P3.full.Narr/2:(mp.P3.full.Narr/2-1)).'*mp.P3.full.dx;
                    %mp.P3.compact.xs = (-mp.P3.compact.Narr/2:(mp.P3.compact.Narr/2-1)).'*mp.P3.compact.dx;
                end
        
                if(mp.flagPlot)
                    figure(504); imagesc(padOrCropEven(mp.P3.full.mask,length(mp.P1.full.mask)) + mp.P1.full.mask); axis xy equal tight; colorbar;
                end
        end
        
        mp.P3.full.Narr= length(mp.P3.full.mask);
%         %%
%         fac = 1;%15/15.022; %.9985;
%         mp.P3.full.mask = falco_gen_multi_ring_SP(fac*mp.rEdgesLeft,fac*mp.rEdgesRight,mp.P2.full.dx,mp.P2.D,mp.centering);
%         mp.P3.full.Narr= length(mp.P3.full.mask);
%         A1 = mp.P3.full.mask; A1pad = padOrCropEven(A1,1000); 
%         A2 = fitsread('~/Downloads/out_RSPLC1D_A_maxTrPH_2848Dpup9960_33WA228_38LS83_10FPres5_BW10N11_c100_Nring10_2D_N1000even_sizedForFull.fits');
%         figure; imagesc(A2-A1pad); colorbar;
%         %%
%         keyboard
        %mp.P3.compact.Nbeam = mp.P1.compact.Narr;  %--Number of pixels across the array containing the SP pupil in the compact model  

        mp.P3.full.dx = mp.P2.full.dx;
        mp.P3.compact.dx = mp.P2.compact.dx;
        mp.P3.compact.Narr = length(mp.P3.compact.mask);

        
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
        
        if(nnz(strcmp(mp.whichPupil,{'LUVOIRA5','LUVOIR_B_offaxis','HabEx_B_offaxis'}))>0 && mp.flagApod)
            % Full aperture stop 
            mp.P3.full.Narr = 2^(nextpow2(mp.P1.full.Nbeam));
%             mp.P3.full.dx = mp.P2.full.dx;

            if(strcmp(mp.P3.apodType,'Simple'))
                inputs.Nbeam = mp.P1.full.Nbeam;     % number of points across incoming beam 
                inputs.Npad = mp.P3.full.Narr;% 
                inputs.OD = mp.P3.ODnorm;
                inputs.ID = mp.P3.IDnorm;
                inputs.num_strut = 0;
                inputs.strut_angs = [];%Angles of the struts 
                inputs.strut_width = 0;% spider width (fraction of the pupil diameter)

                mp.P3.full.mask= falco_gen_pupil_Simple( inputs );

% <<<<<<< thin_film
% %             mp.P3.full.mask = falco_gen_multi_ring_SP([0.20, 0.55]/2, [0.50, 0.95]/2, mp.P2.full.dx, mp.P2.D, mp.centering);
% %             mp.P3.full.mask = falco_gen_multi_ring_SP(mp.P3.IDvec, mp.P3.ODvec, mp.P2.full.dx, mp.P2.D, mp.centering);
%             mp.P3.full.mask = padOrCropEven(mp.P3.full.mask, mp.P3.full.Narr);
            
            
%             % Compact aperture stop 
%             inputs.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
%             inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));% 
            
%             mp.P3.compact.Narr = 2^(nextpow2(mp.P1.compact.Nbeam));
% % %            mp.P3.compact.dx = mp.P1.compact.dx;
%             mp.P3.compact.mask = falco_gen_pupil_Simple( inputs );
% %             mp.P3.compact.mask = falco_gen_multi_ring_SP(mp.P3.IDvec,mp.P3.ODvec, mp.P2.compact.dx, mp.P2.D, mp.centering);
%             mp.P3.compact.mask = padOrCropEven(mp.P3.compact.mask, mp.P3.compact.Narr);
% =======
                % Compact aperture stop 
                inputs.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
                inputs.Npad = 2^(nextpow2(mp.P1.compact.Nbeam));% 

                mp.P3.compact.Narr = 2^(nextpow2(mp.P1.compact.Nbeam));
    %            mp.P3.compact.dx = mp.P1.compact.dx;
                mp.P3.compact.mask = falco_gen_pupil_Simple( inputs );
            else
                load(mp.P3.apodType)
                mp.P3.full.mask = padOrCropEven(APOD,mp.P3.full.Narr);
                %mp.P3.full.mask = propcustom_2FT(mp.P3.full.mask);
                mp.P3.compact.Narr = mp.P3.full.Narr;
                mp.P3.compact.mask = mp.P3.full.mask;
            end
            
% >>>>>>> master
        else
            disp('Using vortex without apodizer or aperture stop.')
        end
	case{'APLC'}
        
        if(~strcmpi(mp.centering,'pixel'))
            error('Use pixel centering for APLC');
        end
        if(mp.P1.full.Nbeam~=mp.P1.compact.Nbeam)
            error('APLC currently requires Nbeam for the full and compact.');
        end

        mp.P3.full.mask = falco_gen_tradApodizer(mp.P1.full.mask,mp.P1.full.Nbeam,mp.F3.Rin,(1+mp.fracBW/2)*mp.F3.Rout,mp.useGPU);
        mp.P3.full.Narr = length(mp.P3.full.mask);
        
        mp.P3.compact.mask = mp.P3.full.mask;
        mp.P3.compact.Narr = length(mp.P3.compact.mask);
        
        mp.P3.full.dx = mp.P2.full.dx;
        mp.P3.compact.dx = mp.P2.compact.dx;
        
        if(mp.flagPlot)
            figure(504); imagesc(padOrCropEven(mp.P3.full.mask,length(mp.P1.full.mask)).*mp.P1.full.mask); axis xy equal tight; colorbar;
        end
        
%         %--Generate lower-resolution SP for the compact model
%         mp.P3.compact.mask = falco_gen_multi_ring_SP(mp.rEdgesLeft,mp.rEdgesRight,mp.P2.compact.dx,mp.P2.D,mp.centering);

    otherwise
        error([mp.coro,' is not a valid option for mp.coro.']);
end
% figure;
% imagesc(mp.P1.full.mask+mp.P3.full.mask);
% axis image;
% 
% asdf


end %--END OF FUNCTION



