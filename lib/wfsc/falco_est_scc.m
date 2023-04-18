% Copyright 2018, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute the electric estimate using a self-coherent camera (SCC).
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% ev : structure of estimation variables

function ev = falco_est_scc(mp)

%     fprintf('Estimating electric field with the SCC ...\n'); tic;
    %mp.flagPlot = true;

    ev.Eest = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
    ev.Eest_full = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);
    ev.IincoEst = zeros(mp.Fend.corr.Npix, mp.jac.Nmode);
    ev.Im = zeros(mp.Fend.Neta, mp.Fend.Nxi);
    ev.imageArray = zeros(mp.Fend.Neta, mp.Fend.Nxi, 1, mp.Nsbp);
    imageShape0 = [mp.Fend.Neta, mp.Fend.Nxi];
    imageShapeMax = max(imageShape0);

    for iStar = 1:mp.compact.star.count

        modvar = ModelVariables;
        modvar.starIndex = iStar;
        modvar.whichSource = 'star';

        for iSubband = 1:mp.Nsbp

            modvar.sbpIndex = iSubband;
            modeIndex = (iStar-1)*mp.Nsbp + iSubband;
%             fprintf('Wavelength: %u/%u,    ', iSubband, mp.Nsbp);
%             fprintf('Mode: %u/%u \n ', modeIndex, mp.jac.Nmode);  

            image0 = falco_get_sbp_image(mp, iSubband);
            image0Vec = image0(mp.Fend.corr.maskBool);

            if iStar == 1 % Image already includes all stars, so don't sum over star loop
                ev.imageArray(:, :, 1, iSubband) = image0;
                ev.Im = ev.Im + mp.sbp_weights(iSubband)*image0; % subband-averaged image for plotting
            end

            % Consider replacing the butterworth filter with just
            % convolving the dark hole software mask with a soft-edge
            % kernel.
            
            imageShapeMax = max([mp.Fend.Neta, mp.Fend.Nxi]);
            ctr_x = 0.;
            ctr_y = -0.5*( mp.Fend.corr.Rout + mp.Fend.corr.Rin)*mp.Fend.res;
            Y0 = (mp.Fend.corr.Rout - mp.Fend.corr.Rin +2)*mp.Fend.res/2;
            X0 = (mp.Fend.corr.Rout+2)*mp.Fend.res;
            expx = 4;%3.8;
            expy = 7.9; %7.6;%2.8;
            
            expxcirc = 3.3;%5;
            %filter_butterworth =falco_butterworth(imageShapeMax, imageShapeMax, ctr_x, ctr_y, expx,expy, X0, Y0);
            filter_butterworth = falco_butterworth(imageShapeMax, imageShapeMax, 0, -500, expx,expy, 800, 375, 'squared').*...
                falco_butterworth(imageShapeMax, imageShapeMax, 0, 0, expxcirc,expxcirc, 325, 325, 'circle');

            %filter_butterworth = (1-falco_butterworth(imageShapeMax, imageShapeMax, mp.scc.butterworth_exponent, mp.Fend.corr.Rin*mp.Fend.res)) .* ...
             %   falco_butterworth(imageShapeMax, imageShapeMax, 0.5, mp.Fend.corr.Rout*mp.Fend.res);
                %falco_butterworth(imageShapeMax, imageShapeMax, mp.scc.butterworth_exponent, mp.Fend.corr.Rout*mp.Fend.res);
            
            imageFilt = pad_crop(image0, [imageShapeMax, imageShapeMax]);
            imageFilt = imageFilt .* filter_butterworth; % The filter makes it harder to find the side lobe center.


            imageFilt = ifftshift(imageFilt);
            pupil_scc = fft2(imageFilt);%/imageShapeMax;
            pupil_scc = fftshift(pupil_scc);
            %pupil = pad_crop(pupil, SCCpupsubwindowsize);
            
            if mp.flagPlot
                figure(1001); imagesc(log10(abs(image0))); axis xy equal tight; colorbar; drawnow;
                %figure(1002); imagesc(log10(abs(fftshift(imageFilt)))); axis xy equal tight; colorbar; drawnow;

                %figure(1003); imagesc(log10(abs(pupil_scc))); axis xy equal tight; colorbar; drawnow;
                %figure(1004); imagesc(angle(pupil_scc)); axis xy equal tight; colorbar; drawnow;
            end
        
            % Make template
            
%             size_central_peak = 2*mp.tb.sciCam.subwindowsize/mp.Fend.res/1.22;
%             pix_per_mm = size_central_peak/2/scc.diameter_lyot;
%             size_lateral_peak = (scc.diameter_lyot+scc.diameter_pinhole)*pix_per_mm;
%             size_uncertainty = size_lateral_peak/scc.uncertainty_lateral_peak;
%             distance_lateral_peak = scc.distance_pinhole*pix_per_mm;
%             x0_lateral_peak = mp.tb.sciCam.subwindowsize/2 + distance_lateral_peak*cos(scc.theta_0*pi/180);
%             y0_lateral_peak = mp.tb.sciCam.subwindowsize/2 + distance_lateral_peak*sin(scc.theta_0*pi/180);
%             pupil_scc = offcenter_crop(pupil_scc, y0_lateral_peak, x0_lateral_peak, scc.pupsubwindowsize, scc.pupsubwindowsize);
%             inputs.Nbeam = size_lateral_peak + size_uncertainty; % number of points across usable pupil   

            pupil_scc = offcenter_crop(pupil_scc, mp.scc.pupil_center_row, mp.scc.pupil_center_col,...
                imageShapeMax, imageShapeMax);
%                 mp.scc.pupil_subwindow_array_width, mp.scc.pupil_subwindow_array_width);
            inputs.OD = 1;
            inputs.Nbeam = mp.scc.pupil_subwindow_diameter;
            inputs.Npad = imageShapeMax; %mp.scc.pupil_subwindow_array_width;
            window2D = falco_gen_pupil_Simple(inputs);
            pupil_scc = pupil_scc .* window2D;

            if mp.flagPlot
                figure(1018); imagesc(pad_crop(power(abs(pupil_scc), 0.25),100)); axis xy equal tight;
            end

            pupil_scc = fftshift(pupil_scc);
            estimate_scc = ifft2(pupil_scc);%/imageShapeMax;
            estimate_scc = ifftshift(estimate_scc);
            estimate_scc = pad_crop(estimate_scc, imageShape0);%./filter_butterworth;
            estimate_scc = -1j*(estimate_scc); % Not necessary, but makes closer to "true" estimate
            Eest = estimate_scc(mp.Fend.corr.maskBool);

            if mp.flagPlot
                figure(1005); imagesc(log10(abs(estimate_scc))); axis xy equal tight; colorbar; drawnow;
                %figure(1006); imagesc(angle(estimate_scc)); axis xy equal tight; colorbar; drawnow;
            end
            
            ev.Eest(:, modeIndex) = Eest;
            ev.Eest_full(:, :, 1, modeIndex) = estimate_scc;
            
            % NOTE: Incoherent estimate is incorrect because of the scale
            % factor between the SCC estimate and the true E-field.
            ev.IincoEst(:, modeIndex) =  image0Vec - abs(Eest).^2; % incoherent light
            
        end
    end
    
end %--END OF FUNCTION
