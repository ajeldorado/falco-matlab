% hcst_polarizationContrastDiffWvls.m
%
% Compute polarization contrast, i.e. ratio of light between going thru the
% VVC and the fused silica, as a function of wavelength
%
% Jorge Llop - DEc 23, 2020

%% Check polarization contrast, i.e. the ratio between going thru the empty slot and going thru the VVC


% Go through the Zernike mask
hcst_FPM_move(bench,[24.5 , 3, nan]);

tint = 1e-2;
hcst_andor_setExposureTime(bench,tint);
im0 = hcst_andor_getImage(bench) - dark0;
peak_crossed = max(im0(:))/tint;

% filter position
hcst_FW_setPos(bench,4);% Focus viewing mode
NDfilter_cal = 21.2938 * 4.4506;

% Move to VVC mask
hcst_FPM_move(bench,[1 ,1, bench.FPM.VORTEX_F0]);
indtint = 2;
while 1
    tint = tint_arr2(indtint);
    hcst_andor_setExposureTime(bench,tint);

    im0 = hcst_andor_getImage(bench) - dark0;

    peak = max(im0(:));

    figure(112);
    imagesc(log10_4plot(im0/2^16));
    axis image; 
    colorbar;%caxis([-4 0]);
%     title([num2str(count)],'FontSize',32)
    text(10,3,num2str(peak),'Color','w','FontSize',32);
    drawnow;
    if peak>60e3
        if indtint>1
            indtint = indtint-1;
        else
            break
        end
        disp('Exposure too high!')
    elseif peak<500
        if indtint<numtint
            indtint = indtint+1;
        end
        disp('Exposure too low!')
    else
        peak_vvc = peak/tint;
        break
    end
end

polContrast = peak_crossed/(peak_vvc*NDfilter_cal);

disp(['Polarization contrast: ',num2str(polContrast)])
