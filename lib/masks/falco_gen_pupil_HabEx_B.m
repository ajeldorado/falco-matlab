% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [mask] = falco_gen_pupil_HabEx_B(inputs)
%falco_gen_pupil_HabEx_B Summary of this function goes here
%   Detailed explanation goes here

    apRad = inputs.Nbeam/2;
    N = inputs.Npad;
    aperture_num = inputs.aperture_num;
    gap_size = inputs.gap_size;
    
    [X2,Y2] = meshgrid(-N/2:N/2-1);

    % rough guess at center
    centerrow = 1887;
    centercol = 1983;

    flnm = ['HabExB/OAP',num2str(aperture_num),'_',num2str(gap_size),'mmGap'];
    if(inputs.fivemmRadii)
        flnm = [flnm,'_5mmRadii'];
    end
        
	mask = imread([flnm,'.png']);
    
    [rows,cols,dims] = size(mask);
    [X1,Y1] = meshgrid(-cols/2:cols/2-1,-rows/2:rows/2-1);

    mask = mean(mask,3)/255;
    mask = 1 - mask;

    if(aperture_num==1)
        if(gap_size == 6)
            centerrow = centerrow+71;
            centercol = centercol+4;
        elseif(gap_size ==10)
            if(inputs.fivemmRadii)
                centerrow = centerrow+71.5;
                centercol = centercol+3.5;
            else
                centerrow = centerrow+48.6;
                centercol = centercol+8.2;
            end
        elseif(gap_size ==25)
            centerrow = centerrow+46.9;
            centercol = centercol+10.9;
        elseif(gap_size ==50)
            centerrow = centerrow+86;
            centercol = centercol+7.1;
        end
        scal_fac = apRad/3779*2;
    elseif(aperture_num==2)
        if(gap_size ==10)
            centerrow = centerrow+53;
            centercol = centercol-9.5;
        elseif(gap_size ==25)
            centerrow = centerrow+52;
            centercol = centercol-1;
        elseif(gap_size ==50)
            centerrow = centerrow+52.5;
            centercol = centercol-8.5;
        end
        scal_fac = apRad/3783*2;
    end

    shiftrows = rows/2+1 - centerrow;
    shiftcols = cols/2+1 - centercol;
    
    shiftrows = shiftrows/1000*apRad;
    shiftcols = shiftcols/1000*apRad;
    
    scal_fac = scal_fac*0.98;
    
    mask = interp2(X1,Y1,mask,(X2-shiftcols)/scal_fac,(Y2-shiftrows)/scal_fac,'linear',0);
end

