% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Inputs: 
% Nbeam - Beam diameter in samples
% Npad  - Size of the full, padded array (NpadxNpad)

function PUPIL = falco_gen_pupil_dst_LUVOIR_B( Nbeam, Npad )
%PUPIL = falco_gen_pupil_dst_LUVOIR_B( Nbeam, Npad )
%   Generates the DST LUVOIR B mask  

    flnm = 'dst_LUVB_seg100um.fits';
    
    % Load the high-res pupil mask file. Unzip it if neccessary.  
    try
        in = fitsread(flnm);
    catch
        unzip([flnm,'.zip']);
        try
            in = fitsread(flnm);
        catch
            error('Cant find the DST LUVOIR B pupil mask file.');
        end
    end
    
    Nbeam0 = 2*2084; % beam diameter in the saved file 
    
    magn = Nbeam/Nbeam0; % Magnification factor 
    
    % Scale the beam by the magnification factor 
    N = length(in);
    [X,Y] = meshgrid(-N/2:N/2-1);
    out = interp2(X,Y,in,X/magn,Y/magn,'linear',0);
   
    % Pad or crop to Npad
    PUPIL = padOrCropEven(out,Npad);
    
end