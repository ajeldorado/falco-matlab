% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ EPsupport, apRad ] = falco_gen_pupil_SCDA( apType, spiderType, spiderSize, Narray, centering )
%falco_gen_pupil_SCDA Generates the SCDA pupil 
%   Detailed explanation goes here

    apDir = ['util',filesep,'aperture_files',filesep];

    apDia = '2000';

    primary = fitsread([apDir,apType,'_',apDia,'pix.fits']); 

    if(strcmp(spiderType,'none'))
        EPsupport = primary;
    else
        spiders = fitsread([apDir,spiderType,'_',apDia,'pix_',spiderSize,'cm.fits']);  
        EPsupport = primary.*spiders;
    end

    [Narray0,~] = size(primary); % size of the file 
    apRad = (Narray/Narray0)*str2num(apDia)/2; % Aperture radius (units of samples)

    % Downsample, also correct for inter-pixel centering
    if( strcmp(centering,'interpixel') || strcmp(centering,'even') )
        [X1,Y1] = meshgrid(-Narray0/2:Narray0/2-1); % Grids with Cartesian (x,y) coordinates 
        [X2,Y2] = meshgrid(linspace(-Narray0/2+0.5,Narray0/2-0.5,Narray));
        EPsupport = interp2(X1,Y1,EPsupport,X2,Y2,'linear',0);
    else
        [X1,Y1] = meshgrid(-Narray0/2:Narray0/2-1); % Grids with Cartesian (x,y) coordinates 
        [X2,Y2] = meshgrid(linspace(-Narray0/2,Narray0/2-1,Narray));
        EPsupport = interp2(X1,Y1,EPsupport,X2,Y2,'linear',0);
    end
    
end

