% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%--Function to compute a stack of normalized, 2-D Zernike mode maps.
%
%--Requires the PROPER library and the modified PROPER function
% propcustom_zernikes.m, which allows for interpixel centering.
%
%--REQUIRED INPUTS
% Nbeam: Number of points across the beam diameter
% centering: pixel centering of the beam ('pixel' or 'interpixel').
% indsZnoll: vector containing noll indices of Zernikes to compute maps for
%
%--OUTPUTS
%  ZmapCube: Datacube of normalized Zernike maps
%  
%
%--VERSION CHANGE HISTORY
% Written on 2018-03-03 by A.J. Riggs.

function ZmapCube = falco_gen_norm_zernike_maps(Nbeam,centering,indsZnoll)

    %--Set array size as minimum width to contain the beam.
    if(strcmpi(centering,'interpixel') )
        Narray = Nbeam; %--No zero-padding needed if beam is centered between pixels
    else
        Narray = Nbeam + 2; %--number of points across output array. Requires two more pixels when pixel centered.
    end

    %--PROPER setup values
    Dbeam = 1;                 %--Diameter of aperture, normalized to itself
    wl   = 1e-6;              % wavelength (m); Dummy value--no propagation here, so not used.
    beamdiamfrac = Narray/Nbeam;

    %--Initialize wavefront structure in PROPER
    bm = prop_begin(Dbeam, wl, Narray,'beam_diam_fraction',beamdiamfrac);
    % figure(2); imagesc(abs(bm.wf)); axis xy equal tight; colorbar;

    %--Use modified PROPER function to generate the Zernike datacube
    Nzern = length(indsZnoll);
    ZmapCube = zeros(Narray,Narray,Nzern);
    
    bm.centering = centering;
    for iz = 1:Nzern
        [~, ZmapCube(:,:,iz)] = propcustom_zernikes( bm, indsZnoll(iz), 1 );
    end

    % %--DEBUGGING: Verify centering is correct.
    % Zmap1 = ZmapCube(:,:,1);
    % if(strcmpi(centering,'pixel'))
    %     Zmap1 = Zmap1(2:end,2:end);
    % end
    % figure(1); imagesc(Zmap1); axis xy equal tight; colorbar;
    % figure(11); imagesc(Zmap1-rot90(Zmap1,2)); axis xy equal tight; colorbar;



end %--END OF FUNCTION
