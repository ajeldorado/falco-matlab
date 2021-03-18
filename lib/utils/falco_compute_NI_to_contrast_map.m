% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% conversionMap = falco_compute_NI_to_contrast_map(mp, Ndim)
%
% Compute the 2-D conversion map from normalized intensity to contrast.
%
% INPUTS
% ------
% mp : structure of all model parameters
% Ndim : number of axes from which to compute the 2-D map. Must be 1 or 2.
%
% OUTPUTS
% -------
% conversionMap : 2-D map the size of the image array

function conversionMap = falco_compute_NI_to_contrast_map(mp, Ndim)

    conversionMap = zeros(mp.Fend.Neta, mp.Fend.Nxi);

    % Define points at which offset the star
    if Ndim == 1

    %     maxSep = max(mp.Fend.RHOS(:));
        maxSep = max([abs(mp.Fend.xisDL(:)); abs(mp.Fend.etasDL(:))]);
        coords = 0:(1/(mp.Fend.res)):maxSep;

        xis = coords;
        etas = coords;
        rhos = sqrt(2)*coords;

    %     if mp.Fend.Nxi > mp.Fend.Neta
    % %         coords = mp.Fend.xisDL;
    % %         coords = coords(coords >= 0);
    %         xis = coords;
    %         etas = zeros(size(coords));
    %     else
    % %         coords = mp.Fend.etasDL;
    % %         coords = coords(coords >= 0);
    %         xis = zeros(size(coords));
    %         etas = coords;
    %     end

    elseif Ndim ==2

        [XIS, ETAS] = meshgrid(mp.Fend.xisDL, mp.Fend.etasDL);
        xis = XIS(mp.Fend.corr.maskBool);
        etas = ETAS(mp.Fend.corr.maskBool);

    end

    xiOffsetVec0 = mp.star.xiOffsetVec; % used for FALCO full models [lambda0/D]
    etaOffsetVec0 = mp.star.etaOffsetVec;  % used for FALCO full models [lambda0/D]
    Npoints = length(xis);
    peakVec = zeros(Npoints, 1);

    tic; fprintf('Computing conversion from NI to contrast...')
    for ip = 1:Npoints

        mp.star.xiOffsetVec = xiOffsetVec0 + xis(ip); % used for FALCO full models [lambda0/D]
        mp.star.etaOffsetVec = etaOffsetVec0 + etas(ip);
        if mp.full.flagPROPER
            mp.full.xoffset = xis(ip); % [lambda0/D]
            mp.full.yoffset = etas(ip); % [lambda0/D]
        end
        
        Itemp = falco_get_summed_image(mp);
        peakVec(ip) = max(Itemp(:));

    end
    fprintf('done. Time = %.2f s\n', toc)

    % Reset
    mp.star.xiOffsetVec = xiOffsetVec0; % used for FALCO full models [lambda0/D]
    mp.star.etaOffsetVec = etaOffsetVec0;  % used for FALCO full models [lambda0/D]
    if mp.full.flagPROPER
        mp.full.xoffset = 0; % [lambda0/D]
        mp.full.yoffset = 0; % [lambda0/D]
    end

    if Ndim == 2
        conversionMap(mp.Fend.corr.maskBool) = peakVec;
    elseif Ndim == 1
     
        conversionMap = interp1(rhos, peakVec, mp.Fend.RHOS, 'spline', 0);
    end
end
