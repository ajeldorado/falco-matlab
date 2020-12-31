% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_set_jac_weights(mp)
%
% Function to set the weights for the Jacobian modes.
%
% REVISION HISTORY:
% ----------------
% Created on 2018-09-17 by A.J. Riggs.

function mp = falco_set_jacobian_weights(mp)
    
    %--Which Zernike modes to include in Jacobian. Given as a vector of Noll indices. 1 is the on-axis piston mode.
    mp.jac.Nzern = numel(mp.jac.zerns);
    if(length(mp.jac.zerns) ~= length(unique(mp.jac.zerns)))
       error('The vector of Zernike modes to control cannot contain repeated values.') 
    end
    if(any(mp.jac.zerns==1) == false)
        error('Piston must be included as a controlled Zernike.')
    else
        mp.jac.Zcoef(mp.jac.zerns==1) = 1; % [meters] reset coefficient for piston to 1
    end
    mp.jac.weightMat = zeros(mp.Nsbp, mp.jac.Nzern); %--Initialize weighting matrix of each Zernike-wavelength mode for the controller

    for izern = 1:mp.jac.Nzern
        whichZern = mp.jac.zerns(izern);
        if(whichZern==1)
            mp.jac.weightMat(:, 1) = ones(mp.Nsbp, 1);
        else
            mp.jac.weightMat(1, izern) = 1;
            mp.jac.weightMat(mp.si_ref, izern) = 1;
            mp.jac.weightMat(mp.Nsbp, izern) = 1;
        end
    end

    %--Half-weighting if endpoint wavelengths are used
    if(strcmpi(mp.estimator,'perfect')) %--For design or modeling without estimation: Choose ctrl wvls evenly between endpoints of the total bandpass
        mp.jac.weightMat(1, :) = 0.5*mp.jac.weightMat(1, :);
        mp.jac.weightMat(mp.Nsbp, :) = 0.5*mp.jac.weightMat(mp.Nsbp, :);
    end

    %--Normalize the summed weights of each column separately
    for izern = 1:mp.jac.Nzern
        colSum = sum(mp.jac.weightMat(:, izern));
        mp.jac.weightMat(:, izern) = (1/colSum)*mp.jac.weightMat(:, izern);
    end

    %--Zero out columns for which the RMS Zernike value is zero
    for izern = 1:mp.jac.Nzern
        if(mp.jac.Zcoef(izern) == 0)
            mp.jac.weightMat(:, izern) = 0*mp.jac.weightMat(:, izern);
        end
    end

    mp.jac.weightMatInd = find(mp.jac.weightMat > 0); %--Indices of the non-zero control Jacobian modes in the weighting matrix
    NmodePerStar = length(mp.jac.weightMatInd);
    mp.jac.weights = repmat(mp.jac.weightMat(mp.jac.weightMatInd), [1, mp.compact.star.count]); %--Vector of control Jacobian mode weights
    mp.jac.Nmode = mp.compact.star.count * NmodePerStar; %--Number of (Zernike-wavelength pair-star) modes in the control Jacobian
    
    %--Get the wavelength indices for the nonzero values in the weight matrix. 
    temp = (1:mp.Nsbp).';
    tempMat = repmat(temp, [1, mp.jac.Nzern]);
    mp.jac.sbp_inds = repmat(tempMat(mp.jac.weightMatInd), [1, mp.compact.star.count]);
    
    %--Get the Zernike indices for the nonzero elements in the weight matrix. 
    temp = mp.jac.zerns;
    tempMat = repmat(temp, [mp.Nsbp, 1]);
    mp.jac.zern_inds = repmat(tempMat(mp.jac.weightMatInd), [1, mp.compact.star.count]);
    
    %--Get the star indices for each mode
    mp.jac.star_inds = zeros(1, mp.jac.Nmode);
    for iStar = 1:mp.compact.star.count
        mp.jac.star_inds((iStar-1)*NmodePerStar+1:iStar*NmodePerStar) = iStar*ones(1, NmodePerStar);
    end

end %--END OF FUNCTION