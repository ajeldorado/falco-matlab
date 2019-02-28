% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_jac_weights(mp)
%
% Function to set the weights for the Jacobian modes.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-09-17 by A.J. Riggs.




function mp = falco_config_jac_weights(mp)


%     % mp.jac.maxZnoll = 1;%3; %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include at least 1 for the on-axis piston mode.
%     % mp.jac.Zcoef = 1e-9*[1, 5, 5, 0, 1, 1]; %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
%     % if(isfield(mp.jac,'maxZnoll')==false);   mp.jac.maxZnoll = 1;  end  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include at least 1 for the on-axis piston mode.
    
    
    if(isfield(mp.jac,'zerns')==false);   mp.jac.zerns = 1;  end  %--Which Zernike modes to include in Jacobian. Given as a vector of Noll indices. 1 is the on-axis piston mode.
    mp.jac.Nzern = numel(mp.jac.zerns);

    indZern1 = find(mp.jac.zerns==1); %--Find index in vector of piston (Z1)
    mp.jac.Zcoef(indZern1) = 1; %[1, 1e-9, 1e-9]; % [meters] reset coefficient for piston to 1

    mp.jac.weightMat = zeros(mp.Nsbp, mp.jac.Nzern); %--Initialize weighting matrix of each Zernike-wavelength mode for the controller

    for izern = 1:mp.jac.Nzern
        whichZern = mp.jac.zerns(izern);
        
        if(whichZern==1)
            mp.jac.weightMat(:,1) = ones(mp.Nsbp,1);
        else
            mp.jac.weightMat(1,izern) = 1;
            mp.jac.weightMat(mp.si_ref,izern) = 1;
            mp.jac.weightMat(mp.Nsbp,izern) = 1;
        end
    end

    %--Half-weighting if endpoint wavelengths are used
    if(strcmpi(mp.estimator,'perfect')) %--For design or modeling without estimation: Choose ctrl wvls evenly between endpoints of the total bandpass
        mp.jac.weightMat(1,:) = 0.5*mp.jac.weightMat(1,:);
        mp.jac.weightMat(mp.Nsbp,:) = 0.5*mp.jac.weightMat(mp.Nsbp,:);
    end

    %--Normalize the summed weights of each column separately
    for izern = 1:mp.jac.Nzern
        colSum = sum( mp.jac.weightMat(:,izern) );
        mp.jac.weightMat(:,izern) = (1/colSum)*mp.jac.weightMat(:,izern);
    end

    %--Zero out columns for which the RMS Zernike value is zero
    for izern = 1:mp.jac.Nzern
        if(mp.jac.Zcoef(izern)==0)
            mp.jac.weightMat(:,izern) = 0*mp.jac.weightMat(:,izern);
        end
    end

    mp.jac.weightMat_ele = find(mp.jac.weightMat>0); %--Indices of the non-zero control Jacobian modes in the weighting matrix
    mp.jac.weights = mp.jac.weightMat(mp.jac.weightMat_ele); %--Vector of control Jacobian mode weights
    mp.jac.Nmode = length(mp.jac.weights); %--Number of (Zernike-wavelength pair) modes in the control Jacobian

    %--Get the wavelength indices for the nonzero values in the weight matrix. 
    temp = (1:mp.Nsbp).';
    tempMat = repmat(temp,[1,mp.jac.Nzern]);
    mp.jac.sbp_inds = tempMat(mp.jac.weightMat_ele);
    % mp.Wttlam_si = tempMat(mp.Wttlam_ele);

    %--Get the Zernike indices for the nonzero elements in the weight matrix. 
    temp = mp.jac.zerns; %1:mp.jac.maxZnoll;
    tempMat = repmat(temp,[mp.Nsbp,1]);
    mp.jac.zern_inds = tempMat(mp.jac.weightMat_ele);
    % mp.Wttlam_ti = tempMat(mp.Wttlam_ele);


end %--END OF FUNCTION




%% Previous Method (tip/tilt only; no other Zernikes):

%% Tip/Tilt and Spatial Weighting of the Control Jacobian  #NEWFORTIPTILT
% mp.ti_ref = ceil(mp.Ntt/2);
% mp.mas2lam0D = 1/(mp.lambda0/mp.P1.D*180/pi*3600*1000); %--Conversion factor: milliarcseconds (mas) to lambda0/D
% %--Define the (x,y) values for each tip/tilt offset in units of lambda0/D
% if(mp.Ntt == 5)
%     mp.ttx = mp.mas2lam0D*mp.TToffset*[0,cosd(0),cosd(90),cosd(180),cosd(270)];
%     mp.tty = mp.mas2lam0D*mp.TToffset*[0,sind(0),sind(90),sind(180),sind(270)];
% elseif(mp.Ntt == 4)
%     mp.ttx = mp.mas2lam0D*mp.TToffset*[0,cosd(0),cosd(120),cosd(240)];
%     mp.tty = mp.mas2lam0D*mp.TToffset*[0,sind(0),sind(120),sind(240)];
% elseif(mp.Ntt == 1)
%     mp.ttx = 0;
%     mp.tty = 0;
% else
%     disp('######## ERROR: Number of tip-tilt modes not specified properly. ##########')
%     return    
% end
% 
% mp.Wttlam = zeros(mp.Nsbp,mp.Ntt); %--Initialize weighting matrix of each tip/tilt-wavelength mode for the controller
% if(isinf(mp.NlamForTT))
%     mp.Wttlam = ones(mp.Nsbp,mp.Ntt); %--Full usage and equal weighting for all T/T's and sub-bandpasses.
% elseif(mp.NlamForTT==3)
%     %--Set tip/tilt offsets at the middle and both end sub-bandpasses.
%     mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
%     mp.Wttlam(1,:) = ones(1,mp.Ntt);
%     mp.Wttlam(mp.si_ref,:) = ones(1,mp.Ntt);
%     mp.Wttlam(end,:) = ones(1,mp.Ntt);
% elseif(mp.NlamForTT==2)
%     %--Set tip/tilt offsets at only both end sub-bandpasses.
%     mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
%     mp.Wttlam(1,:) = ones(1,mp.Ntt);
%     mp.Wttlam(end,:) = ones(1,mp.Ntt);
% elseif(mp.NlamForTT==1)
%     %--Set tip/tilt offsets at only the middle sub-bandpass.
%     mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
%     mp.Wttlam(mp.si_ref,:) = ones(1,mp.Ntt);
% elseif(mp.NlamForTT==0)
%     %--Set tip/tilt offsets at no sub-bandpasses.
%     mp.Wttlam(:,1) = mp.Ntt*ones(mp.Nsbp,1);
% end
% 
% mp.Wsum = sum(sum(mp.Wttlam)); %--Sum of all the control Jacobian weights
% 
% mp.Wttlam_ele = find(mp.Wttlam>0); %--Indices of the non-zero control Jacobian modes in the weighting matrix
% mp.jac.weights = mp.Wttlam(mp.Wttlam_ele); %--Vector of control Jacobian mode weights
% mp.jac.Nmode = length(mp.Wttlam_ele); %--Number of modes in the control Jacobian
% 
% %--Get the wavelength indices for the nonzero values in the weight matrix. 
% temp = (1:mp.Nsbp).';
% tempMat = repmat(temp,[1,mp.Ntt]);
% mp.Wttlam_si = tempMat(mp.Wttlam_ele);
% 
% %--Get the tip/tilt indices for the nonzero values in the weight matrix. 
% temp = 1:mp.Ntt;
% tempMat = repmat(temp,[mp.Nsbp,1]);
% mp.Wttlam_ti = tempMat(mp.Wttlam_ele);

