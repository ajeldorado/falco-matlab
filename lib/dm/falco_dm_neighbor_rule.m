% Copyright 2015,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to find neighboring actuators that exceed a specified difference
% in voltage and to scale down those voltages until the rule is met.
%
% INPUTS:
% - Vin = 2-D map of initial DM voltage commands
% - Vlim = maximum difference in voltages between neighbors
%
% OUTPUTS
% - Vout = 2-D map of adjusted DM voltage commands
% - indn = indices (act numbers) of the above
% - delv = values of differential voltages of the actuators violating the neighboring rule
%
% REVISION HISTORY
% - Created on 2015-11-23 by Erkin Sidick.
% - Modified on 2019-03-06 by A.J. Riggs from fun_neighbor_fix.m.
% -----------------------------------------------------------------------

% function [Vout, indn, delv] = falco_dm_neighbor_rule(Vin, Vlim, Nact)
function [Vout, indPair] = falco_dm_neighbor_rule(Vin, Vlim, Nact)

    Vout   = Vin; %--Initialize output voltage map
    indPair  = []; %--Initialize the paired indices list. [Npairs x 2]
    
    kx1 = [0 1; 1 1; 1 0];              % R1-C1
    kx2 = [0 1; 1 1; 1 0; 1 -1];        % R1, C2 - C47
    kx3 = [1 0; 1 -1];                  % R1, C48
    kx4 = [-1 1; 0 1; 1 1; 1 0];        % R2-R47, C1
    kx5 = [-1 1; 0 1; 1 1; 1 0; 1 -1];  % R2-47, C2-47
    kx6 = [1 0; 1 -1];                  % R2-47, C8
    kx7 = [-1 1; 0 1];                  % R48, C1 - C47
    kx8 = [-1 -1];                      % R48, C48
    
for jj = 1:Nact           % Row
    for ii = 1:Nact       % Col
                
        if jj ==1 
            if ii == 1
                kx = kx1;
            elseif ii < Nact
                kx = kx2;
            else
                kx = kx3;
            end
        elseif jj < Nact 
            if ii == 1
                kx = kx4;
            elseif ii < Nact
                kx = kx5;
            else
                kx = kx6;
            end
        else
            if ii < Nact
                kx = kx7;
            else
                kx = kx8;
            end
        end

        kr = jj + kx(:,1);
        kc = ii + kx(:,2);
        nNbr = length(kr); %--Number of neighbors
                
        if nNbr >= 1
            for iNbr = 1:nNbr
                
                a1 = Vout(jj,ii) - Vout(kr(iNbr),kc(iNbr)); %--Compute the delta voltage
                
                if abs(a1) > Vlim %--If neighbor rule is violated
                    
                    indLinCtr = (ii-1)*Nact + jj; %--linear index of center actuator
                    indLinNbr = (kc(iNbr)-1)*Nact + kr(iNbr); %--linear index of neigboring actuator
                    indPair = [indPair; indLinCtr,indLinNbr];

                    fx = (abs(a1) - Vlim) / 2;
                    Vout(jj,ii) = Vout(jj,ii) - sign(a1)*fx;
                    Vout(kr(iNbr),kc(iNbr)) = Vout(kr(iNbr),kc(iNbr)) + sign(a1)*fx;
                end
            end
        end
    end
end
        
return