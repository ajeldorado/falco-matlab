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
% - indn = 
% - delv = 
%
% REVISION HISTORY
% - Created on 2015-11-23 by Erkin Sidick.
% - Modified on 2019-03-06 by A.J. Riggs from fun_neighbor_fix.m.
% -----------------------------------------------------------------------

function [Vout, indn, delv] = falco_dm_neighbor_rule(Vin, Vlim, Nact)

    indn = [];
    Vout   = Vin;
    
    delv = [];
    
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
        ns = length(kr);
        
        %cont = [jj ii]
        
        if ns >= 1
        for j1 = 1:ns
            
            a1 = Vout(jj,ii) - Vout(kr(j1),kc(j1));
            
            indn = [indn; jj ii kr(j1) kc(j1)];
            delv = [delv a1];

            if abs(a1) > Vlim
                fx = (abs(a1) - Vlim) / 2;
                Vout(jj,ii) = Vout(jj,ii) - sign(a1)*fx;
                Vout(kr(j1),kc(j1)) = Vout(kr(j1),kc(j1)) + sign(a1)*fx;
            end
        end
        end
    end
end
        
return
