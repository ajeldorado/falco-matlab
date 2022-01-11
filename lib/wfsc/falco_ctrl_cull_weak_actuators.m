% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Cull weak actuators from the set used for control.
%
% INPUTS
% ------
% mp : structure of model parameters
% cvar : structure containing variables for the controller
% jacStruct : structure containing the Jacobians
%
% OUTPUTS
% -------
% mp : structure of model parameters
% jacStruct : structure containing the Jacobians

function [mp, jacStruct] = falco_ctrl_cull_weak_actuators(mp, cvar, jacStruct)

    %% Initialization of Flag and Which Actuators

    %--If new actuators are added, perform a new cull of actuators.
    if cvar.Itr == 1
        cvar.flagCullAct = true;
    else
        if isfield(mp, 'dm_ind_sched')
            cvar.flagCullAct = ~isequal(mp.dm_ind_sched{cvar.Itr}, mp.dm_ind_sched{cvar.Itr-1});
        else
            cvar.flagCullAct = false;
        end
    end
    mp.flagCullActHist(cvar.Itr) = cvar.flagCullAct;

    %--Before performing new cull, include all actuators again
    if cvar.flagCullAct
        %--Re-include all actuators in the basis set. Need act_ele to be a column vector.
        if any(mp.dm_ind == 1); mp.dm1.act_ele = (1:mp.dm1.NactTotal).'; end
        if any(mp.dm_ind == 2); mp.dm2.act_ele = (1:mp.dm2.NactTotal).'; end
        if any(mp.dm_ind == 8); mp.dm8.act_ele = (1:mp.dm8.NactTotal).'; end
        if any(mp.dm_ind == 9); mp.dm9.act_ele = (1:mp.dm9.NactTotal).'; end
        %--Update the number of elements used per DM
        if any(mp.dm_ind == 1); mp.dm1.Nele = length(mp.dm1.act_ele); else; mp.dm1.Nele = 0; end
        if any(mp.dm_ind == 2); mp.dm2.Nele = length(mp.dm2.act_ele); else; mp.dm2.Nele = 0; end
        if any(mp.dm_ind == 8); mp.dm8.Nele = length(mp.dm8.act_ele); else; mp.dm8.Nele = 0; end
        if any(mp.dm_ind == 9); mp.dm9.Nele = length(mp.dm9.act_ele); else; mp.dm9.Nele = 0; end
    end
        
    %% Cull Weak Actuators
    %--Reduce the number of actuators used based on their relative strength in the Jacobian
    if cvar.flagCullAct && cvar.flagRelin
        fprintf('Weeding out weak actuators from the control Jacobian...\n'); 
        if any(mp.dm_ind == 1)
            %--Crop out very weak-effect actuators. Need mp.dm1.act_ele to be a column vector
            G1intNorm = sum( mean(abs(jacStruct.G1).^2,3), 1).';
            mp.dm1.act_ele = find(G1intNorm/max(G1intNorm(:))>=10^(mp.logGmin));
            clear G1intNorm
        end
        if any(mp.dm_ind == 2)
            G2intNorm = sum( mean(abs(jacStruct.G2).^2,3),1).';
            mp.dm2.act_ele = find(G2intNorm/max(G2intNorm(:))>=10^(mp.logGmin));
            clear G2intNorm
        end
        if any(mp.dm_ind == 8)
            G8intNorm = sum(mean(abs(jacStruct.G8).^2,3),1).';
            mp.dm8.act_ele = find(G8intNorm/max(G8intNorm(:))>=10^(mp.logGmin));
            clear G8intNorm
        end    
        if any(mp.dm_ind == 9)
            G9intNorm = sum(mean(abs(jacStruct.G9).^2,3),1).';
            mp.dm9.act_ele = find(G9intNorm/max(G9intNorm(:))>=10^(mp.logGmin));
            clear G9intNorm
        end

        %--Remove pinned or dead actuators from Jacobian
        if any(mp.dm_ind == 1)
            mp.dm1.act_ele = setdiff(mp.dm1.act_ele, mp.dm1.pinned);
        end
        if any(mp.dm_ind == 2)
            mp.dm2.act_ele = setdiff(mp.dm2.act_ele, mp.dm2.pinned);
        end

        %--Add back in all actuators that are tied (to keep the tied actuator logic correct when a weak actuator is included)
        if any(mp.dm_ind == 1)
           for ti=1:size(mp.dm1.tied,1)
               if(any(mp.dm1.act_ele==mp.dm1.tied(ti,1))==false);  mp.dm1.act_ele = [mp.dm1.act_ele; mp.dm1.tied(ti,1)];  end
               if(any(mp.dm1.act_ele==mp.dm1.tied(ti,2))==false);  mp.dm1.act_ele = [mp.dm1.act_ele; mp.dm1.tied(ti,2)];  end
           end
           mp.dm1.act_ele = sort(mp.dm1.act_ele); %--Need to sort for the logic in model_Jacobian.m
        end
        if any(mp.dm_ind == 2)
           for ti=1:size(mp.dm2.tied,1)
               if(any(mp.dm2.act_ele==mp.dm2.tied(ti,1))==false);  mp.dm2.act_ele = [mp.dm2.act_ele; mp.dm2.tied(ti,1)];  end
               if(any(mp.dm2.act_ele==mp.dm2.tied(ti,2))==false);  mp.dm2.act_ele = [mp.dm2.act_ele; mp.dm2.tied(ti,2)];  end
           end
           mp.dm2.act_ele = sort(mp.dm2.act_ele);
        end
        if any(mp.dm_ind == 8)
           for ti=1:size(mp.dm8.tied,1)
               if(any(mp.dm8.act_ele==mp.dm8.tied(ti,1))==false);  mp.dm8.act_ele = [mp.dm8.act_ele; mp.dm8.tied(ti,1)];  end
               if(any(mp.dm8.act_ele==mp.dm8.tied(ti,2))==false);  mp.dm8.act_ele = [mp.dm8.act_ele; mp.dm8.tied(ti,2)];  end
           end
           mp.dm8.act_ele = sort(mp.dm8.act_ele);
        end
        if any(mp.dm_ind == 9)
           for ti=1:size(mp.dm9.tied,1)
               if(any(mp.dm9.act_ele==mp.dm9.tied(ti,1))==false);  mp.dm9.act_ele = [mp.dm9.act_ele; mp.dm9.tied(ti,1)];  end
               if(any(mp.dm9.act_ele==mp.dm9.tied(ti,2))==false);  mp.dm9.act_ele = [mp.dm9.act_ele; mp.dm9.tied(ti,2)];  end
           end
           mp.dm9.act_ele = sort(mp.dm9.act_ele);
        end

        %--Update the number of elements used per DM
        if any(mp.dm_ind == 1); mp.dm1.Nele = length(mp.dm1.act_ele); end
        if any(mp.dm_ind == 2); mp.dm2.Nele = length(mp.dm2.act_ele); end
        if any(mp.dm_ind == 8); mp.dm8.Nele = length(mp.dm8.act_ele); end
        if any(mp.dm_ind == 9); mp.dm9.Nele = length(mp.dm9.act_ele); end
        %mp.NelePerDMvec = [length(mp.dm1.Nele), length(mp.dm2.Nele), length(mp.dm3.Nele), length(mp.dm4.Nele), length(mp.dm5.Nele), length(mp.dm6.Nele), length(mp.dm7.Nele), length(mp.dm8.Nele), length(mp.dm9.Nele) ];

        if any(mp.dm_ind == 1); fprintf('  DM1: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm1.Nele, mp.dm1.NactTotal,100*mp.dm1.Nele/mp.dm1.NactTotal); end
        if any(mp.dm_ind == 2); fprintf('  DM2: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm2.Nele, mp.dm2.NactTotal,100*mp.dm2.Nele/mp.dm2.NactTotal); end
        if any(mp.dm_ind == 8); fprintf('  DM8: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm8.Nele, mp.dm8.NactTotal,100*mp.dm8.Nele/mp.dm8.NactTotal); end
        if any(mp.dm_ind == 9); fprintf('  DM9: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm9.Nele, mp.dm9.NactTotal,100*mp.dm9.Nele/mp.dm9.NactTotal); end

        %--Crop out unused actuators from the control Jacobian
        if any(mp.dm_ind == 1); jacStruct.G1 = jacStruct.G1(:,mp.dm1.act_ele,:); end
        if any(mp.dm_ind == 2); jacStruct.G2 = jacStruct.G2(:,mp.dm2.act_ele,:); end
        if any(mp.dm_ind == 8); jacStruct.G8 = jacStruct.G8(:,mp.dm8.act_ele,:); end
        if any(mp.dm_ind == 9); jacStruct.G9 = jacStruct.G9(:,mp.dm9.act_ele,:); end
    end  

end %--END OF FUNCTION falco_ctrl_cull

