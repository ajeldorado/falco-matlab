% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to cull weak actuators from the set of used ones.
%
% INPUTS:
% -mp = structure of model parameters
% -cvar = structure containing variables for the controller
% -jacStruct = structure containing the Jacobians
%
% OUTPUTS:
% -mp = structure of model parameters
% -jacStruct = structure containing the Jacobians
%
% REVISION HISTORY:
% --------------
% Moved on 2019-09-30 by A.J. Riggs to be a separate file. Had tried
%   nesting in falco_wfsc_loop to save RAM but that did not work.

function [mp,jacStruct] = falco_ctrl_cull(mp,cvar,jacStruct)

        %% Cull Weak Actuators
        %--MOVE TO A FUNCTION
        %--Reduce the number of actuators used based on their relative strength in the Jacobian
        if(cvar.flagCullAct && cvar.flagRelin)
            fprintf('Weeding out weak actuators from the control Jacobian...\n'); 
            if(any(mp.dm_ind==1))
                %--Crop out very weak-effect actuators. Need mp.dm1.act_ele to be a column vector
                G1intNorm = sum( mean(abs(jacStruct.G1).^2,3), 1).';
                mp.dm1.act_ele = find(G1intNorm/max(G1intNorm(:))>=10^(mp.logGmin));
                clear G1intNorm
            end
            if(any(mp.dm_ind==2))
                G2intNorm = sum( mean(abs(jacStruct.G2).^2,3),1).';
                mp.dm2.act_ele = find(G2intNorm/max(G2intNorm(:))>=10^(mp.logGmin));
                clear G2intNorm
            end
            if(any(mp.dm_ind==8))
                G8intNorm(1:end) = sum(mean(abs(jacStruct.G8).^2,3),1).';
                mp.dm8.act_ele = find(G8intNorm/max(G8intNorm(:))>=10^(mp.logGmin));
                clear G8intNorm
            end    
            if(any(mp.dm_ind==9))
                G9intNorm = sum(mean(abs(jacStruct.G9).^2,3),1).';
                mp.dm9.act_ele = find(G9intNorm/max(G9intNorm(:))>=10^(mp.logGmin));
                clear G9intNorm
            end

           %--Add back in all actuators that are tied (to keep the tied actuator logic correct when a weak actuator is included)
           if(any(mp.dm_ind==1))
               for ti=1:size(mp.dm1.tied,1)
                   if(any(mp.dm1.act_ele==mp.dm1.tied(ti,1))==false);  mp.dm1.act_ele = [mp.dm1.act_ele; mp.dm1.tied(ti,1)];  end
                   if(any(mp.dm1.act_ele==mp.dm1.tied(ti,2))==false);  mp.dm1.act_ele = [mp.dm1.act_ele; mp.dm1.tied(ti,2)];  end
               end
               mp.dm1.act_ele = sort(mp.dm1.act_ele); %--Need to sort for the logic in model_Jacobian.m
           end
           if(any(mp.dm_ind==2))
               for ti=1:size(mp.dm2.tied,1)
                   if(any(mp.dm2.act_ele==mp.dm2.tied(ti,1))==false);  mp.dm2.act_ele = [mp.dm2.act_ele; mp.dm2.tied(ti,1)];  end
                   if(any(mp.dm2.act_ele==mp.dm2.tied(ti,2))==false);  mp.dm2.act_ele = [mp.dm2.act_ele; mp.dm2.tied(ti,2)];  end
               end
               mp.dm2.act_ele = sort(mp.dm2.act_ele);
           end
           if(any(mp.dm_ind==8))
               for ti=1:size(mp.dm8.tied,1)
                   if(any(mp.dm8.act_ele==mp.dm8.tied(ti,1))==false);  mp.dm8.act_ele = [mp.dm8.act_ele; mp.dm8.tied(ti,1)];  end
                   if(any(mp.dm8.act_ele==mp.dm8.tied(ti,2))==false);  mp.dm8.act_ele = [mp.dm8.act_ele; mp.dm8.tied(ti,2)];  end
               end
               mp.dm8.act_ele = sort(mp.dm8.act_ele);
           end
           if(any(mp.dm_ind==9))
               for ti=1:size(mp.dm9.tied,1)
                   if(any(mp.dm9.act_ele==mp.dm9.tied(ti,1))==false);  mp.dm9.act_ele = [mp.dm9.act_ele; mp.dm9.tied(ti,1)];  end
                   if(any(mp.dm9.act_ele==mp.dm9.tied(ti,2))==false);  mp.dm9.act_ele = [mp.dm9.act_ele; mp.dm9.tied(ti,2)];  end
               end
               mp.dm9.act_ele = sort(mp.dm9.act_ele);
           end
            
            %--Update the number of elements used per DM
            if(any(mp.dm_ind==1)); mp.dm1.Nele = length(mp.dm1.act_ele); end
            if(any(mp.dm_ind==2)); mp.dm2.Nele = length(mp.dm2.act_ele); end
            if(any(mp.dm_ind==8)); mp.dm8.Nele = length(mp.dm8.act_ele); end
            if(any(mp.dm_ind==9)); mp.dm9.Nele = length(mp.dm9.act_ele); end
            %mp.NelePerDMvec = [length(mp.dm1.Nele), length(mp.dm2.Nele), length(mp.dm3.Nele), length(mp.dm4.Nele), length(mp.dm5.Nele), length(mp.dm6.Nele), length(mp.dm7.Nele), length(mp.dm8.Nele), length(mp.dm9.Nele) ];

            if(any(mp.dm_ind==1)); fprintf('  DM1: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm1.Nele, mp.dm1.NactTotal,100*mp.dm1.Nele/mp.dm1.NactTotal); end
            if(any(mp.dm_ind==2)); fprintf('  DM2: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm2.Nele, mp.dm2.NactTotal,100*mp.dm2.Nele/mp.dm2.NactTotal); end
            if(any(mp.dm_ind==8)); fprintf('  DM8: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm8.Nele, mp.dm8.NactTotal,100*mp.dm8.Nele/mp.dm8.NactTotal); end
            if(any(mp.dm_ind==9)); fprintf('  DM9: %d/%d (%.2f%%) actuators kept for Jacobian\n', mp.dm9.Nele, mp.dm9.NactTotal,100*mp.dm9.Nele/mp.dm9.NactTotal); end
            
            %--Crop out unused actuators from the control Jacobian
            if(any(mp.dm_ind==1)); jacStruct.G1 = jacStruct.G1(:,mp.dm1.act_ele,:); end
            if(any(mp.dm_ind==2)); jacStruct.G2 = jacStruct.G2(:,mp.dm2.act_ele,:); end
            if(any(mp.dm_ind==8)); jacStruct.G8 = jacStruct.G8(:,mp.dm8.act_ele,:); end
            if(any(mp.dm_ind==9)); jacStruct.G9 = jacStruct.G9(:,mp.dm9.act_ele,:); end
        end  

end %--END OF FUNCTION falco_ctrl_cull

