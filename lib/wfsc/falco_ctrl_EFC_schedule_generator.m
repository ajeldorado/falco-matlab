% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function for 
% -This function 
%
%
%
% REVISION HISTORY:
%--Heavily modified on 2018-07-24 by A.J. Riggs to change how the inputs
%are done.
%--Created on 2018-05-09 by A.J. Riggs.


function [Nitr, relinItrVec, gridSearchItrVec, log10regSched, dm_ind_sched] = falco_ctrl_EFC_schedule_generator(sched_mat)

    %--Number of correction iterations
    Nitr = sum(sched_mat(:,1));


    %--Create the vectors of 
    %  1) iteration numbers at which to relinearize the Jacobian
    %  2) log10(regularization) at each correction iteration
    relinItrVec = []; %--Initialize
    gridSearchItrVec = []; %--Initialize
    log10regSched = zeros(Nitr,1); %--Initialize
    dmIndList = zeros(Nitr,1); %--Initialize this temporary variable
    ItCounter = 1;
    for ii=1:size(sched_mat,1)

        if( (sched_mat(ii,4)==1) ) %--When to re-linearize
            relinItrVec = [relinItrVec, ItCounter];
        end
        
        if( (sched_mat(ii,5)==1) ) %--When to re-do the empirical EFC grid search
            gridSearchItrVec = [gridSearchItrVec, ItCounter];
        end
        
        
        if( (sched_mat(ii,1)~=0) )
            %--Make the vector of regularizations at each iteration
            log10regSched(ItCounter:(ItCounter+sched_mat(ii,1)-1)) = sched_mat(ii,2);
            dmIndList(ItCounter:(ItCounter+sched_mat(ii,1)-1)) = sched_mat(ii,3);
        end

        ItCounter = ItCounter + sched_mat(ii,1);
%         figure(201); plot(1:Nitr,log10regSched,'-bo');
%         figure(202); plot(1:Nitr,dmIndList,'-bo');
%         keyboard
    end
    %figure(503); plot(1:Nitr,regSched);

    %--Store DM number index vectors as cells since they can vary in length
    for Itr = 1:Nitr
        dm_ind = [];
        numAsStr = num2str(dmIndList(Itr));
        for ii=1:length(numAsStr)
            dm_ind = [dm_ind,str2num(numAsStr(ii))];
        end        
        dm_ind_sched{Itr} = dm_ind;
    end
    



end %--END OF FUNCTION
