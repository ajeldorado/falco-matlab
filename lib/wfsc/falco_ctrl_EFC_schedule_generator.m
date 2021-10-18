% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Decode the controller scheduling matrix into its components.
%
% CONTROL SCHEDULE. Columns of sched_mat are: 
% Column 1: # of iterations, 
% Column 2: log10(regularization), 
% Column 3: which DMs to use (12, 128, 129, or 1289) for control
% Column 4: flag (0 = false, 1 = true), whether to re-linearize
%   at that iteration.
% Column 5: flag (0 = false, 1 = true), whether to perform an
%   EFC parameter grid search to find the set giving the best
%   contrast .
% The imaginary part of the log10(regularization) in column 2 is
%  replaced for that iteration with the optimal log10(regularization)
% A row starting with [0, 0, 0, 1...] is for relinearizing only at that time


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

    end

    %--Store DM number index vectors as cells since they can vary in length
    for Itr = 1:Nitr
        dm_ind = [];
        numAsStr = num2str(dmIndList(Itr));
        for ii=1:length(numAsStr)
            dm_ind = [dm_ind,str2num(numAsStr(ii))];
        end        
        dm_ind_sched{Itr} = dm_ind;
    end

end
