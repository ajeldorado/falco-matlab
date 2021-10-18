% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Convert the 2-D array specifying tied actuators into a list of actuator
% pairs. The tied actuator matrix is used by the ConstrainDM class.
%
% INPUTS
% ------
% tieMat : 2-D array specifying tied actuators
%
% OUTPUTS
% -------
% outList : Nx2 array specifying linear indices of N tied actuator pairs.

function tiePairs = falco_convert_dm_tie_matrix_into_pairs(tieMat)

    % Check whether a tie matrix satisfies the format constraints
    isValid = ConstrainDM.checktie(tieMat);
    if ~isValid
        error('tied actuator matrix has invalid values.');
    end
    
    tiePairs = zeros(0, 2);
    for tieBlockNum = 1:max(max(tieMat)) % ignore -1 (dead) and 0 (untied) values
        
        tempPairs = zeros(0, 2);
        inds = find(tieMat == tieBlockNum);
        nInds = length(inds);
        if nInds > 1
            tempPairs = inds(1)*ones(nInds-1, 2);
            tempPairs(:, 2) = inds(2:end);
        end
        
        tiePairs = [tiePairs; tempPairs];
    end
    
end
