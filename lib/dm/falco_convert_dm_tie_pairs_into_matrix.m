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
% outList : Nx2 array specifying linear indices of N tied actuator pairs.
%
% OUTPUTS
% -------
% tieMat : 2-D array specifying tied actuators

function tieMat = falco_convert_dm_tie_pairs_into_matrix(tiePairs, nact)

    % Input checks
    if ~isempty(tiePairs)
        Check.two_dim_array(tiePairs)
        if size(tiePairs, 2) ~= 2
            error('The array of tied actuator pairs must have exactly two columns.')
        end
    end
    Check.positive_scalar_integer(nact)
   
    tieMat = zeros(nact, nact); % initialize
    
    if ~isempty(tiePairs)
        
        % Sort values within each row, and then sort rows.
        tiePairs = sortrows(sort(tiePairs, 2));
                
        tieBlockNum = 1; % initialize tie block ID number. 
        for iRow = 1:size(tiePairs, 1)
            
            isRepeatedValue = (tiePairs(iRow, 1) == tiePairs(iRow, 2));
            
            if ~isRepeatedValue
                
                maxTieMatVal = max(tieMat(tiePairs(iRow, :)));
                minTieMatVal = min(tieMat(tiePairs(iRow, :)));
                useAnExistingBlock = max(tieMat(tiePairs(iRow, :))) > 0;

                if useAnExistingBlock && (minTieMatVal > 0) 
                    tieMat(tieMat == maxTieMatVal) = minTieMatVal;
                    tieMat(tiePairs(iRow, :)) = minTieMatVal;
                
                elseif useAnExistingBlock
                    tieMat(tiePairs(iRow, :)) = maxTieMatVal;
                                                            
                else % assign a new block number
                    tieMat(tiePairs(iRow, :)) = tieBlockNum;
                    tieBlockNum = tieBlockNum + 1;
                end
                
            end
                
        end
        
    end
    
    % Verify that tieMat satisfies the format constraints
    isValid = ConstrainDM.checktie(tieMat);
    if ~isValid
        error('tied actuator matrix has invalid values.');
    end

end
