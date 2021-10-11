% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Check if a variable is a 2-D square array.
%
% INPUTS
% x: variable to check
%
% OUTPUTS
% -------
% None
%

function check_two_dim_square_array(x)
    
    if ~isnumeric(x)
        error('check_two_dim_square_array:InputMustBeNumeric', 'input value is not numeric')
    end
    
    if numel(x) == 0
        error('check_two_dim_square_array:InputMustBeArray', 'input is empty')
    elseif numel(x) == 1
        error('check_two_dim_square_array:InputMustBeArray', 'input is a scalar')
    end
    
    if length(size(x)) ~= 2
        error('check_two_dim_square_array:InputMustBeTwoDim', 'input is not two dimensional')
    end
    
    if size(x, 1) ~= size(x, 2)
        error('check_two_dim_square_array:InputMustBeSquare', 'input is not a square array')
    end
    
end
