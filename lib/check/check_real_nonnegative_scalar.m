% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Check if a variable is a nonnegative, finite, real-valued scalar.
%
% INPUTS
% x: variable to check
%
% OUTPUTS
% -------
% None
%

function check_real_nonnegative_scalar(x)
    
    if ~isnumeric(x)
        error('check_real_nonnegative_scalar:InputMustBeNumeric', 'input value is not numeric')
    end
    
    if numel(x) ~= 1
        error('check_real_nonnegative_scalar:InputMustBeScalar', 'input is not a scalar')
    end
    
    if isinf(x)
        error('check_real_nonnegative_scalar:InputMustBeFinite', 'input is infinite')
    end
    
    if ~isreal(x)
        error('check_real_nonnegative_scalar:InputMustBeReal', 'input is not real-valued')
    end
    
    if x < 0
        error('check_real_nonnegative_scalar:InputMustBeNonnegative', 'input is negative')
    end
    
end
