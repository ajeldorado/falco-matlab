% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Check if a variable is a finite, real-valued scalar.
%
% INPUTS
% x: variable to check
%
% OUTPUTS
% -------
% None
%

function check_real_scalar(x)
    
    if ~isnumeric(x)
        error('check_real_scalar:InputMustBeNumeric', 'input value is not numeric')
    end
    
    if numel(x) ~= 1
        error('check_real_scalar:InputMustBeScalar', 'input is not a scalar')
    end
    
    if isinf(x)
        error('check_real_scalar:InputMustBeFinite', 'input is infinite')
    end
    
    if ~isreal(x)
        error('check_real_scalar:InputMustBeReal', 'input is not real-valued')
    end
    
end
