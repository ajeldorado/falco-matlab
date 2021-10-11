% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Check if a variable is an real-valued scalar integer.
%
% INPUTS
% x: variable to check
%
% OUTPUTS
% -------
% None
%

function check_scalar_integer(x)
    
    if ~isnumeric(x)
        error('check_scalar_integer:InputMustBeNumeric', 'input value is not numeric')
    end
    
    if numel(x) ~= 1
        error('check_scalar_integer:InputMustBeScalar', 'input is not a scalar')
    end
    
    if isinf(x)
        error('check_scalar_integer:InputMustBeFinite', 'input is infinite')
    end
    
    if ~isreal(x)
        error('check_scalar_integer:InputMustBeReal', 'input is not real-valued')
    end
    
    if x ~= floor(x)
        error('check_scalar_integer:InputMustBeIntegral', 'input is not a integer')
    end
    
end
