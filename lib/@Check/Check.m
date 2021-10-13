% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Checks on function inputs.
%
% These functions return nothing if input is valid. If input is invalid,
% then an error is thrown.

classdef Check
   
    methods (Static)

        
        function real_scalar(x)
            % Check if a variable is a real-valued scalar.
            if ~isnumeric(x)
                error('ValueError:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('ValueError:InputMustBeScalar', 'input is not a scalar')
            end

            if ~isreal(x)
                error('ValueError:InputMustBeReal', 'input is not real-valued')
            end
        end
        
        
        function real_nonnegative_scalar(x)
            % Check if a variable is a nonnegative, real-valued scalar.
            Check.real_scalar(x)
            if x < 0
                error('ValueError:InputMustBeNonnegative', 'input is negative')
            end
        end
        
        
        function real_positive_scalar(x)
            % Check if a variable is a positive, real-valued scalar.
            Check.real_scalar(x)
            if x <= 0
                error('ValueError:InputMustBePositive', 'input is negative')
            end
        end
        
        
        function scalar_integer(x)
            % Check if a variable is a scalar integer.
            if ~isnumeric(x)
                error('ValueError:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('ValueError:InputMustBeScalar', 'input is not a scalar')
            end

            if isinf(x)
                error('ValueError:InputMustBeFinite', 'input is infinite')
            end

            if ~isreal(x)
                error('ValueError:InputMustBeReal', 'input is not real-valued')
            end

            if x ~= floor(x)
                error('ValueError:InputMustBeIntegral', 'input is not a integer')
            end
        end

        
        function nonnegative_scalar_integer(x)
            % Check if a variable is a nonnegative scalar integer.
            Check.scalar_integer(x)
            if x < 0
                error('ValueError:InputMustBeNonnegative', 'input is negative')
            end
        end        
        
        
        function positive_scalar_integer(x)
            % Check if a variable is a positive scalar integer.
            Check.scalar_integer(x)
            if x <= 0
                error('ValueError:InputMustBePositive', 'input is not positive')
            end
        end

        
        function two_dim_array(x)
            % Check if a variable is a 2-D array.
            if ~isnumeric(x)
                error('ValueError:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) == 0
                error('ValueError:InputMustBeArray', 'input is empty')
            elseif numel(x) == 1
                error('ValueError:InputMustBeArray', 'input is a scalar')
            end

            if length(size(x)) ~= 2
                error('ValueError:InputMustBeTwoDim', 'input is not two dimensional')
            end            
        end
        
        
        function two_dim_square_array(x)
            % Check if a variable is a square-shaped, 2-D array.
            Check.two_dim_array(x)
            if size(x, 1) ~= size(x, 2)
                error('ValueError:InputMustBeSquare', 'input is not a square array')
            end
        end
        
        
    end
   
    
end


