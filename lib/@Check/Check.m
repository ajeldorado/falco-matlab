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
            % Check if a variable is a finite, real-valued scalar.
            if ~isnumeric(x)
                error('real_scalar:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('real_scalar:InputMustBeScalar', 'input is not a scalar')
            end

            if isinf(x)
                error('real_scalar:InputMustBeFinite', 'input is infinite')
            end

            if ~isreal(x)
                error('real_scalar:InputMustBeReal', 'input is not real-valued')
            end
        end
        
        
        function real_nonnegative_scalar(x)
            % Check if a variable is a nonnegative, finite, real-valued scalar.
            if ~isnumeric(x)
                error('real_nonnegative_scalar:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('real_nonnegative_scalar:InputMustBeScalar', 'input is not a scalar')
            end

            if isinf(x)
                error('real_nonnegative_scalar:InputMustBeFinite', 'input is infinite')
            end

            if ~isreal(x)
                error('real_nonnegative_scalar:InputMustBeReal', 'input is not real-valued')
            end
            
            if x < 0
                error('real_nonnegative_scalar:InputMustBeNonnegative', 'input is negative')
            end
        end
        
        
        function scalar_integer(x)
            % Check if a variable is a scalar integer.
            if ~isnumeric(x)
                error('scalar_integer:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('scalar_integer:InputMustBeScalar', 'input is not a scalar')
            end

            if isinf(x)
                error('scalar_integer:InputMustBeFinite', 'input is infinite')
            end

            if ~isreal(x)
                error('scalar_integer:InputMustBeReal', 'input is not real-valued')
            end

            if x ~= floor(x)
                error('scalar_integer:InputMustBeIntegral', 'input is not a integer')
            end
        end

        
        function nonnegative_scalar_integer(x)
            % Check if a variable is a nonnegative scalar integer.
            if ~isnumeric(x)
                error('nonnegative_scalar_integer:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('nonnegative_scalar_integer:InputMustBeScalar', 'input is not a scalar')
            end

            if isinf(x)
                error('nonnegative_scalar_integer:InputMustBeFinite', 'input is infinite')
            end

            if ~isreal(x)
                error('nonnegative_scalar_integer:InputMustBeReal', 'input is not real-valued')
            end

            if x ~= floor(x)
                error('nonnegative_scalar_integer:InputMustBeIntegral', 'input is not a integer')
            end
            
            if x < 0
                error('nonnegative_scalar_integer:InputMustBeNonnegative', 'input is negative')
            end
        end        
        
        
        function positive_scalar_integer(x)
            % Check if a variable is a positive scalar integer.
            if ~isnumeric(x)
                error('positive_scalar_integer:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) ~= 1
                error('positive_scalar_integer:InputMustBeScalar', 'input is not a scalar')
            end

            if isinf(x)
                error('positive_scalar_integer:InputMustBeFinite', 'input is infinite')
            end

            if ~isreal(x)
                error('positive_scalar_integer:InputMustBeReal', 'input is not real-valued')
            end

            if x ~= floor(x)
                error('positive_scalar_integer:InputMustBeIntegral', 'input is not a integer')
            end
            
            if x <= 0
                error('positive_scalar_integer:InputMustBePositive', 'input is not positive')
            end
        end

        
        function two_dim_array(x)
            % Check if a variable is a 2-D array.
            if ~isnumeric(x)
                error('two_dim_array:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) == 0
                error('two_dim_array:InputMustBeArray', 'input is empty')
            elseif numel(x) == 1
                error('two_dim_array:InputMustBeArray', 'input is a scalar')
            end

            if length(size(x)) ~= 2
                error('two_dim_array:InputMustBeTwoDim', 'input is not two dimensional')
            end            
        end
        
        
        function two_dim_square_array(x)
            % Check if a variable is a square-shaped, 2-D array.
            if ~isnumeric(x)
                error('two_dim_square_array:InputMustBeNumeric', 'input value is not numeric')
            end

            if numel(x) == 0
                error('two_dim_square_array:InputMustBeArray', 'input is empty')
            elseif numel(x) == 1
                error('two_dim_square_array:InputMustBeArray', 'input is a scalar')
            end

            if length(size(x)) ~= 2
                error('two_dim_square_array:InputMustBeTwoDim', 'input is not two dimensional')
            end
            
            if size(x, 1) ~= size(x, 2)
                error('two_dim_square_array:InputMustBeSquare', 'input is not a square array')
            end
        end
        
        
    end
   
    
end


