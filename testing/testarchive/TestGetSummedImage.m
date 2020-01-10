classdef TestGetSummedImage < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        param=load('Parameters.mat');
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods(TestMethodSetup)
        function addPath(testCase)
            addpath(genpath([testCase.param.mp.path.falco 'lib/utils']));
        end
    end
    methods(TestMethodTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.param.mp.path.falco 'lib/utils']))
        end
    end
    
%% Tests
%
%  Creates four testrs:
%
% # *test Image Size* with actal size equal to [mp.Fend.Nxi,mp.Fend.Nxi]
% 
    methods (Test)
        function testImageSize(testCase)
            import matlab.unittest.constraints.HasSize
            act = falco_get_summed_image(testCase.param.mp);
            testCase.assertThat(act,HasSize([testCase.param.mp.Fend.Nxi testCase.param.mp.Fend.Nxi]))
        end
         
    end    
end