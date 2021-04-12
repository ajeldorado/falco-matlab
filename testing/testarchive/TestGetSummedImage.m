classdef TestGetSummedImage < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is load to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        mp=ImageParameters();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'lib']))
        end
    end   
    

%% Tests
%
%  Creates tests:
%
% # *test Image Size* with actal size equal to [mp.Fend.Nxi,mp.Fend.Nxi]
% 
    methods (Test)
        
        function testImageSize(testCase)
            import matlab.unittest.constraints.HasSize
            act = falco_get_summed_image(testCase.mp);
            testCase.assertThat(act, HasSize([testCase.mp.Fend.Nxi testCase.mp.Fend.Nxi]))
        end
         
        function  testParallelVsSerialImageResult(testCase)
         
            disp('Generating image without parfor.')
            testCase.mp.flagParfor = false;
            imageSerial = falco_get_summed_image(testCase.mp);
            
            disp('Generating image with parfor.')
            testCase.mp.flagParfor = true;
            imageParallel = falco_get_summed_image(testCase.mp);
            
            absDiff = abs(imageSerial - imageParallel);
            
            testCase.verifyLessThan(absDiff, 1e-12)
            
            
        end
    end    
end