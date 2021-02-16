%% Test falco_set_jacobian_weights.m
%
% We define some tests for falco_set_spectral_properties.m to test  
% responses to different input parameters. 
classdef TestSetSpectralProperties < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        mp=Parameters();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco 'lib/utils']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco 'lib/utils']))
        end
    end    
    
%% *Properties*
%
% We defined a function called Parameters.m which when called it generates
% FALCO default parameters as well as some required parameters to run the
% main example in ../falco-matlab/main/EXAMPLE_main_WFIRST_LC.m.
    properties
        mp=Parameters();
    end
    
%% *Tests*
%
% *Creates five tests:*
%
% # *testSetSpectralPropertiesInput1* Input parameters for this test are 
%         mp.lambda0 = 1e-6, mp.flagSim = false, mp.Nsbp = 1, mp.Nwpsbp = 1, 
%         and mp.fracBW = 0.20. The results when calling the function with
%         these inputs should be mp.sbp_weights = 1, mp.sbp_centers = 1e-6, 
%         mp.full.lambda_weights_all = 1, and mp.full.lambdas = 1e-6. This 
%         test verifies that this are the 
%         actual results.    
% # *testSetSpectralPropertiesInput2* Input parameters for this test are 
%         mp.lambda0 = 1e-6, mp.flagSim = false, mp.Nsbp = 1, mp.Nwpsbp = 3, 
%         and mp.fracBW = 0.20. The results when calling the function with
%         these inputs should be mp.sbp_weights = 1, 
%         mp.sbp_centers = 1e-6, mp.full.lambda_weights_all = [0.25; 0.5; 0.25], 
%         mp.full.lambdas = linspace(0.9, 1.1, 3).'*1e-6. This test verifies 
%         that this are the actual results.
% # *testSetSpectralPropertiesInput3* Input parameters for this test are 
%         mp.lambda0 = 1e-6, mp.flagSim = true, mp.Nsbp = 5, mp.Nwpsbp = 1,
%         mp.fracBW = 0.20. The results when calling the function with
%         these inputs should be 
%         mp.sbp_weights = [0.1250; 0.25; 0.25; 0.25; 0.125],
%         mp.sbp_centers = linspace(0.9, 1.1, 5)*1e-6,
%         mp.full.lambda_weights_all = [0.1250; 0.25; 0.25; 0.25; 0.125],
%         and mp.full.lambdas = linspace(0.9, 1.1, 5)*1e-6. This test 
%         verifies that this are the actual results.  
% # *testSetSpectralPropertiesInput4* Input parameters for this test are 
%         mp.lambda0 = 1e-6, mp.flagSim = false, mp.Nsbp = 5, 
%         mp.Nwpsbp = 1, mp.fracBW = 0.20. The results when calling the 
%         function with these inputs should be 
%         mp.sbp_weights = 0.2*ones(5,1), mp.sbp_centers = (9.2:0.4:10.8)*1e-6
%         mp.full.lambda_weights_all = [0.05; 0.1*ones(9,1) 0.05],
%         mp.full.lambdas = linspace(0.9, 1.1, 11).'*1e-6. This test 
%         verifies that this are the actual results.
% # *testSetSpectralPropertiesInput5* Input parameters for this test are 
%         mp.lambda0 = 1e-6, mp.flagSim = false, mp.Nsbp = 5, mp.Nwpsbp = 3,
%         mp.fracBW = 0.20. The results when calling the function with
%         these inputs should be 
%         mp.sbp_weights = 0.2*ones(5,1), mp.sbp_centers = (9.2:0.4:10.8)*1e-7,
%         mp.full.lambda_weights_all = [0.05; 0.1*ones(9,1); 0.05], and
%         mp.full.lambdas = linspace(0.9, 1.1, 11).'*1e-6. This test 
%         verifies that this are the actual results.
%
    methods (Test)     
        function testSetSpectralPropertiesInput1(testCase)
            mp=testCase.mp;
            mp.lambda0 = 1e-6;
            mp.flagSim = false; 
            mp.Nsbp = 1;
            mp.Nwpsbp = 1;
            mp.fracBW = 0.20;
            mp = falco_set_spectral_properties(mp);
            testCase.verifyEqual(1,mp.sbp_weights);
            testCase.verifyEqual(1e-6,mp.sbp_centers); 
            testCase.verifyEqual(1,mp.full.lambda_weights_all); 
            testCase.verifyEqual(1e-6,mp.full.lambdas); 
        end
        function testSetSpectralPropertiesInput2(testCase)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            mp=testCase.mp;
            mp.lambda0 = 1e-6;
            mp.flagSim = false; 
            mp.Nsbp = 1;
            mp.Nwpsbp = 3;
            mp.fracBW = 0.20;
            mp = falco_set_spectral_properties(mp);
            testCase.verifyEqual(1,mp.sbp_weights);
            testCase.verifyEqual(1e-6,mp.sbp_centers); 
            testCase.verifyEqual([0.25; 0.5; 0.25],mp.full.lambda_weights_all); 
            testCase.verifyThat(linspace(0.9, 1.1, 3).'*1e-6,IsEqualTo(mp.full.lambdas,'Within',RelativeTolerance(0.001)))
        end             
        function testSetSpectralPropertiesInput3(testCase)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            mp=testCase.mp;
            mp.lambda0 = 1e-6;
            mp.flagSim = true;
            mp.Nsbp = 5;
            mp.Nwpsbp = 1;
            mp.fracBW = 0.20;
            mp = falco_set_spectral_properties(mp);
            testCase.verifyEqual([0.1250; 0.25; 0.25; 0.25; 0.125],mp.sbp_weights);
            testCase.verifyThat(linspace(0.9, 1.1, 5)*1e-6,IsEqualTo(mp.sbp_centers,'Within',RelativeTolerance(0.001)))
            testCase.verifyEqual([0.1250; 0.25; 0.25; 0.25; 0.125],mp.full.lambda_weights_all); 
            testCase.verifyThat(linspace(0.9, 1.1, 5).'*1e-6,IsEqualTo(mp.full.lambdas,'Within',RelativeTolerance(0.001)))
        end   
        function testSetSpectralPropertiesInput4(testCase)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            mp=testCase.mp;
            mp.lambda0 = 1e-6;
            mp.flagSim = false;
            mp.Nsbp = 5;
            mp.Nwpsbp = 1;
            mp.fracBW = 0.20;
            mp = falco_set_spectral_properties(mp);
            testCase.verifyEqual(0.2*ones(5,1),mp.sbp_weights);
            testCase.verifyThat(linspace(0.92, 1.08, 5)*1e-6,IsEqualTo(mp.sbp_centers,'Within',RelativeTolerance(0.001)))
            testCase.verifyEqual(0.2*ones(5,1),mp.full.lambda_weights_all); 
            testCase.verifyThat(linspace(0.92, 1.08, 5).'*1e-6,IsEqualTo(mp.full.lambdas,'Within',RelativeTolerance(0.001)))
        end 
        function testSetSpectralPropertiesInput5(testCase)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            mp=testCase.mp;
            mp.lambda0 = 1e-6;
            mp.flagSim = false;
            mp.Nsbp = 5;
            mp.Nwpsbp = 3;
            mp.fracBW = 0.20;
            mp = falco_set_spectral_properties(mp);
            testCase.verifyEqual(0.2*ones(5,1),mp.sbp_weights);
            testCase.verifyThat((9.2:0.4:10.8)*1e-7,IsEqualTo(mp.sbp_centers,'Within',RelativeTolerance(0.001)))
            testCase.verifyEqual([0.05; 0.1*ones(9,1); 0.05],mp.full.lambda_weights_all); 
            testCase.verifyThat(linspace(0.9, 1.1, 11).'*1e-6,IsEqualTo(mp.full.lambdas,'Within',RelativeTolerance(0.001)))
        end         
    end    
end





