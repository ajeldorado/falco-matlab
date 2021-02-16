%% Test falco_set_jacobian_weights.m
%
% We define some tests for falco_set_jacobian_weights.m to test  
% error exceptions, and responses to a different input parameters. The 
% test assumes that the user has set the path to FALCO functionality. 
classdef TestJacWeights < matlab.unittest.TestCase
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
    
%% Tests
%
%  Creates four tests:
%
% # *testJacWeightserror1* The parameter mp.jac.zerns repeats a zernike in 
%         position two as follows zernike=[1,1,3]. When this happens a call 
%         to the function generates an error. This test verifies that when
%         we do this there is in fact an error exception thrown.
% # *testJacWeightserror2* When the parameter mp.jac.zerns is missing piston
%         a call to the function generates an error exception. This test 
%         verifies this is true.
% # *testJacWeightsinputs1* Inputs1 test case where mp.estimator = 'perfect',
%         mp.jac.zerns = 1, mp.jac.Zcoef = 1e-9, mp.Nsbp = 1,
%         mp.si_ref = ceil(mp.Nsbp/2), and mp = falco_set_jacobian_weights(mp).
%         These inputs generate the following results:  mp.jac.Nzern = 1, 
%         mp.jac.zerns = 1, mp.jac.Zcoef = Zcoef: 1, mp.jac.weights = 1,
%         mp.jac.Nmode = 1, mp.jac.weightMat = 1, mp.jac.sbp_inds = 1,
%         mp.jac.zern_inds = 1. This test verifies this is the case.
% # *testJacWeightsinputs2* Inputs2 test case where mp.estimator =
%         'perfect', mp.jac.zerns = [1, 5, 6], mp.jac.Zcoef = [1, 1,1]*1e-9,
%         mp.Nsbp = 3, and mp.si_ref = ceil(mp.Nsbp/2). The results should
%         be: mp.jac.Nzern = 3, mp.jac.zerns = [1, 5, 6], 
%         mp.jac.Zcoef = Zcoef: [1 1.0000e-09 1.0000e-09],
%         mp.jac.weights = repmat([.25; 0.5; 0.25], [3,1]),
%         mp.jac.Nmode = 9, mp.jac.weightMat = [0.25*ones(1,3), 
%         0.5*ones(1,3); 0.25*ones(1,3)], mp.jac.sbp_inds = [1, 2, 3, 1,
%         2, 3, 1, 2, 3].', and mp.jac.zern_inds = [1, 1, 1, 5, 5, 5, 6,6,6].'.
% # *testJacWeightsinputs3* Inputs test case where mp.estimator = 'pwp-bp', 
%         mp.jac.zerns = [1, 5, 6], mp.jac.Zcoef = [1, 1, 1]*1e-9, 
%         mp.Nsbp= 3, mp.si_ref = ceil(mp.Nsbp/2). The results should be: 
%         mp.jac.Nzern = 3, mp.jac.zerns = [1, 5, 6], mp.jac.Zcoef = [1 
%         1.0000e-09 1.0000e-09], mp.jac.weights = (1/3)*ones(9,1), 
%         mp.jac.Nmode = 9, mp.jac.weightMat = (1/3)*ones(3,3),
%         mp.jac.sbp_inds = [1, 2, 3, 1, 2, 3, 1, 2, 3].',
%         mp.jac.zern_inds = [1, 1, 1, 5, 5, 5, 6, 6, 6].'. This test verifies 
%         that this is the case.  
    methods (Test)
        function testJacWeightserror1(testCase)
            import matlab.unittest.constraints.Throws
            mp=testCase.mp;
            mp.estimator = 'perfect';
            mp.jac.zerns = [1, 1, 3];
            mp.jac.Zcoef = [1, 1, 1]*1e-9;
            mp.Nsbp = 2;
            mp.si_ref = ceil(mp.Nsbp/2);
            testCase.verifyThat(@() falco_set_jacobian_weights(mp),Throws(?MException))
        end
        function testJacWeightserror2(testCase)
            import matlab.unittest.constraints.Throws
            mp=testCase.mp;
            mp.estimator = 'perfect';
            mp.jac.zerns = [1, 2, 3, 4];
            mp.jac.Zcoef = [1, 1, 1]*1e-9;
            mp.Nsbp = 2;
            mp.si_ref = ceil(mp.Nsbp/2);
            testCase.verifyThat(@() falco_set_jacobian_weights(mp),Throws(?MException))
        end       
        function testJacWeightsinputs1(testCase)
            mp=testCase.mp;
            mp.estimator = 'perfect';
            mp.jac.zerns = 1;
            mp.jac.Zcoef = 1e-9;
            mp.Nsbp = 1;
            mp.si_ref = ceil(mp.Nsbp/2);
            mp = falco_set_jacobian_weights(mp);

            testCase.verifyEqual(1,mp.jac.Nzern);  
            testCase.verifyEqual(1,mp.jac.zerns); 
            testCase.verifyEqual(1,mp.jac.Zcoef); 
            testCase.verifyEqual(1,mp.jac.weights); 
            testCase.verifyEqual(1,mp.jac.Nmode);
            testCase.verifyEqual(1,mp.jac.weightMat);
            testCase.verifyEqual(1,mp.jac.sbp_inds);
            testCase.verifyEqual(1,mp.jac.zern_inds);
        end
        function testJacWeightsinputs2(testCase)
            mp=testCase.mp;
            mp.estimator = 'perfect';
            mp.jac.zerns = [1, 5, 6];
            mp.jac.Zcoef = [1, 1, 1]*1e-9;
            mp.Nsbp = 3;
            mp.si_ref = ceil(mp.Nsbp/2);
            mp = falco_set_jacobian_weights(mp);
            % Expected Outputs:    
            actSolution1= 3;
            actSolution2= [1, 5, 6];
            actSolution3= [1 1.0000e-09 1.0000e-09];
            actSolution4= repmat([.25; 0.5; 0.25], [3,1]);
            actSolution5= 9;
            actSolution6=[0.25*ones(1,3); 0.5*ones(1,3); 0.25*ones(1,3)];
            actSolution7= [1, 2, 3, 1, 2, 3, 1, 2, 3].';
            actSolution8= [1, 1, 1, 5, 5, 5, 6, 6, 6].';

            testCase.verifyEqual(actSolution1,mp.jac.Nzern); 
            testCase.verifyEqual(actSolution2,mp.jac.zerns); 
            testCase.verifyEqual(actSolution3,mp.jac.Zcoef); 
            testCase.verifyEqual(actSolution4,mp.jac.weights); 
            testCase.verifyEqual(actSolution5,mp.jac.Nmode); 
            testCase.verifyEqual(actSolution6,mp.jac.weightMat);
            testCase.verifyEqual(actSolution7,mp.jac.sbp_inds);
            testCase.verifyEqual(actSolution8,mp.jac.zern_inds);
        end
        function testJacWeightsinputs3(testCase)
            mp=testCase.mp;
            mp.estimator = 'pwp-bp' ;
            mp.jac.zerns = [1, 5, 6];
            mp.jac.Zcoef = [1, 1, 1]*1e-9;
            mp.Nsbp = 3;
            mp.si_ref = ceil(mp.Nsbp/2);
            mp = falco_set_jacobian_weights(mp);
            % Expected Outputs:    
            actSolution1= 3;
            actSolution2= [1, 5, 6];
            actSolution3= [1 1.0000e-09 1.0000e-09];
            actSolution4= (1/3)*ones(9,1);
            actSolution5= 9;
            actSolution6= (1/3)*ones(3,3);
            actSolution7= [1, 2, 3, 1, 2, 3, 1, 2, 3].';
            actSolution8= [1, 1, 1, 5, 5, 5, 6, 6, 6].';

            testCase.verifyEqual(actSolution1,mp.jac.Nzern); 
            testCase.verifyEqual(actSolution2,mp.jac.zerns); 
            testCase.verifyEqual(actSolution3,mp.jac.Zcoef); 
            testCase.verifyEqual(actSolution4,mp.jac.weights); 
            testCase.verifyEqual(actSolution5,mp.jac.Nmode); 
            testCase.verifyEqual(actSolution6,mp.jac.weightMat);
            testCase.verifyEqual(actSolution7,mp.jac.sbp_inds);
            testCase.verifyEqual(actSolution8,mp.jac.zern_inds);
         end
    end    
end