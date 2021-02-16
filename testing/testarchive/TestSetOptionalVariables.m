%% Test User Defined Sinc Function
%
% We apply usual tests such as xx= 0, pi/6, pi/4/, pi/2, pi, 3pi/2, etc.
% The test assumes that the user sets path to FALCO functionality but there
% is an extra test which verifies that the sinc.m is actually a file in the
% ../lib/utils in the FALCO path.
classdef TestSetOptionalVariables < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        mp=DefaultParameters();
    end
     
%% Tests
%
%  Creates tests:
%
% # *testPaths:* verify falco_set_optional_variables.m 
%     creates mp fields for breifs, ws, maps, jac, images, DM, and
%     wsInProgress. Also verify that these fields are true paths in local
%     folder running FALCO's tests.
% # *testBooleanFlags:* verify that the output from
%     falco_set_optional_variables.m contains the boolean fields
%     flagSaveWS, flagSaveEachItr, flagSVD, flagUseLearnedJac, flagUseJac,
%     and flagUseModel. In addition we verify that the value of all these
%     fields is false.
%
    methods (Test)
        function testPaths(testCase)
            import matlab.unittest.constraints.IsFolder;
            import matlab.unittest.constraints.HasField

            mp = falco_set_optional_variables(testCase.mp);
            testCase.verifyThat(mp.path.falco,IsFolder)
            testCase.verifyThat(mp.path.proper,IsFolder)
            testCase.verifyThat(mp.path.config,IsFolder)
            testCase.verifyThat(mp.path.ws,IsFolder)
            testCase.verifyThat(mp.path.falcoaps,IsFolder)
            testCase.verifyThat(mp.path.jac,IsFolder)
            testCase.verifyThat(mp.path.images,IsFolder)
            testCase.verifyThat(mp.path.dm,IsFolder)
            %testCase.verifyThat(mp.path.wsInProgress,IsFolder)
            
            testCase.verifyThat(mp.path, HasField('falco'))
            testCase.verifyThat(mp.path, HasField('proper'))
            testCase.verifyThat(mp.path, HasField('config'))
            testCase.verifyThat(mp.path, HasField('ws'))
            testCase.verifyThat(mp.path, HasField('falcoaps'))
            testCase.verifyThat(mp.path, HasField('jac'))
            testCase.verifyThat(mp.path, HasField('images'))
            testCase.verifyThat(mp.path, HasField('dm'))
            %testCase.verifyThat(mp.path, HasField('wsInProgress'))
        end
        function testBooleanFlags(testCase)            
            import matlab.unittest.constraints.HasField
            mp = falco_set_optional_variables(testCase.mp);
            testCase.verifyThat(mp, HasField('flagSaveWS'))
            testCase.verifyThat(mp, HasField('flagSaveEachItr'))
            testCase.verifyThat(mp, HasField('flagSVD'))
            testCase.verifyThat(mp, HasField('flagUseLearnedJac'))
            testCase.verifyThat(mp.est, HasField('flagUseJac'))
            testCase.verifyThat(mp.ctrl, HasField('flagUseModel'))
            
            testCase.verifyFalse(mp.flagSaveWS,false)
            testCase.verifyFalse(mp.flagSaveEachItr,false)
            testCase.verifyFalse(mp.flagSVD,false)
            testCase.verifyFalse(mp.flagUseLearnedJac,false)
            testCase.verifyFalse(mp.est.flagUseJac,false)
            testCase.verifyFalse(mp.ctrl.flagUseModel,false)
        end        
        
    end    
end