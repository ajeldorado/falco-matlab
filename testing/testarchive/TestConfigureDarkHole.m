%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test falco_gen_SW_mask.m
%
% We define some tests for falco_configure_dark_hole_region.m to test responses to
% different input parameters. 
classdef TestConfigureDarkHole < matlab.unittest.TestCase    
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
%     properties
%         mp=Parameters();
%     end

%% Setup and Teardown Methods
%
%  Add and remove path to library functions to be tested.

    methods (TestClassSetup)
        function addPath(testCase)
            pathToFalco = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
            addpath(genpath([pathToFalco filesep 'lib']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            pathToFalco = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
            rmpath(genpath([pathToFalco filesep 'lib']));
        end
    end
    
%% Tests
                               
    methods (Test) 
        
        function testCircularSoftwareMaskArea(testCase)
            
            %--Correction and scoring region definition
            mp.Fend.corr.Rin = 2;   % inner radius of dark hole correction region [lambda0/D]
            mp.Fend.corr.Rout  = 10;  % outer radius of dark hole correction region [lambda0/D]
            mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]
            mp.Fend.score = mp.Fend.corr;
            
            mp.centering = 'pixel';
            mp.Fend.sides = 'lr';
            mp.Fend.shape = 'circle';
            mp.Fend.res = 10;
            % mp.Fend.xiOffset = 6;
            
            mp.Fend.eval.res = 20;
            mp.flagFiber = false;
            mp.thput_eval_x = 7;
            mp.thput_eval_y = 0;
            
            mp = falco_configure_dark_hole_region(mp);
            
            area = sum(sum(mp.Fend.corr.maskBool));
            
            areaExpected = pi*(mp.Fend.corr.Rout^2 - mp.Fend.corr.Rin^2)*(2*mp.Fend.corr.ang/360)*(mp.Fend.res^2);
            
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(area, IsEqualTo(areaExpected,'Within', RelativeTolerance(0.001)))
        end

        function testTwoZoneMask(testCase)
            
            %--Correction and scoring region definition
            mp.Fend.corr.Rin = [2, 2];   % inner radius of dark hole correction region [lambda0/D]
            mp.Fend.corr.Rout  = [5, 5];  % outer radius of dark hole correction region [lambda0/D]
            mp.Fend.corr.ang  = [150, 180];  % angular opening of dark hole correction region [degrees]
            mp.Fend.score = mp.Fend.corr;
            
            mp.centering = 'pixel';
            mp.Fend.sides = {'lr', 'lr'};
            mp.Fend.shape = {'circle', 'square'};
            mp.Fend.res = 10;
            mp.Fend.xiOffset = [0, 20];
%             mp.Fend.xiFOV = 30;
%             mp.Fend.etaFOV = 30;
            
            mp.Fend.eval.res = 20;
            mp.flagFiber = false;
            mp.thput_eval_x = 7;
            mp.thput_eval_y = 0;
            
            mp = falco_configure_dark_hole_region(mp);
            
            area = sum(sum(mp.Fend.corr.maskBool));
            
            areaExpected = pi*(mp.Fend.corr.Rout(1)^2 - mp.Fend.corr.Rin(1)^2)*(2*mp.Fend.corr.ang(1)/360)*(mp.Fend.res^2) + ...
                        (4*mp.Fend.corr.Rout(2)^2 - pi*mp.Fend.corr.Rin(2)^2)*(mp.Fend.res^2);
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(area, IsEqualTo(areaExpected,'Within', RelativeTolerance(0.02)))
        end        
        
    end    
end