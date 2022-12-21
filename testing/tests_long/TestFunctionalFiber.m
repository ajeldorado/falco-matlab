%% Functional Test Vortex
%
% The test script will perform the Wavefront Sensing and Control first, then it will
% test verify the values of the outputs to expected values define for each
% test.
classdef TestFunctionalFiber < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is loaded to be used
% by methods. In this case we use the mp.path.falco to addpath to the
% function being tested.
    properties
        mp=ConfigurationVortex();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco filesep 'models']));
            addpath(genpath([testCase.mp.path.falco filesep 'setup']));
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'models']))
            rmpath(genpath([testCase.mp.path.falco filesep 'setup']))
            rmpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end  
    
%% *Tests*

    methods (Test)     
        function testFunctionalVortex(testCase)
            mp = testCase.mp;
            
            %% Step 3: Perform the Wavefront Sensing and Control
            mp.runLabel = 'test_fiber';
            
            % % Fiber parameters
            mp.flagFiber = true;
            mp.flagLenslet = false;

            mp.Fend.x_fiber = [6,5];%[5.3405 -2.6702 -2.6702]; %Fiber core center positions in lambda_0/D
            mp.Fend.y_fiber = [0,-3];%[0 4.625 -4.625];
            mp.Fend.Nfiber = numel(mp.Fend.x_fiber);

            mp.fiber.a = 0.507;%0.875;%0.66; %Radius of the fiber core in lambda_0/D
            mp.fiber.a_phys = 1.75e-6; %Physical radius of the fiber core in meters
            mp.fiber.NA = 0.12; %Numerical aperture of the fiber
            
            %--Use just 1 wavelength
            mp.fracBW = 0.1;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
            mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
            mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

            mp.Nitr = 1; %--Number of wavefront control iterations
            mp.estimator = 'pairwise';

            [mp, out] = falco_flesh_out_workspace(mp);
            [mp, out] = falco_wfsc_loop(mp, out);

            %% Tests:
            Iend = out.InormFiberHist(end, 1); % 1.9032e-13
            testCase.verifyGreaterThan(Iend, 1e-13)
            testCase.verifyLessThan(Iend, 3e-13)
            Iend = out.InormFiberHist(end, 2); % 2.9146e-9
            testCase.verifyGreaterThan(Iend, 2.0e-9)
            testCase.verifyLessThan(Iend, 4.0e-9)
            
            dm1pv = out.dm1.Spv(end); % 1.5057e-08
            testCase.verifyGreaterThan(dm1pv, 3.0e-11)
            testCase.verifyLessThan(dm1pv, 4.0e-11)
            
            thput = out.thput(end, 1); % [0.3198, 0.3143]
            testCase.verifyGreaterThan(thput, 0.315)
            testCase.verifyLessThan(thput, 0.325)
            thput = out.thput(end, 2); % [0.3198, 0.3143]
            testCase.verifyGreaterThan(thput, 0.310)
            testCase.verifyLessThan(thput, 0.320)
            
            import matlab.unittest.constraints.EveryElementOf
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(out.log10regHist, IsEqualTo([-2.5])) 
        end       
    end    
end





