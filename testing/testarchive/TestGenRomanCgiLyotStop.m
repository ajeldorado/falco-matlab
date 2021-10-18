%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test falco_gen_Roman_CGI_lyot_stop_symm_fillet.m
%
% Test that falco_gen_Roman_CGI_lyot_stop_symm_fillet.m produces a
% similar mask to falco_gen_pupil_Roman_CGI_20200513.m.

classdef TestGenRomanCgiLyotStop < matlab.unittest.TestCase  
%% Properties
%
% A presaved file with FALCO parameters was saved and is loaded to be used
% by methods.
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
            addpath(genpath([pathToFalco filesep 'lib_external']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            pathToFalco = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
            rmpath(genpath([pathToFalco filesep 'lib']));
            addpath(genpath([pathToFalco filesep 'lib_external']));
        end
    end

%% Tests

    methods (Test)   
 
        function testObscurationCoverage(testCase)
            Nbeam = 250;
            ID = 0.50;
            OD = 0.80;
            wStrut = 3.6/100.;
            rocFillet = 0.001;
            upsampleFactor = 100;
            centering = 'pixel';
            lyotRounded = falco_gen_Roman_CGI_lyot_stop_symm_fillet(Nbeam, ID, OD, wStrut, rocFillet, upsampleFactor, centering);

            changes.flagLyot = true;
            changes.ID = ID;
            changes.OD = OD;
            changes.wStrut = wStrut;
            lyotSharp = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);

            nArray = ceil_odd(length(lyotRounded));
            lyotRounded = pad_crop(lyotRounded, nArray);
            lyotSharp = pad_crop(lyotSharp, nArray);

            % Symmetrize the sharp-cornered Lyot stop since the rounded one
            % is symmetric about the vertical axis.
            lyotRounded = floor(lyotRounded);
            lyotSharp = (lyotSharp + fliplr(lyotSharp))/2;
            lyotSharp = floor(lyotSharp);

            diff = lyotSharp - lyotRounded;
            summedDiff = sum(sum(abs(diff)));

            testCase.verifyEqual(summedDiff, 0) 
        end
        
        function testAreaRatio(testCase)
            Nbeam = 200;
            ID = 0.50;
            OD = 0.80;
            wStrut = 3.6/100.;
            rocFillet = 0.001;
            upsampleFactor = 100;
            centering = 'pixel';
            lyotRounded = falco_gen_Roman_CGI_lyot_stop_symm_fillet(Nbeam, ID, OD, wStrut, rocFillet, upsampleFactor, centering);

            changes.flagLyot = true;
            changes.ID = ID;
            changes.OD = OD;
            changes.wStrut = wStrut;
            lyotSharp = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);

            areaRatio = sum(sum(lyotRounded))/sum(sum(lyotSharp));
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(areaRatio, IsEqualTo(1, 'Within', RelativeTolerance(0.005)))           
        end
        
    end
end
