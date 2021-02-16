%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test DM Apodized V Simple Design
%
% We test pupil Luvoir B (translation and rotation).
classdef TestPupilRomanCGI < matlab.unittest.TestCase   
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
% # *testPupilRomanCGITranslation* Verify the pupil Luvoir B translation
%                            meet given constraints.
% # *testPupilRomanCGITranslationandRotation*  Verify the pupil Luvoir B translation 
%                                        and rotation meet given constraints.
% # *testPupilRomanCGIfromFile*  Verify the pupil Luvoir B translation 
%                                        and rotation meet given constraints.
%
    methods (Test)    
        function testPupilRomanCGITranslation(testCase)
            Nbeam = 100;
            centering = 'pixel';
            pupil = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering);
            
            % Translation test
            changes.xShear = -11/100;
            changes.yShear = 19/100;
            pupilOffset = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);
            diff = pad_crop(pupil, size(pupilOffset)) - circshift(pupilOffset, -Nbeam*[changes.yShear, changes.xShear]);
            testCase.verifyLessThan(sum(abs(diff(:))),1e-8)
        end
        function testPupilRomanCGITranslationandRotation(testCase)
            Nbeam = 100;
            centering = 'pixel';
            pupil = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering);
            
            % Translation test
            changes.xShear = -11/100;
            changes.yShear = 19/100;
            pupilOffset = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);
            
            
            % Test rotation (and translation)
            changes.clock_deg = 90;
            pupilRotOffset = falco_gen_pupil_Roman_CGI_20200513(Nbeam, centering, changes);
            % figure(5); imagesc(pupilRotOffset); axis xy equal tight; colorbar; drawnow;
            pupilRot = zeros(size(pupil));
            pupilRot(2:end, 2:end) = rot90(pupil(2:end, 2:end),-1);
            diff = pad_crop(pupilRot, size(pupilOffset)) - circshift(pupilRotOffset, -Nbeam*[changes.yShear, changes.xShear]);
            
            testCase.verifyLessThan(sum(abs(diff(:))/sum(pupilRotOffset(:))),1e-4)
        end
        function testPupilRomanCGIfromFile(testCase)
            pupil0 = imread('pupil_CGI-20200513_8k_binary_noF.png');
            pupilFromFile = double(pupil0)/double(max(pupil0(:)));
            Narray = size(pupilFromFile,1);
            dsfac = 16; % downsampling factor
            Nbeam = 2*4027.25;
            NbeamDS = Nbeam/dsfac;
            NarrayDS = Narray/dsfac;
            pupilFromFileDS = falco_bin_downsample(pupilFromFile, dsfac);
            
            % Generate pupil representation in FALCO
            shift = dsfac/2-0.5;
            changes.xShear = -0.5/Nbeam - shift/Nbeam;
            changes.yShear = -52.85/Nbeam - shift/Nbeam;
            pupilFromFALCODS = falco_gen_pupil_Roman_CGI_20200513(NbeamDS, 'pixel', changes);
            pupilFromFALCODS = pad_crop(pupilFromFALCODS, NarrayDS);
            diff = pupilFromFileDS - pupilFromFALCODS;
            percentvalue=sum(sum(abs(diff)))/sum(sum(pupilFromFileDS))*100;
            testCase.verifyLessThan(percentvalue,0.1)
        end
    end    
end