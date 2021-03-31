%% Functional Test ZWFS
%
% The test script will first run the Zernike wavefront sensor, and then it will
% verify that the RMS sensing error is below some level.

classdef TestFunctionalZWFS < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is lodaded to be used
% by methods. In this case we use the mp.path.falco to addpath to the
% function being tested.
    properties
        mp=ConfigurationLC();
    end

%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath([testCase.mp.path.falco filesep 'setup']));
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'setup']))
            rmpath(genpath([testCase.mp.path.falco filesep 'lib']));
        end
    end    

    
%% *Tests*
%
% *Creates tests:*
%
% # *testFunctionalZWFS* Input parameters are predefined in the
% ConfigurationLC.m function which is called by the properties of the
% class and passed in to the test methods.
%
    methods (Test)     
        function testFunctionalZWFS(testCase)
            mp=testCase.mp;

            mp.runLabel = 'testing_ZWFS';
            mp.flagWFS = true; % Activate the WFS mode 
            mp.wfs.flagSim = true; % Simulates WFS images, if true 

            %% Pupil definitions

            mp.flagApod = false;
            mp.whichPupil = 'LUVOIR_B_offaxis';

            mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.

            mp.P1.full.Nbeam = 500; %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
            mp.P4.full.Nbeam = mp.P1.full.Nbeam; % P4 must be the same as P1 for Vortex. 

            mp.P1.compact.Nbeam = 500;
            mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; % P4 must be the same as P1 for Vortex.
            mp.P1.wGap = 0.01; % Fractional gap width
            mp.P4.padFacPct = 0; 
            mp = falco_gen_pupil_LUVOIR_B_with_phase(mp);

            %%- segmented mirror errors
            numSegments = hexSegMirror_numSegments(4); % Number of segments in "full" hex aperture
            % LUVOIR B has four rings, but ignores some corner segments 


            %% ZWFS Mask Properties
            mp.wfs.lambda0 = 425e-9;%--Central wavelength of the whole spectral bandpass [meters]
            mp.wfs.fracBW = 0.10; %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
            mp.wfs.Nsbp = 3; %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control

            %%- ZWFS mask properties 
            mp.wfs.mask.type = 'transmissive';
            maskDepth_m  = 213e-9;%m
            maskMaterial = 'FS';%Fused Silica 
            maskRadius_lamOverD = 1.06/2;%(maskDiam_m/2)/fNum/mp.lambda0;

            mp.wfs.mask.material = maskMaterial; % Required for transmissive mask 
            mp.wfs.mask.Rin  = maskRadius_lamOverD;% Radius of the ZWFS dimple 
            mp.wfs.mask.Rout = Inf;% Outer 
            mp.wfs.mask.depth = maskDepth_m; % Depth of the Zernike dimple
            mp.wfs.mask.FPMampFac = 1.0;% Transmission of the Zernike dimple
            mp.wfs.mask.res = 20;

            %%- ZWFS camera properties 
            mp.wfs.cam.Nbeam = mp.P1.full.Nbeam;% beam size at WFS camera 
            mp.wfs.cam.Narr = 512; % array size at WFS camera 

            mp.wfs.cam.Npix = 512; % Number of camera pixels across the Narr samples 
            mp.wfs.cam.centerPixOffset = [0,0];% Offset of camera wrt beam center (pixels) 


            %% Flesh out ws
            [mp, out] = falco_flesh_out_workspace(mp);

            %% Generate ZWFS images 
            rng(1)

            disp('***** ZWFS demo *****');
            disp('Generating ZWFS calib image...');

            Ical = falco_zwfs_getCalibrationImage(mp);
            b = falco_zwfs_getReferenceWave(mp);
            IZ0 = falco_zwfs_sim_image(mp);

            % IZ image mask 
            mask = imerode(logical(mp.P1.full.mask),strel('disk', 1));
            mask = padOrCropEven(mask,mp.wfs.cam.Narr);


            %%-- Apply first set of errors
            disp('Applying WFE to primary mirror...');
            mp.P1.pistons = randn(1,numSegments)/100;% Segment piston in waves 
            mp.P1.tiltxs  = randn(1,numSegments)/50;% %Tilts on segments in horiz direction (waves/apDia)
            mp.P1.tiltys  = randn(1,numSegments)/50;% %Tilts on segments in vert direction (waves/apDia)

            mp = falco_gen_pupil_LUVOIR_B_with_phase(mp);
            mp = falco_compute_entrance_pupil_coordinates(mp);

            actual_phz1 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

            disp('Generating ZWFS image...');
            IZ1 = falco_zwfs_sim_image(mp);

            %%-- Apply second set of errors
            disp('Applying new WFE to primary mirror...');
            mp.P1.pistons = mp.P1.pistons + randn(1,numSegments)/2000;% Segment piston in waves 
            mp.P1.tiltxs  = mp.P1.tiltxs + randn(1,numSegments)/1000;% %Tilts on segments in horiz direction (waves/apDia)
            mp.P1.tiltys  = mp.P1.tiltys + randn(1,numSegments)/1000;% %Tilts on segments in vert direction (waves/apDia)

            mp = falco_gen_pupil_LUVOIR_B_with_phase(mp);
            mp = falco_compute_entrance_pupil_coordinates(mp);

            actual_phz2 = angle(mp.P1.compact.E(:,:,ceil(mp.Nsbp/2)));

            disp('Generating ZWFS image...');
            IZ2 = falco_zwfs_sim_image(mp);

            %% Reconstruct phases

            disp('Reconstructing the wavefront...');
            theta = 2*pi*mp.wfs.mask.depth*(mp.wfs.mask.n(mp.wfs.lambda0)-1)/mp.wfs.lambda0;
            phz1 = falco_zwfs_reconstructor(Ical, IZ1, mask, b, theta, 'f','subbias');
            phz2 = falco_zwfs_reconstructor(Ical, IZ2, mask, b, theta, 'f','subbias');

            phz1 = circshift(rot90(phz1,2),[1 1]);
            phz2 = circshift(rot90(phz2,2),[1 1]);
            phz1(~mask) = 0;
            phz2(~mask) = 0;

            diffphz_meas = phz1 - phz2;
            diffphz_meas(~mask) = NaN;

            diffphz_actual = actual_phz1 - actual_phz2;
            diffphz_actual = padOrCropEven(diffphz_actual,length(diffphz_meas));
            diffphz_actual = diffphz_actual - mean(diffphz_actual(mask));
            diffphz_actual(~mask) = NaN;

            %% result

            rmsError = falco_rms(diffphz_actual(mask)*mp.lambda0-diffphz_meas(mask)*mp.wfs.lambda0)/4/pi; % meters
            
            testCase.verifyLessThan(rmsError, 50e-12); % less than 50 picometers sensing error.
        end
    end
end