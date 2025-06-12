classdef TestDownstreamImageShiftFull < matlab.unittest.TestCase
%% Properties
%
% A presaved file with FALCO parameters was saved and is load to be used
% by methods. In this case we only use the mp.path.falco + lib/utils to
% addpath to utils functions to be tested.
    properties
        mp=ConfigurationLC();
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

%% Tests
%  - Backwards compatible--still runs with nothing defined.
%  - Full, Nsbp = 1, Nwpsbp = 1: x is shifted by 1 lambda/D, y is shifted by -2 lambda/D.
%  - Full, Nsbp = 3, Nwpsbp = 1: Works for scalar shifts; NrelayFend 0 or 1.
%  - Full, Nsbp = 3, Nwpsbp = 1: Works for vector of shifts
%  - Full, Nsbp = 3, Nwpsbp = 3: Works for scalar shifts
%  - Full, Nsbp = 3, Nwpsbp = 3: Works for vector of shifts

    methods (Test)

        function testBackwardsCompatible(testCase)
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            testCase.mp.runLabel = 'test_backwards_compatible';
            
            image0 = falco_get_sbp_image(testCase.mp, 1);
            image1 = falco_get_summed_image(testCase.mp);
            testCase.verifyTrue(true) % Got to this point
        end

        function testFullNsbp1Nwpsbp1NrelayFend0(testCase)
            testCase.mp.Nsbp = 1;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.flagRotation = true;
            testCase.mp.NrelayFend = 0;

            testCase.mp.Fend.full.xiOffset = 0;
            testCase.mp.Fend.full.etaOffset = 0;
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';
            
            % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1
                

                image0 = falco_get_sbp_image(testCase.mp, iSubband);
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;

            end

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = 4/5; % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2; % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            for iSubband = testCase.mp.Nsbp:-1:1

                image1a = falco_get_sbp_image(testCase.mp, iSubband);
                image1b = circshift(image1a, -testCase.mp.Fend.res*[testCase.mp.Fend.full.etaOffset, testCase.mp.Fend.full.xiOffset]);
                image1b = pad_crop(image1b, Ncrop);

                absDiff = abs(image0_cube(:, :, iSubband) - image1b);
                testCond = all(absDiff(:) <= atol);

                % figure(1); imagesc(log10(image0_cube(:, :, iSubband))); axis xy equal tight; colorbar;
                % figure(2); imagesc(log10(image1a)); axis xy equal tight; colorbar;
                % figure(3); imagesc(log10(image1b)); axis xy equal tight; colorbar;
                % figure(4); imagesc(absDiff); axis xy equal tight; colorbar;
                % drawnow;

                testCase.verifyTrue(testCond)
            end
        end

        function testFullNsbp1Nwpsbp1NrelayFend1(testCase)
            testCase.mp.Nsbp = 1;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.flagRotation = true;
            testCase.mp.NrelayFend = 1;

            testCase.mp.Fend.full.xiOffset = 0;
            testCase.mp.Fend.full.etaOffset = 0;
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';
            
            % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1
                

                image0 = falco_get_sbp_image(testCase.mp, iSubband);
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;

            end

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = 4/5; % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2; % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            for iSubband = testCase.mp.Nsbp:-1:1

                image1a = falco_get_sbp_image(testCase.mp, iSubband);
                image1b = circshift(image1a, -testCase.mp.Fend.res*[testCase.mp.Fend.full.etaOffset, testCase.mp.Fend.full.xiOffset]);
                image1b = pad_crop(image1b, Ncrop);

                absDiff = abs(image0_cube(:, :, iSubband) - image1b);
                testCond = all(absDiff(:) <= atol);

                % figure(1); imagesc(log10(image0_cube(:, :, iSubband))); axis xy equal tight; colorbar;
                % figure(2); imagesc(log10(image1a)); axis xy equal tight; colorbar;
                % figure(3); imagesc(log10(image1b)); axis xy equal tight; colorbar;
                % figure(4); imagesc(absDiff); axis xy equal tight; colorbar;
                % drawnow;

                testCase.verifyTrue(testCond)
            end
        end

        function testFullNsbp3Nwpsbp1Scalar(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.Fend.full.xiOffset = 0;
            testCase.mp.Fend.full.etaOffset = 0;
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % testCase.mp = falco_set_spectral_properties(testCase.mp);
            % testCase.mp = falco_compute_psf_norm_factor(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            % Unshifted image
            image0 = falco_get_summed_image(testCase.mp);
            % image0 = falco_get_sbp_image(testCase.mp, iSubband);
            image0 = pad_crop(image0, Ncrop);

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = 4/5; % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2; % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % image1a = falco_get_sbp_image(testCase.mp, iSubband);
            image1a = falco_get_summed_image(testCase.mp);
            image1b = circshift(image1a, -testCase.mp.Fend.res*[testCase.mp.Fend.full.etaOffset, testCase.mp.Fend.full.xiOffset]);
            image1b = pad_crop(image1b, Ncrop);

            absDiff = abs(image0 - image1b);
            testCond = all(absDiff(:) <= atol);

            testCase.verifyTrue(testCond)
             
        end

        function testFullNsbp3Nwpsbp1Vector(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.Fend.full.xiOffset = 0;
            testCase.mp.Fend.full.etaOffset = 0;
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % testCase.mp = falco_set_spectral_properties(testCase.mp);
            % testCase.mp = falco_compute_psf_norm_factor(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            % Unshifted image
            image0 = falco_get_summed_image(testCase.mp);
            % image0 = falco_get_sbp_image(testCase.mp, iSubband);
            image0 = pad_crop(image0, Ncrop);

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = 4/5*ones(testCase.mp.Nsbp, 1); % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2*ones(testCase.mp.Nsbp, 1); % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % image1a = falco_get_sbp_image(testCase.mp, iSubband);
            image1a = falco_get_summed_image(testCase.mp);
            image1b = circshift(image1a, round(-testCase.mp.Fend.res*[testCase.mp.Fend.full.etaOffset(1), testCase.mp.Fend.full.xiOffset(1)]));
            image1b = pad_crop(image1b, Ncrop);

            absDiff = abs(image0 - image1b);
            testCond = all(absDiff(:) <= atol);

            testCase.verifyTrue(testCond)
             
        end

        function testFullNsbp3Nwpsbp1VectorSameValue(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.Fend.full.xiOffset = 0;
            testCase.mp.Fend.full.etaOffset = 0;
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1
                image0 = falco_get_sbp_image(testCase.mp, iSubband);
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;
            end

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = 4/5*ones(testCase.mp.Nsbp, 1); % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2*ones(testCase.mp.Nsbp, 1); % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            for iSubband = testCase.mp.Nsbp:-1:1
                image1a = falco_get_sbp_image(testCase.mp, iSubband);
                image1b = circshift(image1a, round(-testCase.mp.Fend.res*[testCase.mp.Fend.full.etaOffset(iSubband), testCase.mp.Fend.full.xiOffset(iSubband)]));
                image1b = pad_crop(image1b, Ncrop);
   
                absDiff = abs(image0_cube(:, :, iSubband) - image1b);
                testCond = all(absDiff(:) <= atol);
                testCase.verifyTrue(testCond)
            end
             
        end

        function testFullNsbp3Nwpsbp1VectorDifferentValues(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.Fend.full.xiOffset = 0;
            testCase.mp.Fend.full.etaOffset = 0;
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1
                image0 = falco_get_sbp_image(testCase.mp, iSubband);
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;
            end

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = [4/5, 2/5, 6/5]; % lambda0/D
            testCase.mp.Fend.full.etaOffset = [0, -2, 2]; % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            for iSubband = testCase.mp.Nsbp:-1:1
                image1a = falco_get_sbp_image(testCase.mp, iSubband);
                image1b = circshift(image1a, round(-testCase.mp.Fend.res*[testCase.mp.Fend.full.etaOffset(iSubband), testCase.mp.Fend.full.xiOffset(iSubband)]));
                image1b = pad_crop(image1b, Ncrop);
   
                absDiff = abs(image0_cube(:, :, iSubband) - image1b);
                testCond = all(absDiff(:) <= atol);
                testCase.verifyTrue(testCond)
            end
             
        end

        function testFullNsbp3Nwpsbp3VectorDifferentValues(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 3;
            testCase.mp.fracBW = 0.20;
            
            testCase.mp.Fend.full.xiOffset = 4/5*ones(3, 1); % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2*ones(3, 1); % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
 
            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1
                image0 = falco_get_sbp_image(testCase.mp, iSubband);
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;
            end

            % Shifted image
            testCase.mp.P4.compact = rmfield(testCase.mp.P4.compact, 'E');
            testCase.mp.P4.full = rmfield(testCase.mp.P4.full, 'E');
            testCase.mp.Fend.full.xiOffset = 4/5*ones(7, 1); % lambda0/D
            testCase.mp.Fend.full.etaOffset = -2*ones(7, 1); % lambda0/D
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            for iSubband = testCase.mp.Nsbp:-1:1
                image1a = falco_get_sbp_image(testCase.mp, iSubband);
                image1b = pad_crop(image1a, Ncrop);
                absDiff = abs(image0_cube(:, :, iSubband) - image1b);

                % figure(1); imagesc(log10(image0_cube(:, :, iSubband))); axis xy equal tight; colorbar;
                % figure(2); imagesc(log10(image1a)); axis xy equal tight; colorbar;
                % figure(3); imagesc(log10(image1b)); axis xy equal tight; colorbar;
                % figure(4); imagesc(absDiff); axis xy equal tight; colorbar;
                % drawnow;

                testCond = all(absDiff(:) <= atol);
                testCase.verifyTrue(testCond)
            end
             
        end


    end    
end
