classdef TestDownstreamImageShiftCompact < matlab.unittest.TestCase
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
            addpath(genpath([testCase.mp.path.falco filesep 'lib']));
            %[testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
        end
        % function fleshOut
        % end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath([testCase.mp.path.falco filesep 'lib']))
        end
    end   

%% Tests
%  - Backwards compatible--still runs with nothing defined.
%  - Compact, Nsbp = 1, Nwpsbp = 1: x is shifted by 1 lambda/D, y is shifted by -2 lambda/D.
%  - Compact, Nsbp = 3, Nwpsbp = 1: Works for scalar shifts
%  - Compact, Nsbp = 3, Nwpsbp = 1: Works for vector of shifts
%  - Compact, Nsbp = 3, Nwpsbp = 3: Works for scalar shifts
%  - Compact, Nsbp = 3, Nwpsbp = 3: Works for vector of shifts

    methods (Test)

        function testBackwardsCompatible(testCase)
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            testCase.mp.runLabel = 'test_backwards_compatible_compact';

            modvar = ModelVariables();
            % modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model
            modvar.sbpIndex = 1;
            modvar.starIndex = 1;
            modvar.whichSource = 'star';

            Eout = model_compact(testCase.mp, modvar);
            image0 = abs(Eout).^2;
            
            testCase.verifyTrue(true) % Got to this point
        end


        function testFullNsbp1Nwpsbp1(testCase)
            testCase.mp.Nsbp = 1;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            % testCase.mp = falco_set_spectral_properties(testCase.mp);
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';
            
            modvar = ModelVariables();
            modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model
            modvar.starIndex = 1;
            modvar.whichSource = 'star';

            for iSubband = 1:testCase.mp.Nsbp

                modvar.sbpIndex = iSubband;

                % Unshifted image
                testCase.mp.Fend.compact.xiOffset = 0;
                testCase.mp.Fend.compact.etaOffset = 0;
                % image0 = falco_get_sbp_image(testCase.mp, iSubband);
                E0 = model_compact(testCase.mp, modvar);
                image0 = abs(E0).^2;
                image0 = pad_crop(image0, Ncrop);

                % Shifted image
                testCase.mp.Fend.compact.xiOffset = 4/5; % lambda0/D
                testCase.mp.Fend.compact.etaOffset = -2; % lambda0/D
                E1 = model_compact(testCase.mp, modvar);
                image1a = abs(E1).^2;
                image1b = circshift(image1a, -testCase.mp.Fend.res*[testCase.mp.Fend.compact.etaOffset, testCase.mp.Fend.compact.xiOffset]);
                image1b = pad_crop(image1b, Ncrop);

                absDiff = abs(image0 - image1b);
                testCond = all(absDiff(:) <= atol);

                % figure(1); imagesc(log10(image0)); axis xy equal tight; colorbar;
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

            modvar = ModelVariables();
            modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model
            modvar.starIndex = 1;
            modvar.whichSource = 'star';
            
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            for iSubband = 1:testCase.mp.Nsbp

                modvar.sbpIndex = iSubband;

                % Unshifted image
                testCase.mp.Fend.compact.xiOffset = 0;
                testCase.mp.Fend.compact.etaOffset = 0;
                E0 = model_compact(testCase.mp, modvar);
                image0 = abs(E0).^2;
                image0 = pad_crop(image0, Ncrop);
    
                % Shifted image
                testCase.mp.Fend.compact.xiOffset = 4/5; % lambda0/D
                testCase.mp.Fend.compact.etaOffset = -2; % lambda0/D
                E1 = model_compact(testCase.mp, modvar);
                image1a = abs(E1).^2;
                image1b = circshift(image1a, -testCase.mp.Fend.res*[testCase.mp.Fend.compact.etaOffset, testCase.mp.Fend.compact.xiOffset]);
                image1b = pad_crop(image1b, Ncrop);
    
                absDiff = abs(image0 - image1b);
                testCond = all(absDiff(:) <= atol);
    
                testCase.verifyTrue(testCond)
            end
             
        end


        function testFullNsbp3Nwpsbp1Vector(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % testCase.mp = falco_set_spectral_properties(testCase.mp);
            % testCase.mp = falco_compute_psf_norm_factor(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            modvar = ModelVariables();
            modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model
            modvar.starIndex = 1;
            modvar.whichSource = 'star';

            for iSubband = 1:testCase.mp.Nsbp

                modvar.sbpIndex = iSubband;
                
                % Unshifted image
                testCase.mp.Fend.compact.xiOffset = 0;
                testCase.mp.Fend.compact.etaOffset = 0;
                E0 = model_compact(testCase.mp, modvar);
                image0 = abs(E0).^2;
                image0 = pad_crop(image0, Ncrop);
    
                % Shifted image
                testCase.mp.Fend.compact.xiOffset = 4/5*ones(testCase.mp.Nsbp, 1); % lambda0/D
                testCase.mp.Fend.compact.etaOffset = -2*ones(testCase.mp.Nsbp, 1); % lambda0/D
                E1 = model_compact(testCase.mp, modvar);
                image1a = abs(E1).^2;
                image1b = circshift(image1a, round(-testCase.mp.Fend.res*[testCase.mp.Fend.compact.etaOffset(1), testCase.mp.Fend.compact.xiOffset(1)]));
                image1b = pad_crop(image1b, Ncrop);
    
                absDiff = abs(image0 - image1b);
                testCond = all(absDiff(:) <= atol);
    
                testCase.verifyTrue(testCond)
            end
             
        end

        function testFullNsbp3Nwpsbp1VectorDifferentValues(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 1;
            testCase.mp.fracBW = 0.20;
            
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % testCase.mp = falco_set_spectral_properties(testCase.mp);
            % testCase.mp = falco_compute_psf_norm_factor(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            modvar = ModelVariables();
            modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model
            modvar.starIndex = 1;
            modvar.whichSource = 'star';

             % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1

                modvar.sbpIndex = iSubband;

                testCase.mp.Fend.compact.xiOffset = 0;
                testCase.mp.Fend.compact.etaOffset = 0;
                E0 = model_compact(testCase.mp, modvar);
                image0 = abs(E0).^2;
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;
            end

            % Shifted image
            testCase.mp.Fend.compact.xiOffset = [4/5, 2/5, 6/5]; % lambda0/D
            testCase.mp.Fend.compact.etaOffset = [0, -2, 2]; % lambda0/D
            for iSubband = testCase.mp.Nsbp:-1:1

                modvar.sbpIndex = iSubband;

                E1 = model_compact(testCase.mp, modvar);
                image1a = abs(E1).^2;
                image1b = circshift(image1a, round(-testCase.mp.Fend.res*[testCase.mp.Fend.compact.etaOffset(iSubband), testCase.mp.Fend.compact.xiOffset(iSubband)]));
                image1b = pad_crop(image1b, Ncrop);
   
                absDiff = abs(image0_cube(:, :, iSubband) - image1b);

                % figure(1); imagesc(log10(image0_cube(:, :, iSubband))); axis xy equal tight; colorbar;
                % figure(2); imagesc(log10(image1a)); axis xy equal tight; colorbar;
                % figure(3); imagesc(log10(image1b)); axis xy equal tight; colorbar;
                % figure(4); imagesc(absDiff); axis xy equal tight; colorbar;
                % drawnow;
                % pause(3)

                testCond = all(absDiff(:) <= atol);
                testCase.verifyTrue(testCond)
            end
             
        end


        function testFullNsbp3Nwpsbp3VectorDifferentValues(testCase)
            testCase.mp.Nsbp = 3;
            testCase.mp.Nwpsbp = 3;
            testCase.mp.fracBW = 0.20;
            
            [testCase.mp, ~] = falco_flesh_out_workspace(testCase.mp);
            % testCase.mp = falco_set_spectral_properties(testCase.mp);
            % testCase.mp = falco_compute_psf_norm_factor(testCase.mp);

            Ncrop = ceil_even(2*testCase.mp.Fend.res*(testCase.mp.Fend.FOV-2));
            atol = 1e-16;

            testCase.mp.runLabel = 'test_ds_image_shift';

            testCase.mp.Fend.compact.xiOffset = 4/5*ones(3, 1); % lambda0/D
            testCase.mp.Fend.compact.etaOffset = -2*ones(3, 1); % lambda0/D

            modvar = ModelVariables();
            modvar.wpsbpIndex = -1; %--Dummy index since not needed in compact model
            modvar.starIndex = 1;
            modvar.whichSource = 'star';

             % Unshifted image
            for iSubband = testCase.mp.Nsbp:-1:1

                modvar.sbpIndex = iSubband;

                % testCase.mp.Fend.compact.xiOffset = 0;
                % testCase.mp.Fend.compact.etaOffset = 0;
                % image0 = falco_get_summed_image(testCase.mp);
                E0 = model_compact(testCase.mp, modvar);
                image0 = abs(E0).^2; 
                image0 = pad_crop(image0, Ncrop);
                image0_cube(:, :, iSubband) = image0;
            end

            % Shifted image
            testCase.mp.Fend.compact.xiOffset = 4/5*ones(3, 1); % lambda0/D
            testCase.mp.Fend.compact.etaOffset = -2*ones(3, 1); % lambda0/D
            for iSubband = testCase.mp.Nsbp:-1:1

                modvar.sbpIndex = iSubband;

                E1 = model_compact(testCase.mp, modvar);
                image1a = abs(E1).^2;
                %image1a = falco_get_summed_image(testCase.mp);
                %image1b = circshift(image1a, round(-testCase.mp.Fend.res*[testCase.mp.Fend.compact.etaOffset(iSubband), testCase.mp.Fend.compact.xiOffset(iSubband)]));
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
