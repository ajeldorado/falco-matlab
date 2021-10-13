%% Test input checks
%
% Unit tests of the methods in the Check class.
%
classdef TestConstrainDM < matlab.unittest.TestCase
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    properties
        increasing_data = reshape(0:24, [5, 5]).';
        data_to_fail = ones(5) + eye(5);
        data_ones = ones(5);
        data_end = zeros(12, 12);
        bigeye = 5.0 * eye(5);
        tie_no_tied_dead = zeros(3, 3);
        volts_no_tied_dead = reshape(1:9, [3, 3]).';
        tie_tied = [0, 1, 0; 0, 1, 3; 2, 2, 0];
        volts_tied = [1, 2, 3; 4, 2, 6; 7, 7, 9];
        tie_dead = [-1, 0, -1; 0, 0, 0; 0, 0, 0];
        volts_dead = [0, 2, 0; 4, 5, 6; 7, 8, 9];
        tie_tied_dead = [-1, 1, -1; 0, 1, 3; 2, 2, 0];
        volts_tied_dead = [0, 2, 0; 4, 2, 6; 7, 7, 9];
        
        rng_seed = 2021;
        nact = 48;
        flatmap_flat = zeros(48);
        vmin = 0.;
        vmax = 100.;
        vmid = 50
        vlat = 50.;
        vdiag = 75.;
        vquant = 100 / 2^16;
        maxiter = 1000;
        
        vneighbor = 30.;
        vcorner = 30.;
        dm0 = 50*ones(48);
        uniform_flat = 50*ones(48);
    end

    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath('../../lib'));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath('../../lib'));
        end
    end
    
%% Unit tests of ConstrainDM

    methods (Test)
                
        % Unit tests of ConstrainDM.check_index_validity
        function test_index_inside(testCase)
            testCase.verifyTrue(ConstrainDM.check_index_validity([2, 3], 2, 3))
        end
        function test_index_outside(testCase)
            testCase.verifyFalse(ConstrainDM.check_index_validity([3, 3], 2, 3))
        end
        
        % Unit tests of ConstrainDM.neighbor_values_plus
        function testPlus33(testCase)
            data = testCase.increasing_data;
            out = sort(ConstrainDM.neighbor_values_plus(data, 3, 3));
            outExpected = sort([7, 13, 17, 11]);
            testCase.verifyTrue(all(out == outExpected))
        end
        function testPlus11(testCase)
            data = testCase.increasing_data;
            out = sort(ConstrainDM.neighbor_values_plus(data, 1, 1));
            outExpected = sort([5, 1]);
            testCase.verifyTrue(all(out == outExpected))
        end
        function testPlus31(testCase)
            data = testCase.increasing_data;
            out = sort(ConstrainDM.neighbor_values_plus(data, 3, 1));
            outExpected = sort([5, 11, 15]);
            testCase.verifyTrue(all(out == outExpected))
        end
        
        % Unit tests of ConstrainDM.neighbor_values_diag
        function testDiag33(testCase)
            data = testCase.increasing_data;
            out = sort(ConstrainDM.neighbor_values_diag(data, 3, 3));
            outExpected = sort([6, 16, 8, 18]);
            testCase.verifyTrue(all(out == outExpected))
        end
        function testDiag11(testCase)
            data = testCase.increasing_data;
            out = sort(ConstrainDM.neighbor_values_diag(data, 1, 1));
            outExpected = 6;
            testCase.verifyTrue(all(out == outExpected))
        end
        function testDiag31(testCase)
            data = testCase.increasing_data;
            out = sort(ConstrainDM.neighbor_values_diag(data, 3, 1));
            outExpected = sort([6, 16]);
            testCase.verifyTrue(all(out == outExpected))
        end
        
        % Unit tests of ConstrainDM.check_high_low_limit
        function testCheckHighLowLimit_tooHigh(testCase)
            data = testCase.increasing_data;
            testCase.verifyFalse(ConstrainDM.check_high_low_limit(...
                data, max(data(:))-1, min(data(:))))
        end
        function testCheckHighLowLimit_tooLow(testCase)
            data = testCase.increasing_data;
            testCase.verifyFalse(ConstrainDM.check_high_low_limit(...
                data, max(data(:)), min(data(:))+1))
        end
        function testCheckHighLowLimit_withinBounds(testCase)
            data = testCase.increasing_data;
            testCase.verifyTrue(ConstrainDM.check_high_low_limit(...
                data, max(data(:)), min(data(:))))
        end
        
        % Unit tests of ConstrainDM.check_valid_neighbors
        function test_check_valid_neighbors_equal_pass(testCase)
            % Uniform array passes as expected with no flatmap.
            data = testCase.data_ones;
            plus_limit = 1;
            diag_limit = 1;
            high_limit = [];
            low_limit = [];
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(data, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        function test_check_valid_neighbors_pass_dmflat(testCase)
            % An array that should fail NR checks passes if the dmflat is the
            % same (i.e. no penalty for making surface phase-flat).
            data = testCase.bigeye;
            plus_limit = 1;
            diag_limit = 1;
            high_limit = [];
            low_limit = [];
            dmflat = testCase.bigeye;
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(data, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        function test_check_valid_neighbors_fail_bounds_high(testCase)
            data = testCase.data_ones;
            plus_limit = 1;
            diag_limit = 1;
            high_limit = 0.5;
            low_limit = [];
            dmflat = [];
            testCase.verifyFalse(ConstrainDM.check_valid_neighbors(data, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        function test_check_valid_neighbors_fail_bounds_low(testCase)
            data = testCase.data_ones;
            plus_limit = 1;
            diag_limit = 1;
            high_limit = [];
            low_limit = 1.5;
            dmflat = [];
            testCase.verifyFalse(ConstrainDM.check_valid_neighbors(data, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        
        % Unit tests of ConstrainDM.check_tie_dead
        function test_success_no_tied_dead(testCase)
            % Succeeds for known-good data.
            testCase.verifyTrue(ConstrainDM.check_tie_dead(testCase.volts_no_tied_dead, testCase.tie_no_tied_dead))
        end
        function test_success_tied(testCase)
            % Succeeds for known-good data.
            testCase.verifyTrue(ConstrainDM.check_tie_dead(testCase.volts_tied, testCase.tie_tied))
        end
        function test_success_dead(testCase)
            % Succeeds for known-good data"""
            testCase.verifyTrue(ConstrainDM.check_tie_dead(testCase.volts_dead, testCase.tie_dead))
        end
        function test_success_tied_dead(testCase)
            % Succeeds for known-good data"""
            testCase.verifyTrue(ConstrainDM.check_tie_dead(testCase.volts_tied_dead, testCase.tie_tied_dead))
        end
        function test_failure_tied(testCase)
            % Fails when tied actuators do not match"""
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_no_tied_dead, testCase.tie_tied))
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_dead, testCase.tie_tied))
        end
        function test_failure_dead(testCase)
            % Fails when dead actuators are not zero.
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_no_tied_dead, testCase.tie_dead))
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_tied, testCase.tie_dead))
        end
        function test_failure_tied_dead(testCase)
            % Fails when tied actuators do not match or dead actuators are not zero.
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_no_tied_dead, testCase.tie_tied_dead))
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_dead, testCase.tie_tied_dead))
            testCase.verifyFalse(ConstrainDM.check_tie_dead(testCase.volts_tied, testCase.tie_tied_dead))
        end
        function test_volts_and_tie_different_sizes(testCase)
            identifier = 'ValueError:SizeMismatch';
            verifyError(testCase, @() ConstrainDM.check_tie_dead(testCase.volts_no_tied_dead(2:end, 2:end),...
                testCase.tie_no_tied_dead), identifier)
        end
        
        
        % Unit tests of ConstrainDM.checktie
        function test_tie_not_made_of_correct_values(testCase)
            % Bad inputs caught as expected.
            % expect consecutive integers 1 -> N, along with 0 or -1
            testCase.verifyFalse(ConstrainDM.checktie([0.5, 0.5; 0.5, 0.5])) % not integers
            testCase.verifyFalse(ConstrainDM.checktie([1, 2; 3, 5])) % not consecutive
            testCase.verifyFalse(ConstrainDM.checktie([0, 1; 2, 4])) % not consecutive w/ 0
            testCase.verifyFalse(ConstrainDM.checktie([-1, 1; 2, 4])) % not consecutive w/ -1
            testCase.verifyFalse(ConstrainDM.checktie([0, -1; 1, 3])) % not consecutive w/ 0 and -1
            testCase.verifyFalse(ConstrainDM.checktie([1, 2; 3, -2])) % not consecutive (negative)
            testCase.verifyFalse(ConstrainDM.checktie([2, 3; 4, 5])) % consecutive but not from 0
        end
        function test_tie_is_made_of_correct_values(testCase)
            % Good inputs allowed through as expected.
            % expect consecutive integers 1 -> N, along with 0 or -1
            testCase.verifyTrue(ConstrainDM.checktie([0, 0; 0, 0])) % all 0
            testCase.verifyTrue(ConstrainDM.checktie([-1, -1; -1, -1])) % all -1
            testCase.verifyTrue(ConstrainDM.checktie([0, -1; -1, 0])) % 0 and -1
            testCase.verifyTrue(ConstrainDM.checktie([1, 1; 1, 1])) % all 1
            testCase.verifyTrue(ConstrainDM.checktie([1, 2; 1, 3])) % consecutive with a group
            testCase.verifyTrue(ConstrainDM.checktie([1, 2; 1, 2])) % consecutive with groups
            testCase.verifyTrue(ConstrainDM.checktie([1, 3; 2, 4])) % consecutive but all individ.
            testCase.verifyTrue(ConstrainDM.checktie([0, 1; 2, 2])) % 1->N with 0
            testCase.verifyTrue(ConstrainDM.checktie([-1, 1; 2, 2])) % 1->N with -1
            testCase.verifyTrue(ConstrainDM.checktie([0, 1; -1, 2])) % 1->N with 0 and -1
        end

        
        % Unit tests of ConstrainDM.checkflat
        function test_good(testCase)
            % Verify good flats return True as expected.
            vmin = 0.0;
            vmax = 100.0;
            tie = [0, 1, 0; 0, 0, 1; -1, 0, 0];
            flatmap = [0, 2, 100; 4, 5, 2; 0, 8, 9];
            testCase.verifyTrue(ConstrainDM.checkflat(flatmap, vmin, vmax, tie))
        end
        function test_bad(testCase)
            % Verify bad flats return False as expected"""
            % all should be >= 0, <= vmax, all ties at same voltage, and all
            % dead actuator at 0V

            vmin = 0.0;
            vmax = 100.0;
            tie = [0, 1, 0; 0, 0, 1; -1, 0, 0];

            % vmin
            flat1 = [vmin-1, 0, 0; 0, 0, 0; 0, 0, 0];
            testCase.verifyFalse(ConstrainDM.checkflat(flat1, vmin, vmax, tie))


            % vmax
            flat2 = [0, 0, vmax+1; 0, 0, 0; 0, 0, 0];
            testCase.verifyFalse(ConstrainDM.checkflat(flat2, vmin, vmax, tie))

            % ties same V
            flat3 = [0, 1, 0; 0, 0, 0; 0, 0, 0];
            testCase.verifyFalse(ConstrainDM.checkflat(flat3, vmin, vmax, tie))

            % dead @ 0V
            flat4 = [0, 0, 0; 0, 0, 0; 1, 0, 0];
            testCase.verifyFalse(ConstrainDM.checkflat(flat4, vmin, vmax, tie))
        end
        
 
        %% Unit tests of ConstrainDM.dmsmooth
        function test_success_flat(testCase)
            % Runs with no issues on a flat input DM, and returns the same output
            dm0a = ConstrainDM.dmsmooth(testCase.dm0, testCase.vmax, testCase.vquant, testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.assertTrue(ConstrainDM.check_valid_neighbors(dm0a, plus_limit, diag_limit, high_limit, low_limit, dmflat))

            testCase.assertTrue(all(all((dm0a == testCase.dm0))))
        end
        
        function test_success_flat_with_uniform_flatmap(testCase)
            % Runs with no issues on a flat input DM, and returns the same output
            % with a uniform flatmap
            dm0a = ConstrainDM.dmsmooth(testCase.dm0, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, testCase.uniform_flat);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = testCase.uniform_flat;
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dm0a, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            testCase.verifyTrue(all(all(dm0a == testCase.dm0)))
        end
        
        function test_success_flat_with_nonuniform_flatmap(testCase)
            % Runs with no issues on a flat input DM, and returns the same output
            % with a non-uniform flatmap which doesn't violate any neighbor rules
            nonuniform_flat = testCase.vmid*ones(testCase.nact) + min([testCase.vneighbor, testCase.vcorner])/10*eye(testCase.nact);
            dmflat = nonuniform_flat;
            dm0a = ConstrainDM.dmsmooth(testCase.dm0, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = nonuniform_flat;
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dm0a, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            testCase.verifyTrue(all(all(dm0a == testCase.dm0)))
        end
        
        function test_success_stress_fixed_a(testCase)
            % Runs with no failure for stressing cases alternating between bounds.
            %
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            dmc = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmc, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_success_stress_fixed_b(testCase)
            % Runs with no failure for stressing case going upward and downward from a midpoint level.
            % 
            % Up-down checkerboard
            dmud = testCase.vmid * ones(testCase.nact);
            dmud(1:2:end, 1:2:end) = testCase.vmin;
            dmud(2:2:end, 2:2:end) = testCase.vmax;
        
            dmu = ConstrainDM.dmsmooth(dmud, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmu, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        
        function test_success_stress_fixed_noquant_a(testCase)
            % Runs with no failure for stressing cases alternating between
            % bounds when vquant is zero.
            %
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            
            dmc = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, 0, testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmc, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_success_stress_fixed_noquant_b(testCase)
            % Runs with no failure for stressing case going upward and
            % downward from a midpoint level when vquant is zero.
            % 
            % Up-down checkerboard
            dmud = testCase.vmid * ones(testCase.nact);
            dmud(1:2:end, 1:2:end) = testCase.vmin;
            dmud(2:2:end, 2:2:end) = testCase.vmax;
            
            dmu = ConstrainDM.dmsmooth(dmud, testCase.vmax, 0, testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmu, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        
        
        function test_success_stress_fixed_with_nonuniform_flatmap_a(testCase)
            % Runs with no failure for stressing cases alternating between
            % bounds in the presence of a nonuniform flatmap.
            
            nonuniform_flat = testCase.vmid*ones(testCase.nact) + min([testCase.vneighbor, testCase.vcorner])/10*eye(testCase.nact);
            dmflat = nonuniform_flat;
            
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            dmc = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmc, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_success_stress_fixed_with_nonuniform_flatmap_b(testCase)
            % Runs with no failure for stressing case going upward and
            % downward from a midpoint level in the presence of a
            % nonuniform flatmap.
            
            nonuniform_flat = testCase.vmid*ones(testCase.nact) + min([testCase.vneighbor, testCase.vcorner])/10*eye(testCase.nact);
            dmflat = nonuniform_flat;
            
            % Up-down checkerboard
            dmud = testCase.vmid * ones(testCase.nact);
            dmud(1:2:end, 1:2:end) = testCase.vmin;
            dmud(2:2:end, 2:2:end) = testCase.vmax;
        
            dmu = ConstrainDM.dmsmooth(dmud, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmu, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_success_stress_fuzz(testCase)
            % Runs with no failure for a set of stressing cases made from uniformly
            % distributing actuator heights between vmin and vmax
            % 
            % Checks both that operation succeeds and that output is valid
            % 
            % Seed is fixed, so these are reproducible; intent is write explicit
            % checks for the edge cases we can think of and use a randomly-derived
            % set to see if we come across anything we didn't think of.
            
            dmflat = [];
            
            % Fuzz
            rng(3621)
            nfuzz = 20;%200;
            dmfuzz = testCase.vmin + (testCase.vmax - testCase.vmin)*rand(testCase.nact, testCase.nact, nfuzz);
            
            for iSlice = 1:nfuzz
                dmf = squeeze(dmfuzz(:, :, iSlice));
                dmout = ConstrainDM.dmsmooth(dmf, testCase.vmax, testCase.vquant, ...
                    testCase.vneighbor, testCase.vcorner, dmflat);
                
                plus_limit=testCase.vneighbor;
                diag_limit=testCase.vcorner;
                high_limit=testCase.vmax;
                low_limit=testCase.vmin;
                testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmout, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            end
            
        end
        
        function test_success_stress_fuzz_with_nonuniform_flat(testCase)
            % Runs with no failure for a set of stressing cases made from uniformly
            % distributing actuator heights between vmin and vmax
            % 
            % Checks both that operation succeeds and that output is valid
            % 
            % Seed is fixed, so these are reproducible; intent is write explicit
            % checks for the edge cases we can think of and use a randomly-derived
            % set to see if we come across anything we didn't think of.
            
            nonuniform_flat = testCase.vmid*ones(testCase.nact) + min([testCase.vneighbor, testCase.vcorner])/10*eye(testCase.nact);
            dmflat = nonuniform_flat;
            
            % Fuzz
            rng(5599)
            nfuzz = 1;%20;%200;
            dmfuzz = testCase.vmin + (testCase.vmax - testCase.vmin)*rand(testCase.nact, testCase.nact, nfuzz);
            
            for iSlice = 1:nfuzz
                dmf = squeeze(dmfuzz(:, :, iSlice));
                dmout = ConstrainDM.dmsmooth(dmf, testCase.vmax, testCase.vquant, ...
                    testCase.vneighbor, testCase.vcorner, dmflat);
                
                plus_limit=testCase.vneighbor;
                diag_limit=testCase.vcorner;
                high_limit=testCase.vmax;
                low_limit=testCase.vmin;
                testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmout, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            end
            
        end
        
        function test_idempotence_fixed_a(testCase)
            % For fixed stressing case, the output from the first call does
            % not change when passed through again.
            
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            dmc1 = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner);
            dmc2 = ConstrainDM.dmsmooth(dmc1, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner);
            
            testCase.verifyTrue(isequal(dmc1, dmc2))
        end
        
        function test_idempotence_fixed_b(testCase)
            % For fixed stressing case, the output from the first call does
            % not change when passed through again.
            
            % Up-down checkerboard
            dmud = testCase.vmid * ones(testCase.nact);
            dmud(1:2:end, 1:2:end) = testCase.vmin;
            dmud(2:2:end, 2:2:end) = testCase.vmax;
        
            dmu1 = ConstrainDM.dmsmooth(dmud, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner);
            dmu2 = ConstrainDM.dmsmooth(dmu1, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner);
            testCase.verifyTrue(isequal(dmu1, dmu2))

        end
        
        function test_idempotence_fixed_with_nonuniform_flat_a(testCase)
            % For fixed stressing case, the output from the first call does
            % not change when passed through again.
            
            nonuniform_flat = testCase.vmid*ones(testCase.nact) + min([testCase.vneighbor, testCase.vcorner])/10*eye(testCase.nact);
            dmflat = nonuniform_flat;
            
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            
            dmc1 = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            dmc2 = ConstrainDM.dmsmooth(dmc1, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            
            testCase.verifyTrue(isequal(dmc1, dmc2))
        end
        
        function test_idempotence_fixed_with_nonuniform_flat_b(testCase)
            % For fixed stressing case, the output from the first call does
            % not change when passed through again.
            
            nonuniform_flat = testCase.vmid*ones(testCase.nact) + min([testCase.vneighbor, testCase.vcorner])/10*eye(testCase.nact);
            dmflat = nonuniform_flat;
            
            % Up-down checkerboard
            dmud = testCase.vmid * ones(testCase.nact);
            dmud(1:2:end, 1:2:end) = testCase.vmin;
            dmud(2:2:end, 2:2:end) = testCase.vmax;
        
            dmu1 = ConstrainDM.dmsmooth(dmud, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            dmu2 = ConstrainDM.dmsmooth(dmu1, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vcorner, dmflat);
            
            testCase.verifyTrue(isequal(dmu1, dmu2))

        end
        
        function test_idempotence_fuzz(testCase)
            % For randomly-derived stressing cases, the output from the first run,
            % when fed as input to a second run, produces the same result
            
            dmflat = [];
            
            % Fuzz
            rng(3621)
            nfuzz = 20;%200;
            dmfuzz = testCase.vmin + (testCase.vmax - testCase.vmin)*rand(testCase.nact, testCase.nact, nfuzz);
            
            for iSlice = 1:nfuzz
                dmf = squeeze(dmfuzz(:, :, iSlice));
                dmout1 = ConstrainDM.dmsmooth(dmf, testCase.vmax, testCase.vquant, ...
                    testCase.vneighbor, testCase.vcorner, dmflat);
                dmout2 = ConstrainDM.dmsmooth(dmout1, testCase.vmax, testCase.vquant, ...
                    testCase.vneighbor, testCase.vcorner, dmflat);
                
                testCase.verifyTrue(isequal(dmout1, dmout2))

            end
            
        end
        
        function test_margin_fixed(testCase)
            % Verify that all outputs are at least 2 LSB from their min/max, and all
            % all NR gaps are at least 2 LSB beyond the accompanying rule.

            margin = 2.*testCase.vquant;
            
            %TODO
            
        end
        
        function test_margin_fuzz(testCase)
            % Verify that all outputs are at least 2 LSB from their min/max, and all
            % all NR gaps are at least 2 LSB beyond the accompanying rule.

            margin = 2.*testCase.vquant;
            
            %TODO
            
        end
        
        function test_vneighbor_vcorner_differ_a(testCase)
            % Check checkerboard returns successfully if the neighbor rules are
            % different between diagonal and lateral.
            
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            dmc = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vneighbor+1);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vneighbor + 1;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmc, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_vneighbor_vcorner_differ_b(testCase)
            % Check checkerboard returns successfully if the neighbor rules are
            % different between diagonal and lateral.
            
            % Up-down checkerboard
            dmud = testCase.vmid * ones(testCase.nact);
            dmud(1:2:end, 1:2:end) = testCase.vmin;
            dmud(2:2:end, 2:2:end) = testCase.vmax;
        
            dmu = ConstrainDM.dmsmooth(dmud, testCase.vmax, testCase.vquant, ...
                testCase.vneighbor, testCase.vneighbor+1);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vneighbor + 1;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmu, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_nonr(testCase)
            % Check checkerboard returns successfully if the neighbor rules are
            % all larger than the valid voltage range (i.e. disabled)
            
            rangeplus = testCase.vmax - testCase.vmin + 1;
            
            % Checkerboard
            dmcheck = testCase.vmin * ones(testCase.nact);
            dmcheck(1:2:end, 1:2:end) = testCase.vmax;
            dmcheck(2:2:end, 2:2:end) = testCase.vmax;
            dmc = ConstrainDM.dmsmooth(dmcheck, testCase.vmax, testCase.vquant, ...
                rangeplus, rangeplus);
            
            plus_limit=rangeplus;
            diag_limit=rangeplus;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmc, plus_limit, diag_limit, high_limit, low_limit, dmflat))
        end
        
        function test_cap_high_noquant(testCase)
            % Verify high voltage is capped as expected (with no vquant).
            dmhigh = (testCase.vmax + 1) * ones(testCase.nact);
            dmcap = (testCase.vmax) * ones(testCase.nact);
            
            dmout = ConstrainDM.dmsmooth(dmhigh, testCase.vmax, 0, testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmout, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            
            testCase.verifyTrue(isequal(dmout, dmcap))
            
        end
        
        function test_cap_low_noquant(testCase)
            % Verify low voltage is capped as expected (with no vquant).
            dmlow = (testCase.vmin - 1) * ones(testCase.nact);
            dmcap = (testCase.vmin) * ones(testCase.nact);
            
            dmout = ConstrainDM.dmsmooth(dmlow, testCase.vmax, 0, testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmout, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            
            testCase.verifyTrue(isequal(dmout, dmcap))
            
        end
        
        function test_oor_checker(testCase)
            % Check checkerboard that extends out-of-range completes successfully.
            
            % Out-of-range checkerboard            
            dmoor = (testCase.vmin - 1) * ones(testCase.nact);
            dmoor(1:2:end, 1:2:end) = testCase.vmax + 1;
            dmoor(2:2:end, 2:2:end) = testCase.vmax + 1;
            
            dmc = ConstrainDM.dmsmooth(dmoor, testCase.vmax, testCase.vquant, testCase.vneighbor, testCase.vcorner);
            
            plus_limit=testCase.vneighbor;
            diag_limit=testCase.vcorner;
            high_limit=testCase.vmax;
            low_limit=testCase.vmin;
            dmflat = [];
            testCase.verifyTrue(ConstrainDM.check_valid_neighbors(dmc, plus_limit, diag_limit, high_limit, low_limit, dmflat))
            
        end
        
        function test_horizontal_noquant(testCase)
            % Verify analytical result for horizontal NR violations (no vquant).
            tol = 1e-13;

            dmnr = [20, 20, 20; 20, 40, 0; 20, 20, 20];
            dmexpect = [20, 20, 20; 20, 35, 5; 20, 20, 20];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_vertical_noquant(testCase)
            % Verify analytical result for vertical NR violations (no vquant).
            tol = 1e-13;

            dmnr = [20, 20, 20; 20, 40, 20; 20, 0, 20];
            dmexpect = [20, 20, 20; 20, 35, 20; 20, 5, 20];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_diag1_noquant(testCase)
            % Verify analytical result for 1st diag NR violations (no vquant).
            tol = 1e-13;

            dmnr = [20, 20, 20; 20, 40, 20; 20, 20, 0];
            dmexpect = [20, 20, 20; 20, 35, 20; 20, 20, 5];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_diag2_noquant(testCase)
            % Verify analytical result for 2nd diag NR violations (no vquant).
            tol = 1e-13;

            dmnr = [20, 20, 20; 20, 40, 20; 0, 20, 20];
            dmexpect = [20, 20, 20; 20, 35, 20; 5, 20, 20];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_multiNR_noquant(testCase)
            % Verify analytical result (no vquant) with more than 1 NR violation.
            tol = 1e-13;

            dmnr = [40, 0, 20; 20, 20, 0; 20, 40, 20];
            dmexpect = [35, 5, 20; 20, 20, 5; 20, 35, 20];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_check_odd_noquant(testCase)
            % Verify analytical result is correct for checkerboard pattern (nside odd).
            tol = 1e-13;

            dmnr = [0, 100, 0; 100, 0, 100; 0, 100, 0];
            dmexpect = [35, 65, 35; 65, 35, 65; 35, 65, 35];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_check_even_noquant(testCase)
            % Verify analytical result is correct for checkerboard pattern (nside even).
            tol = 1e-13;

            dmnr = [0, 100, 0, 100; 100, 0, 100, 0; 0, 100, 0, 100; 100, 0, 100, 0];
            dmexpect = [35, 65, 35, 65; 65, 35, 65, 35; 35, 65, 35, 65; 65, 35, 65, 35];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_hstripe_odd_noquant(testCase)
            % Verify analytical result is correct for horizontal stripe pattern (nside odd).
            tol = 1e-13;

            dmnr = [0, 0, 0; 100, 100, 100; 0, 0, 0];
            dmexpect = [35, 35, 35; 65, 65, 65; 35, 35, 35];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_hstripe_even_noquant(testCase)
            % Verify analytical result is correct for horizontal stripe pattern (nside even).
            tol = 1e-13;

            dmnr = [0, 0, 0, 0; 100, 100, 100, 100; 0, 0, 0, 0; 100, 100, 100, 100];
            dmexpect = [35, 35, 35, 35; 65, 65, 65, 65; 35, 35, 35, 35; 65, 65, 65, 65];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_vstripe_odd_noquant(testCase)
            % Verify analytical result is correct for vertical stripe pattern (nside odd).
            tol = 1e-13;

            dmnr = [0, 0, 0; 100, 100, 100; 0, 0, 0].';
            dmexpect = [35, 35, 35; 65, 65, 65; 35, 35, 35].';

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_vstripe_even_noquant(testCase)
            % Verify analytical result is correct for vertical stripe pattern (nside even).
            tol = 1e-13;

            dmnr = [0, 0, 0, 0; 100, 100, 100, 100; 0, 0, 0, 0; 100, 100, 100, 100].';
            dmexpect = [35, 35, 35, 35; 65, 65, 65, 65; 35, 35, 35, 35; 65, 65, 65, 65].';

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_vcorner_less_than_vneighbor(testCase)
            % Verify analytical result is correct when vcorner is less than vneighbor.
            tol = 1e-13;

            dmnr = [20, 40, 0; 30, 20, 20; 20, 0, 20];
            dmexpect = [20, 35, 5; 27.5, 20, 20; 20, 2.5, 20];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 30, 25);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_vcorner_greater_than_vneighbor(testCase)
            % Verify analytical result is correct when vcorner is greater than vneighbor.
            tol = 1e-13;

            dmnr = [20, 30, 0; 40, 20, 20; 20, 0, 20];
            dmexpect = [20, 27.5, 2.5; 35, 20, 20; 20, 5, 20];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 25, 30);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        function test_analytic_nonuniform_flatmap(testCase)
            % Verify analytical result is correct with a nonuniform flatmap.
            tol = 1e-13;

            dmnr = [100, 0, 100; 0, 100, 0; 100, 0, 100];
            dmflat = [48, 0, 48; 0, 48, 0; 48, 0, 48];
            dmexpect = [99, 1, 99; 1, 99, 1; 99, 1, 99];

            dmout = ConstrainDM.dmsmooth(dmnr, testCase.vmax, 0, 50, 50, dmflat);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end
        
        
        %% Unit tests of ConstrainDM.constrain_dm
        function test_success_basic(testCase)
            % Verify a setting equal to its flatmap, with no ties, completes.
            volts = testCase.flatmap_flat;
            flatmap = testCase.flatmap_flat;
            tie = zeros(testCase.nact);
            ConstrainDM.constrain_dm(volts, flatmap, tie, testCase.vmax, testCase.vlat, testCase.vdiag, testCase.vquant, testCase.maxiter);
        end
        
        function test_constrain_dm_success_flat(testCase)
            % Verify a setting with smoothing and tying completes with a uniform flatmap.
            rng(testCase.rng_seed);
            volts = rand(48, 48);
            flatmap = testCase.flatmap_flat;

            tie = zeros(size(volts));
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4; % single-element group should still work
            tie(21, 21) = -1;

            ConstrainDM.constrain_dm(volts, flatmap, tie, testCase.vmax, testCase.vlat, ...
                testCase.vdiag, testCase.vquant, testCase.maxiter);
        end

        function test_success_nonflat(testCase)
            % Verify a setting with smoothing and tying completes with a non-uniform flatmap.
            rng(testCase.rng_seed);
            volts = rand(48, 48);
            
            tie = zeros(48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4; % single-element group should still work
            tie(21, 21) = -1;
            
            flatmap_nonflat = 50 + 1.0*randn(48, 48);
            flatmap_nonflat(tie == 1) = 50;
            flatmap_nonflat(tie == 2) = 51;
            flatmap_nonflat(tie == 3) = 49;
            flatmap_nonflat(tie == 4) = 50.5;
            flatmap_nonflat(tie == -1) = 0;
            
            ConstrainDM.constrain_dm(volts, flatmap_nonflat, tie, ...
                testCase.vmax, testCase.vlat, testCase.vdiag, testCase.vquant, testCase.maxiter);
        end
        

        function test_nominal(testCase)
            % Verify method returns same array when input is uniform (except dead).
            dmuniform = ones(testCase.nact);
            
            tie = zeros(48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4; % single-element group should still work
            tie(21, 21) = -1;
            
            dmout = ConstrainDM.constrain_dm(dmuniform, testCase.flatmap_flat, tie, ...
                testCase.vmax, testCase.vlat, testCase.vdiag, testCase.vquant, testCase.maxiter);
            testCase.verifyTrue(all(dmuniform(tie ~= -1) == dmout(tie ~= -1)));
        end

        function test_only_smooth(testCase)
            % Analytically test dmsmooth by itself (no ties).
            tol = 1e-13;

            % Made for 30V neighbor rules
            dmnr = [40, 0, 20; 20, 20, 0; 20, 40, 20];
            dmexpect = [35, 5, 20; 20, 20, 5; 20, 35, 20];

            flatmap = zeros(3, 3);
            tie = zeros(3, 3); % no ties

            % Stub out vquant so we can do an analytical test (no margin)
            vmax = testCase.vmax;
            vlat = 30;
            vdiag = 30;
            vquant = 0;
            maxiter = testCase.maxiter;
            dmout = ConstrainDM.constrain_dm(dmnr, flatmap, tie, vmax, vlat, vdiag, vquant, maxiter);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end

        function test_smooth_and_ties(testCase)
            % Analytically test dmsmooth used alongside tie_with_matrix.
            tol = 1e-13;

            % Made for 30V neighbor rules
            dmnr = [40, 0, 20; 20, 20, 0; 20, 40, 20];
            dmexpect = [27.5, 5, 27.5; 20, 20, 5; 20, 35, 20];

            flatmap = zeros(3, 3);
            tie = zeros(3, 3);
            tie(1, 1) = 1;
            tie(1, 3) = 1;

            % Stub out vquant so we can do an analytical test (no margin)
            vmax = testCase.vmax;
            vlat = 30;
            vdiag = 30;
            vquant = 0;
            maxiter = testCase.maxiter;
            dmout = ConstrainDM.constrain_dm(dmnr, flatmap, tie, vmax, vlat, vdiag, vquant, maxiter);

            testCase.verifyTrue(max(max(abs(dmout - dmexpect))) < tol)
        end

        function test_tie_smooth_converge(testCase)
            % Verify that tie and smooth converge.
            checkerboard = reshape(0:49*49-1, [49, 49]);
            checkerboard = mod(checkerboard, 2);
            checkerboard = checkerboard(1:end-1, 1:end-1);
            checkerboard = checkerboard * testCase.vmax;

            tie = zeros(48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4; % single-element group should still work
            tie(21, 21) = -1;

            dmout = ConstrainDM.constrain_dm(checkerboard, testCase.flatmap_flat, tie, ...
                testCase.vmax, testCase.vlat, testCase.vdiag, testCase.vquant, testCase.maxiter);

            % Check still tied (calling tie_with_matrix should change nothing)
            tie_check = ConstrainDM.tie_with_matrix(dmout, tie);
          testCase.verifyTrue(isequal(dmout, tie_check))
            % testCase.verifyTrue(max(max(abs(dmout - tie_check))) < 10*eps)

            % Check still smoothed (calling dmsmooth should change nothing)
            smooth_check = ConstrainDM.dmsmooth(dmout, testCase.vmax, testCase.vquant, ...
                testCase.vlat, testCase.vdiag, testCase.flatmap_flat);
            testCase.verifyTrue(isequal(dmout, smooth_check))

        end
        
        % @patch('howfsc.util.constrain_dm.tie_with_matrix')
        % function test_maxiter(testCase, mock_tie_with_matrix):
        %     """Verify that method will not repeat indefinitely"""
        %     mock_tie_with_matrix.return_value = testCase.checkerboard % not smooth
        %     with testCase.verifyRaises(ConstrainDMException):
        %         constrain_dm(volts=testCase.checkerboard,
        %                      flatmap=testCase.flatmap_flat,
        %                      tie=np.zeros_like(testCase.flatmap_flat),
        %                      maxiter=testCase.maxiter)
        % 
        % @patch('howfsc.util.constrain_dm.tie_with_matrix')
        % function test_maxiter_default(testCase, mock_tie_with_matrix):
        %     """Verify that method will not repeat indefinitely with default"""
        %     mock_tie_with_matrix.return_value = testCase.checkerboard % not smooth
        %     with testCase.verifyRaises(ConstrainDMException):
        %         constrain_dm(volts=testCase.checkerboard,
        %                      flatmap=testCase.flatmap_flat,
        %                      tie=np.zeros_like(testCase.flatmap_flat))      
        
        
        %% Unit tests of ConstrainDM.tie_with_matrix
        
        function test_tie_with_matrix_success(testCase)
            % Basic call with good inputs completes successfully.
            rng(2021)
            dmin = rand(48, 48);
            dmin_untied = dmin;
            
            % tie everything in advance, dmin_untied should produce this
            dmin(1:3, 1:3) = mean(mean(dmin(1:3, 1:3)));
            dmin(6:11, 6:11) = mean(mean(dmin(6:11, 6:11)));
            dmin(7:9, 1:3) = mean(mean(dmin(7:9, 1:3)));
            dmin(21, 21) = 0; % dead
            
            tie = zeros(48, 48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4;
            tie(21, 21) = -1;
            
            ConstrainDM.tie_with_matrix(dmin_untied, tie);
        end
        
        function test_tie_with_matrix_output_matches(testCase)
            % Verify known input matches known output
            rng(2021)
            dmin = rand(48, 48);
            dm_untied = dmin;
            
            % tie everything in advance, dmin_untied should produce this
            dm_tied = dmin;
            dm_tied(1:3, 1:3) = mean(mean(dmin(1:3, 1:3)));
            dm_tied(6:11, 6:11) = mean(mean(dmin(6:11, 6:11)));
            dm_tied(7:9, 1:3) = mean(mean(dmin(7:9, 1:3)));
            dm_tied(21, 21) = 0; % dead
            
            tie = zeros(48, 48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4;
            tie(21, 21) = -1;
            
            out = ConstrainDM.tie_with_matrix(dm_untied, tie);
            testCase.verifyTrue(isequal(out, dm_tied))
        end
        
        function test_tie_with_matrix_idempotence(testCase)
            % Verify that a tied matrix produces itself.
            rng(2021)
            dmin = rand(48, 48);
            dm_untied = dmin;
            
            % tie everything in advance, dmin_untied should produce this
            dm_tied = dmin;
            dm_tied(1:3, 1:3) = mean(mean(dmin(1:3, 1:3)));
            dm_tied(6:11, 6:11) = mean(mean(dmin(6:11, 6:11)));
            dm_tied(7:9, 1:3) = mean(mean(dmin(7:9, 1:3)));
            dm_tied(21, 21) = 0; % dead
            
            tie = zeros(48, 48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4;
            tie(21, 21) = -1;
            
            out = ConstrainDM.tie_with_matrix(dm_tied, tie);            
            testCase.verifyTrue(max(max(abs(out - dm_tied))) < eps)
        end
        
        function test_actually_tied_as_expected(testCase)
            % Verify that the tied regions are all identical and at the mean of the
            % input
            % 
            % Except dead regions are 0, and untied regions are unchanged
            
            rng(2021)
            dmin = rand(48, 48);
            dm_untied = dmin;
            
            % tie everything in advance, dmin_untied should produce this
            dm_tied = dmin;
            dm_tied(1:3, 1:3) = mean(mean(dmin(1:3, 1:3)));
            dm_tied(6:11, 6:11) = mean(mean(dmin(6:11, 6:11)));
            dm_tied(7:9, 1:3) = mean(mean(dmin(7:9, 1:3)));
            dm_tied(21, 21) = 0; % dead
            
            tie = zeros(48, 48);
            tie(1:3, 1:3) = 1;
            tie(6:11, 6:11) = 2;
            tie(7:9, 1:3) = 3;
            tie(19, 19) = 4;
            tie(21, 21) = -1;
            
            out = ConstrainDM.tie_with_matrix(dm_untied, tie);

            % get a list of tie elements
            tienums = unique(tie).';
            for t = tienums
                if t == 0 % not tied
                    % should not change
                    testCase.verifyTrue(all(dm_untied(tie == t) == out(tie == t)))
                    
                elseif t == -1 % dead
                    % should be zeroed
                    testCase.verifyTrue(all(out(tie == t) == 0))

                else
                    % should be mean of input range
                    meanVal = mean(dm_untied(tie == t));
                    testCase.verifyTrue(all(out(tie == t) == meanVal))
                end
            end
        end
        
        
    end
end
