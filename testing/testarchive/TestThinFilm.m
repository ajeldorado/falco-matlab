%---------------------------------------------------------------------------
% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%---------------------------------------------------------------------------
%% Test falco_thin_film_material_def.m
%
% We define some tests for falco_thin_film_material_def.m to test responses to
% different input parameters. 
classdef TestThinFilm < matlab.unittest.TestCase    
%% Setup and Teardown Methods
%
%  Add and remove path to utils functions to be tested.
%
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath('../../lib/thinfilm'));
        end
    end
    methods (TestClassTeardown)
        function removePath(testCase)
            rmpath(genpath('../../lib/thinfilm'));
        end
    end
     
%% Tests
                                      
    methods (Test)    
        function testTransmissionPMGI(testCase)
            % From PMGI_transmission_only.xlsx
            lam = 600e-9;
            d0 = 4*lam;
            aoi = 0;
            t_Ti_base = 0;
            t_Ni_vec = 0;
            t_PMGI_vec = 100e-9;
            pol = 2;
            [tCoef, ~] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
            T_FALCO = abs(tCoef)^2; % Value from Bala: 0.9431006949; Value from FALCO: 0.94313153
            T_Macleod = 0.9431006949;

            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(T_FALCO, IsEqualTo(T_Macleod, 'Within', RelativeTolerance(0.0001)))
        end
        
        function testTransmissionNickel(testCase)
            % From PMGIon95nmNi_aoi10deg_T_sPol.xlsx
            lam = 400e-9;
            d0 = 4*lam;
            aoi = 10;
            t_Ti_base = 0;
            t_Ni_vec = 95e-9;
            t_PMGI_vec = 0;
            pol = 0;
            [tCoef, ~] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
            T_FALCO = abs(tCoef)^2; % Value from Bala: 0.00087848574  Value from FALCO: 0.000878466587
            T_Macleod = 0.00087848574;

            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(T_FALCO, IsEqualTo(T_Macleod, 'Within', RelativeTolerance(0.0001)))         
        end
        
        function testTransmissionPMGIonNickel1of2(testCase)
            % From PMGIon95nmNi_aoi10deg_T_sPol.xlsx
            lam = 450e-9;
            d0 = 4*lam;
            aoi = 10;
            t_Ti_base = 0;
            t_Ni_vec = 95e-9;
            t_PMGI_vec = 30e-9;
            pol = 1;
            [tCoef, ~] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
            T_FALCO = abs(tCoef)^2; % Value from Bala: 0.00118382732,   Value from FALCO: 0.00118379
            T_Macleod = 0.00118382732;

            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(T_FALCO, IsEqualTo(T_Macleod, 'Within', RelativeTolerance(0.0001)))         
        end
        
        function testTransmissionPMGIonNickel2of2(testCase)
            % From PMGIon95nmPMGI_aoi10deg_T_pPol.xlsx
            lam = 550e-9;
            d0 = 4*lam;
            aoi = 10;
            t_Ti_base = 0;
            t_Ni_vec = 95e-9;
            t_PMGI_vec = 600e-9;
            pol = 1;
            [tCoef, ~] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol);
            T_FALCO = abs(tCoef)^2; % Value from Bala: 0.00121675706  Value from FALCO: 0.001216750339
            T_Macleod = 0.00121675706;
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            testCase.verifyThat(T_FALCO, IsEqualTo(T_Macleod, 'Within', RelativeTolerance(0.0001)))         
        end

    end    
end