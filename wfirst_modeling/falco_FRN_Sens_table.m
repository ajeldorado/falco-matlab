% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to compute the table in Sensitivities.csv, which is used to  
% compute the flux ratio noise (FRN) for the WFIRST CGI.
%
%
%  The 21 sensitivities--in order--in Sensitivities.csv:
%--Rows 1 to 10: Z2 to Z11 sensitivities to 1nm RMS of Zernike phase aberrations at entrance pupil.
%--Rows 11 to 17: Gain Z5 to Z11 sensitivities. (errors caused by DM
%  acutator gain uncertainty when correcting Z5-Z11 with DMs)
%--Row 18: Pupil X shear
%--Row 19: Pupil Y shear
%--Row 20: DM Settling
%--Row 21: DM Thermal
%
%
% REVISION HISTORY:
% - Created by A.J. Riggs on 2019-05-15.
% -------------------------------------------------------------------------

function tableSens = falco_FRN_Sens_table(mp)
    
%--First 3 columns (the easy, bookkeeping parts of the table)
Nmode = 21; %--Number of sensitivity types that the FRN calculator uses
Nann = size(mp.eval.Rsens,1); %--Number of annuli
matSens = zeros(Nann*Nmode,4); %--Initialize
matSens(:,1) = 0:(Nann*Nmode-1); %--overall index
matSens(:,2) = repmat((0:(Nmode-1)).',[Nann,1]); %--sensmode
for ii=1:Nann;  matSens((ii-1)*Nmode+1:ii*Nmode,3) = ii-1;  end %--annzone

%% Compute sensitivities to 1 micron of X- and Y- Pupil Shear

mp.full.pupilShearVal = 1e-6; %--Amount of shear in each direction [meters]

shearSens = falco_get_pupil_shear_sensitivities(mp); %--dimensions [2 x Nann]
for ii=1:Nann;  matSens((ii-1)*Nmode+18:(ii-1)*Nmode+19,4) = 1e9*shearSens(:,ii);  end %--Multiply by 1e9 for the FRN calculator. Re-organize into column 4 of Sensitivities table

%% Compute sensitivities to 1nm RMS of Zernikes Z2 to Z11

mp.full.ZrmsVal = 1e-9; %--RMS values for each Zernike specified in vector indsZnoll [meters] 
mp.eval.indsZnoll = 2:11; %--Use tip/tilt through spherical modes

Zsens = falco_get_Zernike_sensitivities(mp); %--dimensions of [Nzern x Nann]
for ii=1:Nann;  matSens((ii-1)*Nmode+1:(ii-1)*Nmode+10,4) = 1e9*Zsens(:,ii);  end %--Multiply by 1e9 for the FRN calculator. Re-organize into column 4 of Sensitivities table

%% Compute sensitivities to sigma=10% DM gain error when applying 1nm RMS of Zernikes Z2 to Z11

% %--DO NOT USE YET. NOT GIVING EXPECTED RESULTS
% dE2mat = falco_get_DMgainZ5to11_sensitivities(mp); %--dimensions of [Nzern x Nann]
% for ii=1:Nann;  matSens((ii-1)*Nmode+11:(ii-1)*Nmode+17,4) = 1e9*dE2mat(:,ii);  end %--Multiply by 1e9 for the FRN calculator. Re-organize into column 4 of Sensitivities table

%% Compute sensitivities to DM thermal
% This is the squared change in E-field from the DM actuators moving from a 1 mK drift. 
% We use a thermal response value of 2.6% of strain (strain = actuator heights) per Kelvin. 
% This is the same for all actuators because we do not know otherwise. 

dE2vec = falco_get_DMthermal_sensitivities(mp); %--dimensions of [Nzern x Nann]
for ii=1:Nann;  matSens((ii-1)*Nmode+21,4) = 1e9*dE2vec(ii);  end %--Multiply by 1e9 for the FRN calculator. Re-organize into column 4 of Sensitivities table


%% Make into a table for printing as a CSV file

tableSens = table(matSens(:,1),matSens(:,2),matSens(:,3),matSens(:,4));
tableSens.Properties.VariableNames = {'index','sensmode','annzone','sens'};


end %--END OF FUNCTION