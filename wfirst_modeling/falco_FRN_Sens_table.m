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
% tableSens

%% Compute sensitivities to 1nm RMS of Zernikes Z2 to Z11
%     dE2mat = zeros;

mp.full.ZrmsVal = 1e-9; %--RMS values for each Zernike specified in vector indsZnoll [meters] 
mp.eval.indsZnoll = 2:11; %--Use tip/tilt through spherical modes

Zsens = falco_get_Zernike_sensitivities(mp); %--dimensions of [Nzern,Nann]
for ii=1:Nann;  matSens((ii-1)*Nmode+1:(ii-1)*Nmode+10,4) = 1e9*Zsens(:,ii);  end %--Multiply by 1e9 to go from m to nm. Re-organize into column 4 of Sensitivities table

%% Compute sensitivities to 1 micron of X- and Y- Pupil Shear


%% Make into a table for printing as a CSV file

tableSens = table(matSens(:,1),matSens(:,2),matSens(:,3),matSens(:,4));
tableSens.Properties.VariableNames = {'index','sensmode','annzone','sens'};


end %--END OF FUNCTION