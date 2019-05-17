% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to compute the annular zone table for the WFIRST CGI flux ratio 
% noise (FRN) calculation).
%
% REVISION HISTORY:
% - Created by A.J. Riggs on 2019-05-16.
% -------------------------------------------------------------------------

function tableAnn = falco_FRN_AnnularZone_table(mp)
    
    %--Annular Zone Table
    Nann = size(mp.eval.Rsens,1);
    matAnn = zeros(2*Nann,3);
    matAnn(:,1) = [(0:(Nann-1)).';(0:(Nann-1)).']; %--First column
    matAnn(Nann+1:end,2) = 1; %--Second column
    matAnn(1:Nann,3) = mp.eval.Rsens(:,1); %--Third column, 1st half
    matAnn(Nann+1:end,3) = mp.eval.Rsens(:,2); %--Third column, 2nd half
    
    %--Make into a table for printing as a CSV file
    tableAnn = table(matAnn(:,1),matAnn(:,2),matAnn(:,3));
    tableAnn.Properties.VariableNames = {'index_zone','min_vs_max','value_lamOverD'};
    
end %--END OF FUNCTION