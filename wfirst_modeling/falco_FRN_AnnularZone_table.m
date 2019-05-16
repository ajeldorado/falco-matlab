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
    tableAnn = zeros(2*Nann,3);
    tableAnn(:,1) = [(0:(Nann-1)).';(0:(Nann-1)).']; %--First column
    tableAnn(Nann+1:end,2) = 1; %--Second column
    tableAnn(1:Nann,3) = mp.eval.Rsens(:,1); %--Third column, 1st half
    tableAnn(Nann+1:end,3) = mp.eval.Rsens(:,2); %--Third column, 2nd half
    
end %--END OF FUNCTION