% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

function mp = falco_init_storage_arrays(mp)
    %--EFC regularization history
    out.log10regHist = zeros(mp.Nitr,1);

    %--Peak-to-Valley DM voltages
    out.dm1.Vpv = zeros(mp.Nitr,1);
    out.dm2.Vpv = zeros(mp.Nitr,1);
    out.dm8.Vpv = zeros(mp.Nitr,1);
    out.dm9.Vpv = zeros(mp.Nitr,1);

    %--Peak-to-Valley DM surfaces
    out.dm1.Spv = zeros(mp.Nitr,1);
    out.dm2.Spv = zeros(mp.Nitr,1);
    out.dm8.Spv = zeros(mp.Nitr,1);
    out.dm9.Spv = zeros(mp.Nitr,1);

    %--RMS DM surfaces
    out.dm1.Srms = zeros(mp.Nitr,1);
    out.dm2.Srms = zeros(mp.Nitr,1);
    out.dm8.Srms = zeros(mp.Nitr,1);
    out.dm9.Srms = zeros(mp.Nitr,1);

    %--Sensitivities Zernike-Mode Perturbations
    Nannuli = size(mp.eval.Rsens,1);
    Nzern = length(mp.eval.indsZnoll);
    out.Zsens = zeros(Nzern,Nannuli,mp.Nitr);

    %--Store the DM commands at each iteration
    if(isfield(mp,'dm1')); if(isfield(mp.dm1,'V'));  out.dm1.Vall = zeros(mp.dm1.Nact,mp.dm1.Nact,mp.Nitr+1);  end; end
    if(isfield(mp,'dm2')); if(isfield(mp.dm2,'V'));  out.dm2.Vall = zeros(mp.dm2.Nact,mp.dm2.Nact,mp.Nitr+1); end; end
    if(isfield(mp,'dm8')); if(isfield(mp.dm8,'V'));  out.dm8.Vall = zeros(mp.dm8.NactTotal,mp.Nitr+1); end; end
    if(isfield(mp,'dm9')); if(isfield(mp.dm9,'V'));  out.dm9.Vall = zeros(mp.dm9.NactTotal,mp.Nitr+1); end; end

end