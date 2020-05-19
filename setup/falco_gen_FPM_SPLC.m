% Copyright 2018,2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% REVISION HISTORY:
% --------------
% Modified on 2019-05-09 by A.J. Riggs to be more generic a function by
% using the flags mp.full.flagGenFPM and mp.compact.flagGenFPM.
% Created on 2018-10-01 by A.J. Riggs by extracting material from
% falco_init_ws.m.
% ---------------

function mp = falco_gen_FPM_SPLC(mp)

mp.F3.dummy = 1; %--Make sure mp.F3 exists
if(isfield(mp.F3,'ang')==false);  mp.F3.ang = 180;  end  % [degrees]

if(mp.full.flagGenFPM)
    %--Generate the FPM amplitude for the full model
    inputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
    inputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
    inputs.ang = mp.F3.ang;  % [degrees]
    inputs.centering = mp.centering;
    inputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
    mp.F3.full.mask.amp = falco_gen_bowtie_FPM(inputs);
end

if(mp.compact.flagGenFPM)
    %--Generate the FPM amplitude for the compact model
    inputs.rhoInner = mp.F3.Rin; % radius of inner FPM amplitude spot (in lambda_c/D)
    inputs.rhoOuter = mp.F3.Rout; % radius of outer opaque FPM ring (in lambda_c/D)
    inputs.ang = mp.F3.ang;  % [degrees]
    inputs.centering = mp.centering;
    inputs.pixresFPM = mp.F3.compact.res;
    mp.F3.compact.mask.amp = falco_gen_bowtie_FPM(inputs);        
end

if(mp.full.flagPROPER==false)
    mp.F3.full.Nxi = size(mp.F3.full.mask.amp,2);
    mp.F3.full.Neta= size(mp.F3.full.mask.amp,1);
end

mp.F3.compact.Nxi = size(mp.F3.compact.mask.amp,2);
mp.F3.compact.Neta= size(mp.F3.compact.mask.amp,1);
        
end %--END OF FUNCTION