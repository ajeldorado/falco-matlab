% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

function [mp] = falco_gen_FPM(mp)

    switch upper(mp.coro)
        case {'LC','APLC'} %--Occulting spot FPM (can be partially transmissive)
            mp = falco_gen_FPM_LC(mp);
        case{'SPLC','FLC'}
            mp = falco_gen_FPM_SPLC(mp);
        case{'RODDIER'}
            mp = falco_gen_FPM_Roddier(mp);  
    end
    
    %% Hybrid FPMs only:
    %--Stash DM8 and DM9 starting commands if they are given in the main script
    if(isfield(mp,'dm8'))
        if(isfield(mp.dm8,'V')); mp.DM8V0 = mp.dm8.V; end
        if(isfield(mp.dm9,'V')); mp.DM9V0 = mp.dm9.V; end
    end
    switch lower(mp.layout)
        case 'fourier'
            switch upper(mp.coro)
                case{'HLC'}
                    switch mp.dm9.inf0name
                        case '3foldZern'
                            mp = falco_setup_FPM_HLC_3foldZern(mp);
                        otherwise
                            mp = falco_setup_FPM_HLC(mp);
                    end
                    mp = falco_gen_FPM_HLC(mp);
                case{'FOHLC'}
                    mp = falco_setup_FPM_FOHLC(mp);
                    mp = falco_gen_FPM_FOHLC(mp);
                    mp.compact.Nfpm = max([mp.dm8.compact.NdmPad,mp.dm9.compact.NdmPad]); %--Width of the FPM array in the compact model.
                    mp.full.Nfpm = max([mp.dm8.NdmPad,mp.dm9.NdmPad]); %--Width of the FPM array in the full model.
                case{'EHLC'}
                    mp = falco_setup_FPM_EHLC(mp);
                    mp = falco_gen_FPM_EHLC(mp);
            end

            %%--Pre-compute the complex transmission of the allowed Ni+PMGI FPMs.
            switch upper(mp.coro)
                case{'EHLC','HLC'}
                    [mp.complexTransCompact,mp.complexTransFull] = falco_gen_complex_trans_table(mp);
            end
    end

end %--END OF FUNCTION