% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Generate FPMs. For hybrid FPMs only ('HLC' or 'EHLC'), which use DM8 and DM9.

function [mp] = falco_gen_fpm(mp)
    
    %--Stash DM8 and DM9 starting commands if they are given in the main script
    if isfield(mp, 'dm8')
        if(isfield(mp.dm8, 'V')) 
            mp.DM8V0 = mp.dm8.V; 
        end
        if(isfield(mp.dm9, 'V'))
            mp.DM9V0 = mp.dm9.V; 
        end
    end
    
    switch lower(mp.layout)
        case 'fourier'
            switch upper(mp.coro)
                
                case{'HLC'}
                        
                    if ~isfield(mp.compact, 'FPMcube')
                        
                        switch upper(mp.dm9.inf0name)
                            case '3FOLDZERN'
                                mp = falco_setup_FPM_HLC_3foldZern(mp);
                            case{'COS', 'COSINE'}
                                mp = falco_setup_FPM_HLC_cosine(mp);
                            otherwise
                                mp = falco_setup_FPM_HLC(mp);
                        end
                        mp = falco_gen_FPM_HLC(mp);

                    end
                    
                case{'EHLC'}
                    mp = falco_setup_FPM_EHLC(mp);
                    mp = falco_gen_FPM_EHLC(mp);
            end

            %%--Pre-compute the complex transmission of the allowed Ni+PMGI FPMs.
            switch upper(mp.coro)
                case{'HLC', 'EHLC'}
                    if ~isfield(mp.compact, 'FPMcube')
                        [mp.complexTransCompact,mp.complexTransFull] = falco_gen_complex_trans_table(mp);
                    end
            end
    end

end %--END OF FUNCTION
