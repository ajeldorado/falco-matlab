% Copyright 2018-2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
% Any commercial use must be negotiated with the Office 
% of Technology Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute FPM coordinates, both in meters and lambda/D.

function mp = falco_compute_fpm_coordinates(mp)

switch upper(mp.coro)
    
    case{'VORTEX','VC'}
        % No coordinates needed in vortex models
    
    otherwise
        
        fLamD = mp.fl*mp.lambda0/mp.P2.D;
        
        %% Compact Model
        mp.F3.compact.dummy = 1;
        if ~isfield(mp.F3.compact, 'Nxi')
            mp.F3.compact.Nxi = size(mp.F3.compact.mask, 2);
        end
        if ~isfield(mp.F3.compact, 'Neta')
            mp.F3.compact.Neta= size(mp.F3.compact.mask, 1);
        end
        
        % Resolution in compact model
        mp.F3.compact.dxi = fLamD/mp.F3.compact.res; % [meters/pixel]
        mp.F3.compact.deta = mp.F3.compact.dxi; % [meters/pixel]
        
        %--FPM coordinates in the compact model [meters]
        % horizontal axis, xis
        if strcmpi(mp.centering, 'interpixel') || mod(mp.F3.compact.Nxi, 2) == 1
            mp.F3.compact.xis  = (-(mp.F3.compact.Nxi-1)/2:(mp.F3.compact.Nxi-1)/2)*mp.F3.compact.dxi;
        elseif strcmpi(mp.centering, 'pixel')
            mp.F3.compact.xis  = (-mp.F3.compact.Nxi/2: (mp.F3.compact.Nxi/2-1))*mp.F3.compact.dxi;
        end
        % vertical axis, etas
        if strcmpi(mp.centering, 'interpixel') || mod(mp.F3.compact.Neta, 2) == 1
            mp.F3.compact.etas = (-(mp.F3.compact.Neta-1)/2:(mp.F3.compact.Neta-1)/2).'*mp.F3.compact.deta;
        elseif strcmpi(mp.centering, 'pixel')
            mp.F3.compact.etas = (-mp.F3.compact.Neta/2:(mp.F3.compact.Neta/2-1)).'*mp.F3.compact.deta;
        end
        
        %--Dimensionless FPM Coordinates in compact model
        mp.F3.compact.xisDL = mp.F3.compact.xis / fLamD;
        mp.F3.compact.etasDL = mp.F3.compact.etas / fLamD;

        
        %% Full Model
        switch lower(mp.layout)
            
            case{'roman_phasec_proper','wfirst_phaseb_proper', 'proper'}
                % Coorinates not needed in the PROPER model
            
            otherwise
                
                if mp.full.flagPROPER == false
                    
                    mp.F3.full.dummy = 1;
                    if ~isfield(mp.F3.full, 'Nxi')
                        mp.F3.full.Nxi = size(mp.F3.full.mask, 2);
                    end
                    if ~isfield(mp.F3.full, 'Neta')
                        mp.F3.full.Neta= size(mp.F3.full.mask, 1);
                    end
                
                    %--Resolution
                    mp.F3.full.dxi = fLamD/mp.F3.full.res; % [meters/pixel]
                    mp.F3.full.deta = mp.F3.full.dxi; % [meters/pixel]

                    %--Coordinates (dimensionless [DL]) for the FPMs in the full model
                    if strcmpi(mp.centering,'interpixel') || mod(mp.F3.full.Nxi, 2) == 1
                        mp.F3.full.xisDL  = (-(mp.F3.full.Nxi-1)/2:(mp.F3.full.Nxi-1)/2)/mp.F3.full.res;
                        mp.F3.full.etasDL = (-(mp.F3.full.Neta-1)/2:(mp.F3.full.Neta-1)/2)/mp.F3.full.res;
                    else
                        mp.F3.full.xisDL  = (-mp.F3.full.Nxi/2:(mp.F3.full.Nxi/2-1))/mp.F3.full.res;
                        mp.F3.full.etasDL = (-mp.F3.full.Neta/2:(mp.F3.full.Neta/2-1))/mp.F3.full.res;
                    end
                    
                end
        end
        
end


end %--END OF FUNCTION
