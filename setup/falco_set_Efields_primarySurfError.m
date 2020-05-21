% Copyright 2019, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

function mp = falco_set_Efields_primarySurfError(mp,surf,varargin)
%mp = falco_set_Efields_primarySurfError(mp,surf,varargin)
%   Converts a surface map into E-field at P1. Fills in the mp.P1.x.E cube
%   with the appropriate phase at each wavelength. Assumes compact model
%   and full model use same sampling. Useful for studying the impact of
%   surface deformations in the telescope. 
%   
%   Inputs: 
%       mp - model parameters structure 
%       surf - primary mirror surface (meters)
%       varargin(1) - 'overwrite' (default), 'overwrite_phz', or 'perturb'
%           'overwrite' - overwrites any previous data in E-field cube.
%                         The result is just the amplitude of 1 and phase 
%                         due to primary mirror OPD.
%           'overwrite_phz' - overwrites previous phase data in E-field cube.
%                         The result is that the amplitude remains the same
%                         and the phase is only due to primary mirror OPD.   
%           'perturb' - Keeps previous E-field information and pertubs it
%                         by the phase due to primary mirror OPD. 
%
%   Outputs: 
%       mp - updated model parameters structure with new E-field cubes.
%
    
    if(nargin>2)
        flag = varargin(1);
    else
        flag = 'overwrite';
    end

    % If the E-field cube doesn't exist, use overwrite mode 
    if(~isfield(mp.P1.compact,'E') || ~isfield(mp.P1.full,'E'))
        flag = 'overwrite';
    end
        
    
    % Compact model has one field per sub-bandpass 
    for sbpIndex = 1:mp.Nsbp
        
        lambda = mp.sbp_centers(sbpIndex);
        phz = 4*pi*surf/lambda;
        if(strcmpi(flag,'overwrite'))%overwrite the whole E-field
            mp.P1.compact.E(:,:,sbpIndex) = exp(1i*phz);
        elseif(strcmpi(flag,'overwrite_phz'))%overwrite the phase
            mp.P1.compact.E(:,:,sbpIndex) = ...
                abs(mp.P1.compact.E(:,:,sbpIndex)).*exp(1i*phz);
        elseif(strcmpi(flag,'perturb'))%perturb the E-field
            mp.P1.compact.E(:,:,sbpIndex) = ...
                mp.P1.compact.E(:,:,sbpIndex).*exp(1i*phz);
        else
            error(['falco_primarySurf2Efields: update method ',flag,' is not defined. Use overwrite, overwrite_phz, or perturb.']);
        end
        
    end

    % Full model can have more than one field per sub-bandpass
    for sbpIndex = 1:mp.Nsbp
        for wpsbpIndex = 1:mp.Nwpsbp
            
            lambda = mp.full.lambdasMat(sbpIndex,wpsbpIndex);
            phz = 4*pi*surf/lambda;
            
            if(strcmpi(flag,'overwrite'))%overwrite the whole E-field
                mp.P1.full.E(:,:,wpsbpIndex,sbpIndex) = exp(1i*phz);
            elseif(strcmpi(flag,'overwrite_phz'))%overwrite the phase
                mp.P1.full.E(:,:,wpsbpIndex,sbpIndex) = ...
                    abs(mp.P1.full.E(:,:,wpsbpIndex,sbpIndex)).*exp(1i*phz);
            elseif(strcmpi(flag,'perturb'))%perturb the E-field
                mp.P1.full.E(:,:,wpsbpIndex,sbpIndex) = ...
                    mp.P1.full.E(:,:,wpsbpIndex,sbpIndex).*exp(1i*phz);
            else
                error(['falco_primarySurf2Efields: update method ',flag,' is not defined. Use overwrite, overwrite_phz, or perturb.']);
            end
            
        end
    end

%     if(mp.flagPlot)
%         figure;imagesc(angle(mp.P1.full.E(:,:,1,1))/2/pi);axis image;colorbar;title('P1 phase at lam1 (waves)');
%     end

end

