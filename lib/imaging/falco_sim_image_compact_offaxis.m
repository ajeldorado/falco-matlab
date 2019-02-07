% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to return the perfect-knowledge E-field and summed intensity for
% the compact model.
%
% REVISION HISTORY:
%--Modified on 2018-08-27 by A.J. Riggs to include the LOWFS parts from
%    Erkin's code.
%--Created on 2018-01-24 by A.J. Riggs.

function Iout = falco_sim_image_compact_offaxis(mp,x_offset,y_offset,varargin)
    
flagEval  = false;    % flag to use a different (usually higher) resolution at final focal plane for evaluation
modvar.whichSource = 'offaxis';    
modvar.x_offset = x_offset; % mp.thput_eval_x;
modvar.y_offset = y_offset; % mp.thput_eval_y;
  
  icav = 0;             % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'eval'}
            flagEval  = true;   % flag to use a different (usually higher) resolution at final focal plane for evaluation 
%       case {'nact', 'n_act_across_pupil'}
%         icav = icav + 1;
%         nAct = varargin{icav};  % number of actuators across pupil
      otherwise
        error('falco_sim_image_compact: Unknown keyword: %s\n', varargin{icav});
    end
  end


  
    if(flagEval)
        Iout = zeros(mp.F4.eval.Neta, mp.F4.eval.Nxi);
    else
        Iout = zeros(mp.F4.Neta, mp.F4.Nxi);            
    end
    
    for si=1:mp.Nsbp
        modvar.sbpIndex = si; 
        modvar.zernIndex = 1;
        %modvar.ttIndex = mp.Wttlam_ti(im);
        modvar.wpsbpIndex = mp.wi_ref;
        % modvar.whichSource = 'star';

        if(flagEval)
            E2D = model_compact(mp, modvar,'eval');
        else
            E2D = model_compact(mp, modvar);
        end
        
        Iout = Iout + (abs(E2D).^2)*mp.jac.weightMat(si,1);
    end


end %--END OF FUNCTION
