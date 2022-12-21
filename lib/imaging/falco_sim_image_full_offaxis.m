% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to return the band-averaged intensity for the full model.

function [Iout, varargout] = falco_sim_image_full_offaxis(mp, x_offset, y_offset)
    
modvar = ModelVariables;
modvar.whichSource = 'offaxis';
modvar.x_offset = x_offset; % mp.thput_eval_x;
modvar.y_offset = y_offset; % mp.thput_eval_y;
  
if(mp.flagFiber);  Ifiber = zeros(mp.F5.Neta, mp.F5.Nxi);  end %--Initialize
Iout = 0; %--Initialize
for si=1:mp.Nsbp
    modvar.sbpIndex = si; 
    modvar.zernIndex = 1;
    modvar.wpsbpIndex = mp.wi_ref;

    if(mp.flagFiber)
        [E2D, Efiber] = model_full(mp, modvar);
    else
        E2D = model_full(mp, modvar);
    end
        
    Iout = Iout + (abs(E2D).^2)*mp.jac.weightMat(si,1);
    if(mp.flagFiber)
        Ifiber = Ifiber + (abs(Efiber).^2)*mp.jac.weightMat(si,1);
        varargout{1} = Ifiber;
    end
end

end %--END OF FUNCTION