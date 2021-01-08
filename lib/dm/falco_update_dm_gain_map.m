% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Update the DM gain map based on the total voltage command.
% Accounts for nonlinear displacement vs voltage.
%
% INPUTS
% ------
% dm : structure of DM parameters. Is either mp.dm1 or mp.dm2
%
% OUTPUTS
% -------
% dm : structure of DM parameters. Is either mp.dm1 or mp.dm2

function dm = falco_update_dm_gain_map(dm)

switch lower(dm.fitType)
    
    case{'linear', 'poly1'}
        % No change to dm.VtoH
        
    case{'quadratic', 'poly2'}
        if ~isfield(dm, 'p1') || ~isfield(dm, 'p2') || ~isfield(dm, 'p3')
            error(sprintf("The fields p1, p2, and p3 must exist when dm.fitType == 'quadratic'."+...
                "\nThose fields satisfy the formula:\n"+...
                "height = p1*V*V + p2*V + p3"))
        end
        Vtotal = dm.V + dm.biasMap;
        dm.VtoH = 2*dm.p1.*Vtotal + dm.p2;
    
    case{'fourier2'}
        if ~isfield(dm, 'a0') || ~isfield(dm, 'a1') || ~isfield(dm, 'a2') || ...
                ~isfield(dm, 'b1') || ~isfield(dm, 'b2') || ~isfield(dm, 'w')
            error(sprintf("The fields a0, a1, a2, b1, b2, and w must exist when dm.fitType == 'fourier2'."+...
                "\nThose fields satisfy the formula:\n"+...
                "height = a0 + a1*cos(V*w) + b1*sin(V*w) + a2*cos(2*V*w) + b2*sin(2*V*w)"))
        end        
        Vtotal = dm.V + dm.biasMap;
        dm.VtoH = dm.w.*(-dm.a1.*sin(Vtotal.*dm.w) + dm.b1.*cos(Vtotal.*dm.w) + ...
               -2*dm.a2.*sin(2*Vtotal.*dm.w) + 2*dm.b2.*cos(2*Vtotal.*dm.w)); 
           
    otherwise
        error('Value of dm.fitType not recognized.')
end

end %--END OF FUNCTION