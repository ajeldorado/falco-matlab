% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% 
% Convert DM commands (usually voltages) into heights (before they are
% convolved with the influence function to make a surface).
% Accounts for nonlinear displacement vs voltage.
%
% INPUTS
% ------
% dm: structure of DM parameters. Is either mp.dm1 or mp.dm2
%
% OUTPUTS
% -------
% heightMap: 2-D map of DM actuator heights

function heightMap = falco_calc_act_height_from_voltage(dm)

switch lower(dm.fitType)
    
    case{'linear', 'poly1'}
        heightMap = dm.VtoH .* dm.V;
        
    case{'quadratic', 'poly2'}
        if ~isfield(dm, 'p1') || ~isfield(dm, 'p2') || ~isfield(dm, 'p3')
            error(sprintf("The fields p1, p2, and p3 must exist when dm.fitType == 'quadratic'."+...
                "\nThose fields satisfy the formula:\n"+...
                "height = p1*V*V + p2*V + p3"))
        end
        Vbias = dm.biasMap;
        Vtotal = dm.V + Vbias;
        heightMapTotal = dm.p1.*Vtotal.^2 + dm.p2.*Vtotal + dm.p3;
        heightMapBias = dm.p1.*Vbias.^2 + dm.p2.*Vbias + dm.p3;
        heightMap = heightMapTotal - heightMapBias;
        
    case{'fourier2'}
        if ~isfield(dm, 'a0') || ~isfield(dm, 'a1') || ~isfield(dm, 'a2') || ...
                ~isfield(dm, 'b1') || ~isfield(dm, 'b2') || ~isfield(dm, 'w')
            error(sprintf("The fields a0, a1, a2, b1, b2, and w must exist when dm.fitType == 'fourier2'."+...
                "\nThose fields satisfy the formula:\n"+...
                "height = a0 + a1*cos(V*w) + b1*sin(V*w) + a2*cos(2*V*w) + b2*sin(2*V*w)"))
        end  
        Vbias = dm.biasMap;
        Vtotal = dm.V + Vbias;
        heightMapTotal = dm.a0 + dm.a1.*cos(Vtotal.*dm.w) + dm.b1.*sin(Vtotal.*dm.w) + dm.a2.*cos(2*Vtotal.*dm.w) + dm.b2.*sin(2*Vtotal.*dm.w);
        heightMapBias = dm.a0 + dm.a1.*cos(Vbias.*dm.w) + dm.b1.*sin(Vbias.*dm.w) + dm.a2.*cos(2*Vbias.*dm.w) + dm.b2.*sin(2*Vbias.*dm.w);
        heightMap = heightMapTotal - heightMapBias;        

    otherwise
        error('Value of dm.fitType not recognized.')
end

end %--END OF FUNCTION
