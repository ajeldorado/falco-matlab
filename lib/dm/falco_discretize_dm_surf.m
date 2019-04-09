function dm = falco_discretize_dm_surf(dm,flagMethod)

% Calculate surface heights to maximum precision
h_cont = dm.VtoH.*dm.V;

% Calculate number of (fractional) steps required to represent the surface
nsteps = h_cont./dm.HminStep;

% Discretize by removing fractional step
switch flagMethod
    case {'round'}
        nsteps = round(nsteps);
    case {'floor'}
        nsteps = floor(nsteps);
    case {'ceil'}
        nsteps = ceil(nsteps);
    case {'fix'}
        nsteps = fix(nsteps);
    otherwise
        error([mfilename, ': method for rounding must be one of - round, floor, ceil, or fix'])
end

% Calculate discretized surface
dm.Vquantized = nsteps.*dm.HminStep./dm.VtoH;

end %--END OF FUNCTION