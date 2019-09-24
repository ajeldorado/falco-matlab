function dm = falco_discretize_dm_surf(dm,flagMethod,varargin)

% optional flag to use testbed instead of simulation  
% if(nargin>0)
%     flag_tb = varargin{1}; % The first variable argument defines which testbed type
% 
%     switch flag_tb
%         case 'tb'
%             HminStep = dm.HminStep_tb;
%         otherwise
%             HminStep = dm.HminStep;
%             disp('Silent error: HminStep_tb not defined. Using HminStep for testbed.');
%     end
% else
    % Use HminStep for dm model 
    HminStep = dm.HminStep;
% end

% Calculate surface heights to maximum precision
h_cont = dm.VtoH.*dm.V;

% Calculate number of (fractional) steps required to represent the surface
nsteps = h_cont./HminStep;

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
dm.Vquantized = nsteps.*HminStep./dm.VtoH;

end %--END OF FUNCTION