% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_load_dm_defaults(mp)
%
% Function to store pupil mask configuration data for dms.
%
% 
% REVISION HISTORY:
% ----------------
% Modified from falco_config_load_pupil_defaults on 2019-02-21 by G. Ruane
% Created on 2018-05-29 by A.J. Riggs.




function mp = falco_config_load_dm_defaults(mp)

    if(isfield(mp,'dm_ind')==false); mp.dm_ind = [1 2]; end% vector of which DMs to use for control.
    if(isfield(mp,'dm_weights')==false); mp.dm_weights = ones(9,1);  end % vector of relative weighting of DMs for EFC

    %--DM1 parameters
    if(isfield(mp.dm1,'Nact')==false); mp.dm1.Nact = 48; end  % number of actuators across DM1
    if(isfield(mp.dm1,'VtoH')==false); mp.dm1.VtoH = 1*1e-9*ones(mp.dm1.Nact); end  % Gains: volts to meters in surface height;
    if(isfield(mp.dm1,'xtilt')==false); mp.dm1.xtilt = 0; end 
    if(isfield(mp.dm1,'ytilt')==false); mp.dm1.ytilt = 0; end 
    if(isfield(mp.dm1,'zrot')==false); mp.dm1.zrot = 0; end  %--clocking angle (degrees)
    if(isfield(mp.dm1,'xc')==false); mp.dm1.xc = (mp.dm1.Nact/2 - 1/2); end  % x-center of DM in actuator widths
    if(isfield(mp.dm1,'yc')==false); mp.dm1.yc = (mp.dm1.Nact/2 - 1/2); end  % x-center of DM in actuator widths
    if(isfield(mp.dm1,'edgeBuffer')==false); mp.dm1.edgeBuffer = 1; end  % Radius (in actuator spacings) outside of pupil to compute influence functions for.

    %--DM2 parameters
    if(isfield(mp.dm2,'Nact')==false); mp.dm2.Nact = 48; end  % number of actuators across DM1
    if(isfield(mp.dm2,'VtoH')==false); mp.dm2.VtoH = 1*1e-9*ones(mp.dm2.Nact); end  % Gains: volts to meters in surface height;
    if(isfield(mp.dm2,'xtilt')==false); mp.dm2.xtilt = 0; end 
    if(isfield(mp.dm2,'ytilt')==false); mp.dm2.ytilt = 0; end 
    if(isfield(mp.dm2,'zrot')==false); mp.dm2.zrot = 0; end  %--clocking angle (degrees)
    if(isfield(mp.dm2,'xc')==false); mp.dm2.xc = mp.dm2.Nact/2 - 1/2; end  % x-center of DM in actuator widths
    if(isfield(mp.dm2,'yc')==false); mp.dm2.yc = mp.dm2.Nact/2 - 1/2; end  % x-center of DM in actuator widths
    if(isfield(mp.dm2,'edgeBuffer')==false); mp.dm2.edgeBuffer = 1; end  % Radius (in actuator spacings) outside of pupil to compute influence functions for.

    %--DM Actuator characteristics
    if(isfield(mp.dm1,'inf_fn')==false); mp.dm1.inf_fn = 'influence_dm5v2.fits'; end  % file name of influence function
    if(isfield(mp.dm1,'inf_sign')==false); mp.dm1.inf_sign = '+'; end  % sign of the influence function [+ or -]

    if(isfield(mp.dm2,'inf_fn')==false); mp.dm2.inf_fn = 'influence_dm5v2.fits'; end  % file name of influence function
    if(isfield(mp.dm2,'inf_sign')==false); mp.dm2.inf_sign = '+'; end  % sign of the influence function [+ or -]

    if(isfield(mp.dm1,'dm_spacing')==false);  mp.dm1.dm_spacing = 1e-3;  end% inter-actuator pitch THAT YOU DEFINE. Does not have to be the true value [meters]
    if(isfield(mp.dm2,'dm_spacing')==false);  mp.dm2.dm_spacing = 1e-3;  end% inter-actuator pitch THAT YOU DEFINE. Does not have to be the true value [meters]


end %--END OF FUNCTION



