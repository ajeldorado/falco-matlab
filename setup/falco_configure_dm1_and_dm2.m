% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

function mp = falco_configure_dm1_and_dm2(mp)

if( any(mp.dm_ind==1) &&  any(mp.dm_ind==2) ) 
    disp(['DM 1-to-2 Fresnel number (using radius) = ',num2str((mp.P2.D/2)^2/(mp.d_dm1_dm2*mp.lambda0))]);
end

%% DM1
% Read the influence function header data from the FITS file
info = fitsinfo(mp.dm1.inf_fn);
[~, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1)); % Get index in cell array
dx1 = info.PrimaryData.Keywords{idef, 2}; % pixel width of the influence function IN THE FILE [meters];
[~, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
pitch1 = info.PrimaryData.Keywords{idef, 2}; % actuator spacing x (m)

mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
mp.dm1.dx_inf0 = mp.dm1.dm_spacing*(dx1/pitch1);

switch lower(mp.dm1.inf_sign(1))
    case{'-','n','m'}
        mp.dm1.inf0 = -1*mp.dm1.inf0;
    otherwise
        %--Leave coefficient as +1
end   

%--Create influence function datacubes for each DM
mp.dm1.centering = mp.centering;
mp.dm1.compact = mp.dm1;
mp.dm1 = falco_gen_dm_poke_cube(mp.dm1, mp, mp.P2.full.dx,'NOCUBE');
if( any(mp.dm_ind==1) )
    mp.dm1.compact = falco_gen_dm_poke_cube(mp.dm1.compact, mp, mp.P2.compact.dx);
else
    mp.dm1.compact = falco_gen_dm_poke_cube(mp.dm1.compact, mp, mp.P2.compact.dx,'NOCUBE');
end
    
if(isfield(mp.dm1,'V')==false); mp.dm1.V = zeros(mp.dm1.Nact,mp.dm1.Nact); end %--Initial DM voltages

%% DM2
% Read the influence function header data from the FITS file
info = fitsinfo(mp.dm2.inf_fn);
[~, idef] = ismember('P2PDX_M', info.PrimaryData.Keywords(:, 1)); % Get index in cell array
dx2 = info.PrimaryData.Keywords{idef, 2}; % pixel width of the influence function IN THE FILE [meters];
[~, idef] = ismember('C2CDX_M', info.PrimaryData.Keywords(:, 1));
pitch2 = info.PrimaryData.Keywords{idef, 2}; % actuator spacing x (m)

mp.dm2.inf0 = fitsread(mp.dm2.inf_fn);
mp.dm2.dx_inf0 = mp.dm2.dm_spacing*(dx2/pitch2);

switch lower(mp.dm2.inf_sign(1))
    case{'-','n','m'}
        mp.dm2.inf0 = -1*mp.dm2.inf0;
    otherwise
        %--Leave coefficient as +1
end

mp.dm2.centering = mp.centering;
mp.dm2.compact = mp.dm2;
mp.dm2.dx = mp.P2.full.dx;
mp.dm2.compact.dx = mp.P2.compact.dx;
mp.dm2 = falco_gen_dm_poke_cube(mp.dm2, mp, mp.P2.full.dx, 'NOCUBE');
if( any(mp.dm_ind==2) ) 
    mp.dm2.compact = falco_gen_dm_poke_cube(mp.dm2.compact, mp, mp.P2.compact.dx);
else
    mp.dm2.compact = falco_gen_dm_poke_cube(mp.dm2.compact, mp, mp.P2.compact.dx,'NOCUBE');
end

if(isfield(mp.dm2,'V')==false); mp.dm2.V = zeros(mp.dm2.Nact,mp.dm2.Nact); end %--Initial DM voltages

end