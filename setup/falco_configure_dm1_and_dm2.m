% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Flesh out the dm1 and dm2 structures
%
% INPUTS
% ------
% mp : structure of model parameters
%
% OUTPUTS
% -------
% mp : structure of model parameters

function mp = falco_configure_dm1_and_dm2(mp)

if any(mp.dm_ind == 1) && any(mp.dm_ind == 2)
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

mp.dm1.compact.dummy = 1;
mdc = mp.dm1.compact;
mp.dm1.compact = mp.dm1; % BE CAREFUL ABOUT OVERWRITING VARIABLES
for fn = fieldnames(mdc)
    mp.dm1.compact.(fn{1}) = mdc.(fn{1});
end

mp.dm1 = falco_gen_dm_poke_cube(mp.dm1, mp, mp.P2.full.dx,'NOCUBE');
if( any(mp.dm_ind==1) )
    mp.dm1.compact = falco_gen_dm_poke_cube(mp.dm1.compact, mp, mp.P2.compact.dx);
else
    mp.dm1.compact = falco_gen_dm_poke_cube(mp.dm1.compact, mp, mp.P2.compact.dx,'NOCUBE');
end
    
if(isfield(mp.dm1,'V')==false); mp.dm1.V = zeros(mp.dm1.Nact,mp.dm1.Nact); end %--Initial DM voltages

%%
switch mp.dm1.basisType
    
    case 'actuator'
        
        mp.dm1.NbasisModes = mp.dm1.NactTotal;
        mp.dm1.basisCube = zeros(mp.dm1.Nact, mp.dm1.Nact, mp.dm1.NbasisModes);
        for iact = 1:mp.dm1.NactTotal
            tempVec = zeros(mp.dm1.NactTotal, 1);
            tempVec(iact) = 1;
            mp.dm1.basisCube(:, :, iact) = reshape(tempVec, [mp.dm1.Nact, mp.dm1.Nact]);
        end
        
        mp.dm2.NbasisModes = mp.dm2.NactTotal;
        mp.dm2.basisCube = zeros(mp.dm2.Nact, mp.dm2.Nact, mp.dm2.NbasisModes);
        for iact = 1:mp.dm2.NactTotal
            tempVec = zeros(mp.dm2.NactTotal, 1);
            tempVec(iact) = 1;
            mp.dm2.basisCube(:, :, iact) = reshape(tempVec, [mp.dm2.Nact, mp.dm2.Nact]);
        end
        
    case 'fourier'
        
        if length(mp.dm1.fourier_basis_xis) ~= length(mp.dm1.fourier_basis_etas)
           error('fourier_basis_xis and fourier_basis_etas must have the same length')
        end
        mp.dm1.NbasisModes = 2 * length(mp.dm1.fourier_basis_xis); % 2 for sin and cos at each freq
        mp.dm1.basisCube = zeros(mp.dm1.Nact, mp.dm1.Nact, mp.dm1.NbasisModes);
        for ibm = 1:mp.dm1.NbasisModes/2
            [mp.dm1.basisCube(:, :, ibm), mp.dm1.basisCube(:, :, ibm+mp.dm1.NbasisModes/2)] = ...
                falco_make_fourier_mode(mp.dm1.Nact, mp.dm1.NactBeam, ...
                mp.dm1.fourier_basis_xis(ibm), mp.dm1.fourier_basis_etas(ibm));            
        end
        
        if length(mp.dm2.fourier_basis_xis) ~= length(mp.dm2.fourier_basis_etas)
           error('fourier_basis_xis and fourier_basis_etas must have the same length')
        end
        mp.dm2.NbasisModes = 2 * length(mp.dm2.fourier_basis_xis); % 2 for sin and cos at each freq
        mp.dm2.basisCube = zeros(mp.dm2.Nact, mp.dm2.Nact, mp.dm2.NbasisModes);
        for ibm = 1:mp.dm2.NbasisModes/2
            [mp.dm2.basisCube(:, :, ibm), mp.dm2.basisCube(:, :, ibm+mp.dm2.NbasisModes/2)] = ...
                falco_make_fourier_mode(mp.dm2.Nact, mp.dm2.NactBeam, ...
                mp.dm2.fourier_basis_xis(ibm), mp.dm2.fourier_basis_etas(ibm));            
        end
        
    otherwise
        error('Value of the DM basis type not allowed.')
end

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


mp.dm2.compact.dummy = 1;
mdc = mp.dm2.compact;
mp.dm2.compact = mp.dm2; % BE CAREFUL ABOUT OVERWRITING VARIABLES
for fn = fieldnames(mdc)'
    mp.dm2.compact.(fn{1}) = mdc.(fn{1});
end

mp.dm2.dx = mp.P2.full.dx;
mp.dm2.compact.dx = mp.P2.compact.dx;
mp.dm2 = falco_gen_dm_poke_cube(mp.dm2, mp, mp.P2.full.dx, 'NOCUBE');
if( any(mp.dm_ind==2) ) 
    mp.dm2.compact = falco_gen_dm_poke_cube(mp.dm2.compact, mp, mp.P2.compact.dx);
else
    mp.dm2.compact = falco_gen_dm_poke_cube(mp.dm2.compact, mp, mp.P2.compact.dx,'NOCUBE');
end

if(isfield(mp.dm2,'V')==false); mp.dm2.V = zeros(mp.dm2.Nact,mp.dm2.Nact); end %--Initial DM voltages

%%
%--Re-include all actuators in the basis set. Need act_ele to be a column vector.
if any(mp.dm_ind == 1); mp.dm1.act_ele = (1:mp.dm1.NactTotal).'; end
if any(mp.dm_ind == 2); mp.dm2.act_ele = (1:mp.dm2.NactTotal).'; end
if any(mp.dm_ind == 8); mp.dm8.act_ele = (1:mp.dm8.NactTotal).'; end
if any(mp.dm_ind == 9); mp.dm9.act_ele = (1:mp.dm9.NactTotal).'; end
%--Update the number of elements used per DM
if any(mp.dm_ind == 1); mp.dm1.Nele = mp.dm1.NbasisModes; else; mp.dm1.Nele = 0; end
if any(mp.dm_ind == 2); mp.dm2.Nele = mp.dm2.NbasisModes; else; mp.dm2.Nele = 0; end
if any(mp.dm_ind == 8); mp.dm8.Nele = length(mp.dm8.act_ele); else; mp.dm8.Nele = 0; end
if any(mp.dm_ind == 9); mp.dm9.Nele = length(mp.dm9.act_ele); else; mp.dm9.Nele = 0; end

end