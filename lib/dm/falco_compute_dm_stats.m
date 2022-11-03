% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Compute, store, and report statistics on the DM commands and surfaces.

function out = falco_compute_dm_stats(mp, out, Itr)

%--ID and OD of pupil used when computing DM actuation stats [units of pupil diameters]
OD_pup = 1.0;

%--Compute the DM surfaces and their change from the last iteration
if any(mp.dm_ind == 1)
    Vnow = mp.dm1.V;
    DM1surf = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);
    Vprev = mp.dm1.V - mp.dm1.dV;
    mp.dm1.V = Vprev;
    DM1surfPrev = falco_gen_dm_surf(mp.dm1, mp.dm1.compact.dx, mp.dm1.compact.Ndm);
    DM1surfDiff = DM1surf - DM1surfPrev;
    mp.dm1.V = Vnow; % reset        
end

if any(mp.dm_ind == 2)
    Vnow = mp.dm2.V;
    DM2surf = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);
    Vprev = mp.dm2.V - mp.dm2.dV;
    mp.dm2.V = Vprev;
    DM2surfPrev = falco_gen_dm_surf(mp.dm2, mp.dm2.compact.dx, mp.dm2.compact.Ndm);
    DM2surfDiff = DM2surf - DM2surfPrev;
    mp.dm2.V = Vnow; % reset
end

%--Compute coordinates
if any(mp.dm_ind == 1)
    Ndm = length(DM1surf);
elseif any(mp.dm_ind == 2)
    Ndm = length(DM2surf);
end
dx_dm = mp.P2.compact.dx / mp.P2.D; %--Normalized dx. Units of pupil diameters.
switch mp.centering
    case 'interpixel'
        xs = ( -(Ndm-1)/2:(Ndm-1)/2 )*dx_dm;
    otherwise
        xs = ( -(Ndm/2):(Ndm/2-1) )*dx_dm;
end
[XS, YS] = meshgrid(xs);
RS = sqrt(XS.^2 + YS.^2);
rmsSurfInd = find(RS>=mp.P1.IDnorm/2 & RS<=OD_pup/2);

%--Calculate and report updated P-V DM voltages.
if any(mp.dm_ind == 1)
    out.dm1.Vpv(Itr) = max(mp.dm1.V(:)) - min(mp.dm1.V(:));
    fprintf(' DM1 P-V in volts: %.3f \n', out.dm1.Vpv(Itr)); 
end
if any(mp.dm_ind == 2)
    out.dm2.Vpv(Itr) = max(mp.dm2.V(:)) - min(mp.dm2.V(:));
    fprintf(' DM2 P-V in volts: %.3f \n', out.dm2.Vpv(Itr)); 
end
if any(mp.dm_ind == 8)
    out.dm8.Vpv(Itr) = (max(max(mp.dm8.V))-min(min(mp.dm8.V)));
    Nrail8 = length(find( (mp.dm8.V <= mp.dm8.Vmin) | (mp.dm8.V >= mp.dm8.Vmax) ));
    fprintf(' DM8 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm8.Vpv(Itr), Nrail8,mp.dm8.NactTotal,100*Nrail8/mp.dm8.NactTotal); 
end
if any(mp.dm_ind == 9)
    out.dm9.Vpv(Itr) = (max(max(mp.dm9.V))-min(min(mp.dm9.V)));
    Nrail9 = length(find( (mp.dm9.V <= mp.dm9.Vmin) | (mp.dm9.V >= mp.dm9.Vmax) ));
    fprintf(' DM9 P-V in volts: %.3f\t\t%d/%d (%.2f%%) railed actuators \n', out.dm9.Vpv(Itr), Nrail9,mp.dm9.NactTotal,100*Nrail9/mp.dm9.NactTotal); 
end

%--Calculate and report updated P-V and RMS DM surfaces, as well as the
% change (Delta) in those for the current iteration.
if any(mp.dm_ind == 1)
    out.dm1.Spv(Itr) = max(DM1surf(:)) - min(DM1surf(:));
    out.dm1.Srms(Itr) = falco_rms(DM1surf(rmsSurfInd));
    fprintf(' RMS DM1 surface = %.1f nm\n', 1e9*out.dm1.Srms(Itr))
    
    out.dm1.DeltaSpv(Itr) = max(DM1surfDiff(:)) - min(DM1surfDiff(:));
    out.dm1.DeltaSrms(Itr) = falco_rms(DM1surfDiff(rmsSurfInd));
end
if any(mp.dm_ind == 2)
    out.dm2.Spv(Itr) = max(DM2surf(:)) - min(DM2surf(:));
    out.dm2.Srms(Itr) = falco_rms(DM2surf(rmsSurfInd));
    fprintf(' RMS DM2 surface = %.1f nm\n', 1e9*out.dm2.Srms(Itr))
    
    out.dm2.DeltaSpv(Itr) = max(DM2surfDiff(:)) - min(DM2surfDiff(:));
    out.dm2.DeltaSrms(Itr) = falco_rms(DM2surfDiff(rmsSurfInd));
end

if any(mp.dm_ind == 1)
    fprintf(' Delta RMS DM1 surface = %.1f nm this iteration\n', 1e9*out.dm1.DeltaSrms(Itr))
end
if any(mp.dm_ind == 2)
    fprintf(' Delta RMS DM2 surface = %.1f nm this iteration\n', 1e9*out.dm2.DeltaSrms(Itr))
end

%--Report pinned and comoving actuators
if any(mp.dm_ind == 1)
    
    Npinned = length(mp.dm1.pinned);
    if Npinned > 0
        fprintf(' DM1 has %d pinned actuators.\n', Npinned)
    end
    
    Ngroups = length(mp.dm1.comovingGroups);
    Ncomoving = 0;
    for iGroup = 1:Ngroups
        Ncomoving = Ncomoving + length(mp.dm1.comovingGroups{iGroup});
    end
    if Ncomoving > 0
        fprintf(' DM1 has %d co-moving actuators in %d groups.\n', Ncomoving, Ngroups)
    end
    
    out.dm1.pinned{Itr} = mp.dm1.pinned;
    out.dm1.Vpinned{Itr} = mp.dm1.Vpinned;
    out.dm1.comovingGroups{Itr} = mp.dm1.comovingGroups;
    out.dm1.Npinned(Itr) = Npinned;
    out.dm1.Ncomoving(Itr) = Ncomoving;
    
end

if any(mp.dm_ind == 2)
    
    Npinned = length(mp.dm2.pinned);
    if Npinned > 0
        fprintf(' DM2 has %d pinned actuators.\n', Npinned)
    end
    
    Ngroups = length(mp.dm2.comovingGroups);
    Ncomoving = 0;
    for iGroup = 1:Ngroups
        Ncomoving = Ncomoving + length(mp.dm2.comovingGroups{iGroup});
    end
    if Ncomoving > 0
        fprintf(' DM2 has %d co-moving actuators in %d groups.\n', Ncomoving, Ngroups)
    end
    
    out.dm2.pinned{Itr} = mp.dm2.pinned;
    out.dm2.Vpinned{Itr} = mp.dm2.Vpinned;
    out.dm2.comovingGroups{Itr} = mp.dm2.comovingGroups;
    out.dm2.Npinned(Itr) = Npinned;
    out.dm2.Ncomoving(Itr) = Ncomoving;

end

end
