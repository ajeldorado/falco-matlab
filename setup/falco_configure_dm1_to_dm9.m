% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

function mp = falco_configure_dm1_to_dm9(mp)

    mp.dm1.NactTotal=0; mp.dm2.NactTotal=0; mp.dm3.NactTotal=0; mp.dm4.NactTotal=0; mp.dm5.NactTotal=0; mp.dm6.NactTotal=0; mp.dm7.NactTotal=0; mp.dm8.NactTotal=0; mp.dm9.NactTotal=0; 
    mp.dm1.Nele=0; mp.dm2.Nele=0; mp.dm3.Nele=0; mp.dm4.Nele=0; mp.dm5.Nele=0; mp.dm6.Nele=0; mp.dm7.Nele=0; mp.dm8.Nele=0; mp.dm9.Nele=0; %--Initialize for Jacobian calculations later. 

    %--Intialize delta DM voltages. Needed for Kalman filters.
    if(any(mp.dm_ind==1));  mp.dm1.dV = 0;  end
    if(any(mp.dm_ind==2));  mp.dm2.dV = 0;  end
    if(any(mp.dm_ind==3));  mp.dm3.dV = 0;  end
    if(any(mp.dm_ind==4));  mp.dm4.dV = 0;  end
    if(any(mp.dm_ind==5));  mp.dm5.dV = 0;  end
    if(any(mp.dm_ind==6));  mp.dm6.dV = 0;  end
    if(any(mp.dm_ind==7));  mp.dm7.dV = 0;  end
    if(any(mp.dm_ind==8));  mp.dm8.dV = 0;  end
    if(any(mp.dm_ind==9));  mp.dm9.dV = 0;  end

    %%--Intialize tied actuator pairs if not already defined. 
    % Dimensions of the pair list is [Npairs x 2]
    if(any(mp.dm_ind==3)); if(isfield(mp.dm3,'tied')==false); mp.dm3.tied = []; end; end
    if(any(mp.dm_ind==4)); if(isfield(mp.dm4,'tied')==false); mp.dm4.tied = []; end; end
    if(any(mp.dm_ind==5)); if(isfield(mp.dm5,'tied')==false); mp.dm5.tied = []; end; end
    if(any(mp.dm_ind==6)); if(isfield(mp.dm6,'tied')==false); mp.dm6.tied = []; end; end
    if(any(mp.dm_ind==7)); if(isfield(mp.dm7,'tied')==false); mp.dm7.tied = []; end; end
    if(any(mp.dm_ind==8)); if(isfield(mp.dm8,'tied')==false); mp.dm8.tied = []; end; end
    if(any(mp.dm_ind==9)); if(isfield(mp.dm9,'tied')==false); mp.dm9.tied = []; end; end

end