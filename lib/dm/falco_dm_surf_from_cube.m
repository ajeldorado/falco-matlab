% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to produce a DM surface (in meters). Uses linear superposition
% of a datacube of influence functions.

function DMsurf = falco_dm_surf_from_cube(dm,dmFullOrCompact)

    DMsurf = zeros(dmFullOrCompact.NdmPad); % Initialize the empty array
    for iact=1:dm.NactTotal 
        if( any(any(dmFullOrCompact.inf_datacube(:,:,iact))) && any(dm.VtoH(iact)) )
            y_box_ind = dmFullOrCompact.xy_box_lowerLeft(1,iact):dmFullOrCompact.xy_box_lowerLeft(1,iact)+dmFullOrCompact.Nbox-1; % x-indices in pupil arrays for the box
            x_box_ind = dmFullOrCompact.xy_box_lowerLeft(2,iact):dmFullOrCompact.xy_box_lowerLeft(2,iact)+dmFullOrCompact.Nbox-1; % y-indices in pupil arrays for the box
            DMsurf(x_box_ind,y_box_ind) = DMsurf(x_box_ind,y_box_ind) + dm.V(iact)*dm.VtoH(iact)*dmFullOrCompact.inf_datacube(:,:,iact);
        end
    end
    
    %--Adjust the orientation if specified
    if(isfield(dm,'fliplr'))
        if(dm.fliplr)
            DMsurf = fliplr(DMsurf);
        end
    end
    
    if(isfield(dm,'flipud'))
        if(dm.flipud)
            DMsurf = flipud(DMsurf);
        end
    end

end