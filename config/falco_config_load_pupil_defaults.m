% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_load_pupil_defaults(mp,DM)
%
% Function to store pupil mask configuration data for coronagraphs.
% -Add more apodizers as desired.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-05-29 by A.J. Riggs.




function mp = falco_config_load_pupil_defaults(mp)


%%--Pupil Masks
switch mp.whichPupil
    case 'Simple' % Can be used to create circular and annular apertures with radial spiders 
        
        
        if(isfield(mp.P1,'D')==false); mp.P1.D = 4; end  %--meters, diameter of telescope (This is like HabEx A)
        if(isfield(mp.P1,'Dfac')==false); mp.P1.Dfac = 1; end  %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 250; end  
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 250; end 
        if(isfield(mp.P4.full,'Nbeam')==false); mp.P4.full.Nbeam = 250; end  
        if(isfield(mp.P4.compact,'Nbeam')==false); mp.P4.compact.Nbeam = 200; end 
        
        if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        if(isfield(mp.P1,'ODnorm')==false); mp.P1.ODnorm = 1; end % Outer diameter (fraction of Nbeam) 
        
        if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.3; end % Inner radius of Lyot stop
        if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.85; end % Outer radius of Lyot stop 
        
        %if(isfield(mp.P4,'IDnorm')==false); mp.P4.IDnorm = 0; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        %if(isfield(mp.P4,'ODnorm')==false); mp.P4.ODnorm = 0.95; end % Outer diameter (fraction of Nbeam) 
        
        if(isfield(mp.P1,'num_strut')==false); mp.P1.num_strut = 0; end % Number of struts 
        if(isfield(mp.P1,'strut_angs')==false); mp.P1.strut_angs = []; end %Array of angles of the radial struts (deg)
        if(isfield(mp.P1,'strut_width')==false); mp.P1.strut_width = []; end  % Width of the struts (fraction of pupil diam.)
        
        if(isfield(mp.P4,'num_strut')==false); mp.P4.num_strut = 0; end % Number of struts 
        if(isfield(mp.P4,'strut_angs')==false); mp.P4.strut_angs = []; end %Array of angles of the radial struts (deg)
        if(isfield(mp.P4,'strut_width')==false); mp.P4.strut_width = []; end  % Width of the struts (fraction of pupil diam.)
                
    case{'WFIRST180718'}
        if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0.303; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)

        if(isfield(mp.P1,'D')==false); mp.P1.D = 2.3631; end  %--meters, diameter of telescope (used only for mas to lambda/D conversion)
        if(isfield(mp.P1,'Dfac')==false); mp.P1.Dfac = 1; end  %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        %if(isfield(mp.P1,'wStrut')==false); mp.P1.wStrut = 76e-3/mp.P1.D; end  %--76.0mm is the new value as of 2018-07-18
        if(isfield(mp.P4,'wStrut')==false); mp.P4.wStrut = 3.9/100.; end  

        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 324; end  %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 324; end 
        if(isfield(mp.P4.full,'Nbeam')==false); mp.P4.full.Nbeam = 250; end  
        if(isfield(mp.P4.compact,'Nbeam')==false); mp.P4.compact.Nbeam = 200; end 
    
    case{'WFIRST20180103','WFIRST_onaxis'}
        if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0.31; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)

        
        if(isfield(mp.P1,'D')==false); mp.P1.D = 2.3631; end  %--meters, diameter of telescope (used only for mas to lambda/D conversion)
        if(isfield(mp.P1,'Dfac')==false); mp.P1.Dfac = 1; end  %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        if(isfield(mp,'pup_strut_width')==false); mp.pup_strut_width = 3.22/100.; end  %--3.22% is the new value as of 2018-01-03
        if(isfield(mp,'LS_strut_width')==false); mp.LS_strut_width = 3.8/100.; end  

        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 324; end  %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 324; end 
        if(isfield(mp.P4.full,'Nbeam')==false); mp.P4.full.Nbeam = 250; end  
        if(isfield(mp.P4.compact,'Nbeam')==false); mp.P4.compact.Nbeam = 200; end 
        
        %if(isfield(mp.P4,'IDnorm')==false); mp.P4.IDnorm = 0.50; end  % value for original pupil was 2.46/8 = 0.3075
        %if(isfield(mp.P4,'ODnorm')==false); mp.P4.ODnorm = 0.80; end 
        
    case{'LUVOIRA5'}
        if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0.10; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        
        if(isfield(mp.P1,'D')==false); mp.P1.D = 15.2;  end %--meters, circumscribing diameter of telescope (used only for mas-to-lambda/D conversion)
        if(isfield(mp.P1,'Dfac')==false); mp.P1.Dfac = 15.2/13.7; end  %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
        
        if(isfield(mp,'LS_strut_width')==false); mp.LS_strut_width = 1.4/100.; end
        if(isfield(mp.P4,'wStrut'));  mp.LS_strut_width = mp.P4.wStrut;  end

        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 1000; end %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 500; end 
        if(isfield(mp.P4.full,'Nbeam')==false); mp.P4.full.Nbeam = 300; end  
        if(isfield(mp.P4.compact,'Nbeam')==false); mp.P4.compact.Nbeam = 200; end 
        

    case 'LUVOIR_B_offaxis'         % Note:  Nbeam needs to be >~500 to properly resolve segment gaps 
        if(isfield(mp.P1,'IDnorm')==false); mp.P1.IDnorm = 0; end % Inner diameter (fraction of Nbeam; zero if you want an off-axis telescope)
        
        if(isfield(mp.P1,'D')==false); mp.P1.D = 7.989; end  %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm.
        
        if(isfield(mp.P1.full,'Nbeam')==false); mp.P1.full.Nbeam = 1000; end  %--Number of pixels across the actual diameter of the beam/aperture (independent of beam centering)
        if(isfield(mp.P1.compact,'Nbeam')==false); mp.P1.compact.Nbeam = 500; end 
        if(isfield(mp.P4.full,'Nbeam')==false); mp.P4.full.Nbeam = 300; end  
        if(isfield(mp.P4.compact,'Nbeam')==false); mp.P4.compact.Nbeam = 200; end 
        
        %if(isfield(mp.P4,'IDnorm')==false); mp.P4.IDnorm = 0; end 
        %if(isfield(mp.P4,'ODnorm')==false); mp.P4.ODnorm = 0.80; end 
        
        
end


%--Resolution must be the same at pupils P1 and P4 for the Vortex and Lyot coronagraphs.
switch lower(mp.coro)
    case{'lc','hlc','vortex','vc','avc'}
        mp.P4.full.Nbeam = mp.P1.full.Nbeam;  % P4 must be the same as P1 for Vortex and H/LC. 
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;  % P4 must be the same as P1 for Vortex and H/LC. 
end


end %--END OF FUNCTION



