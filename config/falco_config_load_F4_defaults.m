% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_load_F4_defaults(mp)
%
% Function to store mask configuration data for apodized coronagraphs.
% -Add more apodizers as desired.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-06-07 by A.J. Riggs.




function mp = falco_config_load_F4_defaults(mp)

%%--Final Focal Plane (F4) Properties

%--Specs for Correction (Corr) region and the Scoring (Score) region.
% if(isfield(mp.F4.corr,'Rin')==false); mp.F4.corr.Rin  = mp.F3.Rin; end  %--lambda0/Dcircumscribed, inner radius of correction region
if(isfield(mp.F4.corr,'Rin')==false)
    if(isfield(mp.F3,'RinA')) 
        mp.F4.corr.Rin  = mp.F3.RinA;
    else
        mp.F4.corr.Rin  = mp.F3.Rin;
    end
end  %--lambda0/D, inner radius of correction region
if(isfield(mp.F4.score,'Rin')==false); mp.F4.score.Rin = mp.F4.corr.Rin; end  %--Needs to be >= that of Correction mask
if(isfield(mp.F4.corr,'Rout')==false); mp.F4.corr.Rout  = min( [floor(mp.dm1.Nact/2*(1-mp.fracBW/2)), mp.F3.Rout ] ); end %--lambda0/Dcircumscribed, outer radius of correction region
if(isfield(mp.F4.score,'Rout')==false); mp.F4.score.Rout = mp.F4.corr.Rout; end  %--Needs to be <= that of Correction mask
if(isfield(mp.F4.corr,'ang')==false); mp.F4.corr.ang  = 180; end  %--degrees per side
if(isfield(mp.F4.score,'ang')==false); mp.F4.score.ang = 180; end  %--degrees per side

if(isfield(mp.F4,'sides')==false); mp.F4.sides = 'both'; end  %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%%--Final Focal Plane (F4) Properties
if(isfield(mp.F4.compact,'res')==false); mp.F4.res = 3; end  %--Pixels per lambda_c/D
if(isfield(mp.F4.full,'res')==false); mp.F4.full.res = 6; end  %--Pixels per lambda_c/D
if(isfield(mp.F4,'FOV')==false); mp.F4.FOV = 1 + mp.F4.corr.Rout; end  % minimum desired field of view (along both axes) in lambda0/D


end %--END OF FUNCTION



