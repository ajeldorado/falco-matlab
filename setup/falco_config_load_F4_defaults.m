% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_load_Fend.defaults(mp)
%
% Function to store mask configuration data for apodized coronagraphs.
% -Add more apodizers as desired.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-06-07 by A.J. Riggs.




function mp = falco_config_load_Fend.defaults(mp)

%%--Final Focal Plane (Fend. Properties

%--Specs for Correction (Corr) region and the Scoring (Score) region.
% if(isfield(mp.Fend.corr,'Rin')==false); mp.Fend.corr.Rin  = mp.F3.Rin; end  %--lambda0/Dcircumscribed, inner radius of correction region
if(isfield(mp.Fend.corr,'Rin')==false)
    if(isfield(mp.F3,'RinA')) 
        mp.Fend.corr.Rin  = mp.F3.RinA;
    else
        mp.Fend.corr.Rin  = mp.F3.Rin;
    end
end  %--lambda0/D, inner radius of correction region
if(isfield(mp.Fend.score,'Rin')==false); mp.Fend.score.Rin = mp.Fend.corr.Rin; end  %--Needs to be >= that of Correction mask
if(isfield(mp.Fend.corr,'Rout')==false); mp.Fend.corr.Rout  = min( [floor(mp.dm1.Nact/2*(1-mp.fracBW/2)), mp.F3.Rout ] ); end %--lambda0/Dcircumscribed, outer radius of correction region
if(isfield(mp.Fend.score,'Rout')==false); mp.Fend.score.Rout = mp.Fend.corr.Rout; end  %--Needs to be <= that of Correction mask
if(isfield(mp.Fend.corr,'ang')==false); mp.Fend.corr.ang  = 180; end  %--degrees per side
if(isfield(mp.Fend.score,'ang')==false); mp.Fend.score.ang = 180; end  %--degrees per side

if(isfield(mp.Fend.'sides')==false); mp.Fend.sides = 'both'; end  %--options: 'left', 'right','top','bottom'; any other values produce an annular region 


%%--Final Focal Plane (Fend. Properties
if(isfield(mp.Fend.compact,'res')==false); mp.Fend.res = 3; end  %--Pixels per lambda_c/D
if(isfield(mp.Fend.full,'res')==false); mp.Fend.full.res = 6; end  %--Pixels per lambda_c/D
if(isfield(mp.Fend.'FOV')==false); mp.Fend.FOV = 1 + mp.Fend.corr.Rout; end  % minimum desired field of view (along both axes) in lambda0/D


end %--END OF FUNCTION



