% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function dm = falco_gen_dm_poke_cube_PROPER(dm,mp)
% 
%--Function to compute the influence functions for a DM and save them as a
%  datacube. Uses PROPER to generate the DM surface.
%
%--Method:
%--STEPS FOR GENERATING THE DM INFLUENCE FUNCTION DATACUBE WITH PROPER
% 1) Determine and save the final, even-valued size of the square DM surface array. 
%    Should fully contain all the actuators (when not rotated). 
% 2) Determine the (even-valued) size, Nbox, of the square sub-array to crop out that will contain the poke.
% 3) Compute the data cube of influence functions
%  a) Use PROPER to generate each normalized poke one at a time.
%  b) Find the 2-D location of the absolute-value peak of the poke.
%  c) Compute the (x,y) indices of the lower-left corner pixel of the box containing the poke.
%  d) Save out the sub-array of the poked actuator.
%
%--INPUTS
% dm: structure of DM parameters
% mp: structure of optical model parameters. Needed to compute the correct
%     padding for angular spectrum propagation.
%
%--OUTPUT
% dm: structure of DM parameters
% 
%--Important, new fields of output dm structure:
% dm.inf_datacube = datacube of influence functions
% dm.Ndm: number of points across the full DM surface with no extra padding
% dm.NdmPad = number of points across the full DM surface. Padded to enable
% sub-indexing of influence functions padded for A.S. propagation.
% dm.xy_box_lowerLeft = lower-left pixel positions of the influence functions in the padded
%   DM surface array. This is essential to adding the influence functions
%   together later via superposition. Units of pixels assuming dm.Ndmpad
%   points across.
%
%--VERSION HISTORY
% -Created falco_gen_dm_poke_cube_PROPER.m on 2017-11-17 by A.J. Riggs. 
%

function dm = falco_gen_dm_poke_cube_PROPER(dm,mp,varargin)

% % -------------------------------------------------------------------- %
%--DEBUGGING: HARD-CODED VALUES
% clear;
% addpath ~/Repos/FALCO/utils
% addpath ~/Repos/FALCO/dm
% addpath ~/Repos/FALCO/proper_v3.0.1_matlab_22aug17/
% dm.inf0 = fitsread('influence_dm5v2.fits');    %  -1*fitsread('inf64.3.fits');                              
% dm.Nact = 48;
% dm.dx = 1/7*1e-3;
% dm.dx_inf0 = 1e-4;
% dm.dm_spacing = 1e-3;
% dm.xc = 23.5;
% dm.yc = 23.5;
% dm.xtilt = 0;
% dm.ytilt = 0;
% dm.zrot = 0;
% dm.centering = 'pixel';
% dm.flagZYX = false;
% mp.d_dm1_dm2 = 3;
% mp.sbp_centers = 550e-9*[0.9 1 1.1];
% % -------------------------------------------------------------------- %

%--Enable the ability to turn off the data cube calculation with a keyword.
flagGenCube = true;
icav = 0;             % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'nocube'}
        flagGenCube  = false; %--Don't compute the data cube of influence functions. (For the full model.)
      otherwise
        error('falco_gen_dm_poke_cube_PROPER: Unknown keyword: %s\n', varargin{icav});
    end
end

%--Set the order of operations
orderOfOps = 'XYZ';
if(isfield(dm,'flagZYX')); 
    if(dm.flagZYX)
        orderOfOps = 'ZYX'; 
    end
end


% 1) Determine and save the final, even-valued size of the square DM surface array. 
%    Fully contains all the actuators (when not rotated). 
Ninf0 = length(dm.inf0); %--Number of points across the influence function at its native resolution
Ndm0 = ceil_even( Ninf0 + (dm.Nact - 1)*(dm.dm_spacing/dm.dx_inf0) ); %--Number of points across the DM surface at native influence function resolution
dm.NdmMin = ceil_even( Ndm0*(dm.dx_inf0/dm.dx))+2; %--Number of points across the (un-rotated) DM surface at new, desired resolution.
dm.Ndm = ceil_even( max(abs([sqrt(2)*cosd(45-dm.zrot),sqrt(2)*sind(45-dm.zrot)]))*Ndm0*(dm.dx_inf0/dm.dx))+2; %--Number of points across the array to fully contain the DM surface at new, desired resolution and z-rotation angle.
% dm.Ndm = ceil_even( Ndm0*(dm.dx_inf0/dm.dx))+2; %--Number of points across the (un-rotated) DM surface at new, desired resolution.


% 2) Determine the (even-valued) sizes of the square sub-arrays containing the poked actuator. 
dm.Ninf = ceil_even( Ninf0*(dm.dx_inf0/dm.dx) ); %--Number of points across the nominal influence function at the new, desired resolution
dm.Nbox = dm.Ninf;
%--Also compute the padded size for the angular spectrum (AS) propagation between DMs 1 and 2. Set as a power of 2 for FFTs.
Nmin = min(mp.sbp_centers)*max(abs([mp.d_P2_dm1, mp.d_dm1_dm2, (mp.d_P2_dm1 + mp.d_dm1_dm2) ]))/dm.dx^2; %--Minimum number of points across for accurate angular spectrum propagation
dm.NboxAS = 2^(nextpow2(max([dm.Nbox,Nmin])));  %--Zero-pad for FFTs in A.S. propagation. Use a larger array if the max sampling criterion for angular spectrum propagation is violated


% 3) Compute the number of points across the DM to account for the padded influence function.
dm.NdmPad = ceil_even( dm.Ndm + (dm.NboxAS-dm.Nbox) ); %--Number of points across the DM surface (with padding for angular spectrum propagation) at new, desired resolution.
%--Compute coordinates (in meters) of the full DM array
if(strcmpi(dm.centering,'pixel')  ) 
    dm.x_pupPad = (-dm.NdmPad/2:(dm.NdmPad/2 - 1))*dm.dx; % meters, coords for the full DM arrays. Origin is centered on a pixel
else
    dm.x_pupPad = (-(dm.NdmPad-1)/2:(dm.NdmPad-1)/2)*dm.dx; % meters, coords for the full DM arrays. Origin is centered between pixels for an even-sized array
end
dm.y_pupPad = dm.x_pupPad;



% 4) Compute the data cube of influence functions
pupil_ratio = 1;
dm.NactTotal = dm.Nact^2;
wl_dummy = 1e-6; %--dummy value needed to initialize wavelength in PROPER (meters)

%--Adjust the centering of the output DM surface. The shift needs to be in
%units of actuators, not meters, for prop_dm.m.
Darray = dm.NdmPad*dm.dx;
Narray = dm.NdmPad;
switch dm.centering % 0 shift for pixel-centered pupil, or -Darray/2/Narray shift for inter-pixel centering
    case {'interpixel'}
        cshift = -Darray/2/Narray/dm.dm_spacing; 
    case {'pixel'}
        cshift = 0;
end

if(flagGenCube)
    %--Initialize vectors and matrices
    dm.xyc_box = zeros(2,dm.NactTotal);
    dm.xy_box_lowerLeft = zeros(2,dm.NactTotal);
    dm.inf_datacube = zeros(dm.Nbox,dm.Nbox,dm.NactTotal);

    fprintf('Influence function padded to 2^nextpow2(%d) = %d for A.S. propagation. \n',2*ceil(1/2*max([dm.Nbox,Nmin])),dm.NboxAS);
    fprintf('Computing datacube of DM influence functions... '); tic;
    
    dm.act_ele = []; %--Indices of actuators to use. 
    
    %--Loop over all actuators
    for iact = 1:dm.NactTotal

        %--Re-initialize DM commands and poke one actuator 
        dU = zeros(dm.Nact);
        dU(iact) = 1;

        % 4a) Use PROPER to generate each normalized poke one at a time.
        bm  = prop_begin(Darray, wl_dummy, Narray, pupil_ratio);
        [~,DMsurf] = prop_dm(bm, dU, dm.xc-cshift, dm.yc-cshift, dm.dm_spacing,'XTILT',dm.xtilt,'YTILT',dm.ytilt,'ZTILT',dm.zrot,orderOfOps);
%         figure(1); imagesc(DMsurf); axis xy equal tight; colorbar; 
%         figure(2); imagesc(log10(abs(DMsurf))); axis xy equal tight; colorbar; 

        % 4b) Find the 2-D location of the absolute-value peak of the poke.
        [~,inds] = max(abs(DMsurf(:)));
        [ymaxind,xmaxind] = ind2sub(size(DMsurf),inds);
        dm.xyc_box(:,iact) = [xmaxind,ymaxind];

        % 4c) Compute the (x,y) indices of the lower-left corner pixel of the box containing the poke.
        dm.xy_box_lowerLeft(:,iact)  = dm.xyc_box(:,iact) - dm.Nbox/2 + 1; % indices of pixel in lower left of the postage stamp within the whole padded DM surface array

        % 4d) Save out the sub-array of the poked actuator.
        y_inds_box = dm.xy_box_lowerLeft(1,iact):dm.xy_box_lowerLeft(1,iact)+dm.Nbox-1; % x-indices in pupil arrays for the box
        x_inds_box = dm.xy_box_lowerLeft(2,iact):dm.xy_box_lowerLeft(2,iact)+dm.Nbox-1; % y-indices in pupil arrays for the box

        dm.act_ele = [dm.act_ele; iact]; %--Keep this actuator
        dm.inf_datacube(:,:,iact) = DMsurf(x_inds_box,y_inds_box);
%         figure(3); imagesc(log10(abs(dm.inf_datacube(:,:,iact)))); axis xy equal tight; colorbar; pause(1/20);

%         %--Don't use the actuator if it goes out of bounds at all.
%         if( (any(x_inds_box<=0)==false) && (any(y_inds_box<=0)==false) )
%             dm.act_ele = [dm.act_ele; iact]; %--Keep this actuator
%             dm.inf_datacube(:,:,iact) = DMsurf(x_inds_box,y_inds_box);
%             % figure(3); imagesc(log10(abs(dm.inf_datacube(:,:,iact)))); axis xy equal tight; colorbar; pause(1/20);
%         end
        
    end
    fprintf('done.\tTime = %.2fs\n',toc);

else
    dm.act_ele = (1:dm.NactTotal).';
end


end %--END OF FUNCTION



