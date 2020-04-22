function [mp] = falco_gen_FPM_PIAACMC(mp)
% Generate the focal plane mask for the PIAACMC coronagraph.  The PIAACMC
% mask consists of a pattern of sub-lambda/D hexes which are at different
% heights and mirror-coated, meaning it should affect phase only


        %--Make or read in focal plane mask (FPM) amplitude for the full model
        FPMgenInputs.pixresFPM = mp.F3.full.res; %--pixels per lambda_c/D
        FPMgenInputs.centering = mp.centering;
        
%         mp.F3.full.mask.sag = zeros(2056);
        mp.F3.full.mask.sag = rot90(fitsread(strcat(mp.path.falco, 'lib/masks/PIAACMC/cmcMask.fits')), 0);
        mp.F3.full.mask.sag = falco_bin_downsample(mp.F3.full.mask.sag, mp.F3.downSampFac);
        mp.F3.compact.mask.sag = mp.F3.full.mask.sag;
        
        mp.F3.full.Nxi = size(mp.F3.full.mask.sag,2);
        mp.F3.full.Neta= size(mp.F3.full.mask.sag,1);
        mp.F3.compact.Nxi = size(mp.F3.compact.mask.sag,2);
        mp.F3.compact.Neta= size(mp.F3.compact.mask.sag,1);
        
end