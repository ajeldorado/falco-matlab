
function Isum = hcst_simZWFSimage(mp)


%--Compute the DM surfaces outside the full model to save lots of time
if(any(mp.dm_ind==1)); mp.dm1.surfM = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,mp.dm1.NdmPad); end

Isum = 0; % Initialize image

for si=1:mp.Nsbp
    modvar.sbpIndex = si;
    
    for wi=1:mp.Nwpsbp
        modvar.whichSource = 'star';
        modvar.wpsbpIndex = wi;
        Etemp = model_hcst_ZWFS(mp, modvar);
        
        Isum = Isum + (abs(Etemp).^2)*mp.sbp_weights(si)*mp.full.lambda_weights(wi);
    end 
    
    if(mp.planetFlag)
        modvar.whichSource = 'exoplanet';
        Eout = model_hcst_ZWFS(mp,modvar);
        ImPlanetC = abs(Eout).^2; % In contrast
        Isum = Isum + ImPlanetC*mp.sbp_weights(si);
    end

end



% %--Loop over wavelengths within each sub-bandpass
% for ilam = 1:mp.full.Nlam  % Add intensities from all incoherent sources separately
%     
%     
%     %--Starlight
%     modvar.ebpIndex = ilam;
%     modvar.lambda = mp.full.lambdas(ilam);
%     modvar.whichSource = 'star';
%     Eout = model_full(mp, modvar);
%     Isum = Isum + mp.full.all_weights(ilam)*(abs(Eout).^2);
% %     Isum = Isum + mp.entireBWimage.Iweights(ilam)/IweightsSum*(abs(Eout).^2);
%     %Isum = Isum + (abs(Eout).^2)*(mp.jac.weights(tsi)/mp.Wsum)/mp.Nwpsbp;
% 
% %         %--Exoplanet light
% %         if(mp.planetFlag)
% %             modvar.whichSource = 'exoplanet';
% %             Eout = model_full(mp,modvar);
% %             ImPlanetC = abs(Eout).^2; % In contrast
% %             %Isum = Isum + ImPlanetC*(mp.jac.weights(tsi)/mp.Wsum)/mp.Nwpsbp;
% %         end
% 
% end
    
end %--END OF FUNCTION

