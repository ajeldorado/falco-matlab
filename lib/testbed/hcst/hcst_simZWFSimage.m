
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

end %--END OF FUNCTION

