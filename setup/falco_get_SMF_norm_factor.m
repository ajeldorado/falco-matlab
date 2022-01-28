% Function to compute the normalization value for the compact and full
% models. 
function mp = falco_get_SMF_norm_factor(mp)
I00Fiber = zeros(mp.Fend.Nfiber,1);
for II = mp.Fend.Nfiber
    I00FiberII = 0;
    for wi = 1:mp.Nsbp        
        modvar.sbpIndex = wi;
        modvar.ttIndex = 1;
        modvar.wpsbpIndex = 1;
        modvar.x_offset = mp.Fend.x_fiber(II);
        modvar.y_offset = mp.Fend.y_fiber(II);
        modvar.whichSource = 'offaxis';
        [~, Efiber] = model_full(mp, modvar);
        Iplanetmax = abs(Efiber).^2;
        I00FiberII = I00FiberII + Iplanetmax;
    end
    I00Fiber(II) = I00FiberII;
end
mp.Fend.full.I00Fiber = I00Fiber;
end