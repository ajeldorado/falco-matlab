
function inff = inf_modify(inf, width_percent, pixel_shift)
    N       = length(inf);
    met.max = max(inf(:)); 

    x_1do       = linspace(-N/2,N/2,N);
    x_1d        = linspace(-N/2+N/2*width_percent/100,N/2-N/2*width_percent/100,N) + pixel_shift(1);
    y_1d        = linspace(-N/2+N/2*width_percent/100,N/2-N/2*width_percent/100,N) + pixel_shift(2);
    inff        = interp2(x_1do, x_1do', inf, x_1d,y_1d','spline',0); 
    aa          = find (isnan(inff));
    inff(aa)    = 0;
    inff        = inff * met.max/max(inff(:));
return
