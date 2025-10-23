
function inff = inf_modify(inf, width_percent, pixel_shift)
    N       = length(inf);
    met.norm = norm(inf); 

    x_1do       = 1:N;
    x_1d        = linspace(1+width_percent/100,N-width_percent/100,N) + pixel_shift(1);
    y_1d        = linspace(1+width_percent/100,N-width_percent/100,N) + pixel_shift(2);
    inff        = interp2(x_1do, x_1do', inf, x_1d,y_1d','spline'); 
    aa          = find (isnan(inff));
    inff(aa)    = 0;
    inff        = inff * met.norm/norm(inff);
return
