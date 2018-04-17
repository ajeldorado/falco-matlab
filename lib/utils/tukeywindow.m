function window = tukeywindow(n,alpha)
%   The function tukeywindow(n,alpha) returns an n-point, 1-D Tukey window.
% 
%   The alpha value, 0 < alpha < 1, specifies the width of the cosine taper.
%  
%  INPUTS:
%  n: number of points in column vector (in and out)
%  alpha: scalar value, 0 < alpha < 1, that defines the tapered part of the
%  window.
%
%  OUTPUT:
%  window: the 1-D window as an n-point column vector.
%    -Note that the Tukey window is centered on the array.
%--Initialize
window = ones(n,1);


if alpha >= 1
    ts = (0:1:(n-1)).';
    window = 1/2*( 1 - cos(2*pi*ts/(n-1)) );
else
    ts = linspace(0,1,n)';
    % Defines period of the taper as 1/2 period of a sine wave.
    period = alpha/2; 
    tl = floor(period*(n-1))+1;
    th = n-tl+1;
    % Ramp up
    window(1:tl) = ((1+cos(pi/period*(ts(1:tl) - period)))/2);
    % Ramp down
    window(th:end) = (1+cos(pi/period*(ts(th:end) - 1 + period)))/2;
    
end