%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

function radi = radpix( nx, ny )

  ix2  = [-floor(nx / 2 ) : floor((nx - 1) / 2)].^2;
  iy2  = [-floor(ny / 2 ) : floor((ny - 1) / 2)].^2;
  [ax2, ay2]   = meshgrid(ix2, iy2);
  radi = sqrt(ax2 + ay2);

end