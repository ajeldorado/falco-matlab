function [ Z, n, m ] = propcustom_gen_zernike( noll_index, apRad, RHO, THETA  )
%propcustom_gen_zernike Generates a Zernike polynomial in 2D array with
%coordinates defined by RHO, THETA (meshgrid in polar coordinates) and the
%beam radius, apRad. 

    % Convert from noll index to Zernike indices
    [ n, m ] = propcustom_zernikes_noll2index( noll_index );

    % Get Zernike polynomial
    Z = zeros(size(RHO));
    Z0 = zernfun(n,m,RHO(RHO<=apRad)/apRad,THETA(RHO<=apRad),'norm');
    Z(RHO<=apRad) = sqrt(pi)*Z0(:,1);

end

