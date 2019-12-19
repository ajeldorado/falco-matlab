function [ n, m ] = propcustom_zernikes_noll2index( noll_index )
%propcustom_zernikes_noll2index Converts the Noll index to Zernike indices.
% i.e. Z_j -> Z_{n,m}
% 
% Helper function for using propcustom_zernikes.m 

    % Convert from noll index to Zernike indices
    n = 0;
    j1 = noll_index-1;
    while(j1 > n)
        n = n + 1;
        j1 = j1 - n;
    end

    m = (-1)^noll_index  * (mod(n,2) + 2 * floor((j1+mod(n+1,2)) / 2.0 ));

end

