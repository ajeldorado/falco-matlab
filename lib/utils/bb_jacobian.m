

function G_bb = bb_jacobian(mp,jacStruct,dm_ind)

    %Computes the Jacobian matrix for broadband Pair-Wise Probing.

    G_tot = rearrange_jacobians(mp,jacStruct,dm_ind);


    Nele_tot = size(G_tot, 2); 

    
    %G_permuted = permute(G_tot, [3 2 1]); %each page is a pixel (one page is Nsbp x Nele_tot)
    %G = reshape(permute(G_permuted, [2 1 3]), size(G_permuted, 2), [])'; %stack each pages on top of the other of G_permuted
    G_permuted = zeros(2 * mp.Nsbp, Nele_tot, mp.Fend.corr.Npix); %each page is a pixel (one page is 2*Nsbp x Nele_tot)
    for sbp = 1:mp.Nsbp_bb 
        real_parts = G_tot(1:2:end, :, sbp);  %odd rows
        imag_parts = G_tot(2:2:end, :, sbp);  %even rows
    
        % Empiler les parties réelle et imaginaire
        G_permuted((sbp-1)*2 + (1:2), :, :) = permute(cat(3, real_parts, imag_parts), [3, 2, 1]);
    end

    G_bb = reshape(permute(G_permuted, [2 1 3]), size(G_permuted, 2), [])'; %stack each pages on top of the other of G_permuted

% 
%     G_reordered = reshape(G_tot, [2 * mp.Fend.corr.Npix * mp.Nsbp_bb, Nele_tot]); %collapses x (pixels) and y (actuators) into a single dimension
% 
%     % Reshape in group of pixel instead of wavelength 
%     i_s = 1:2 * mp.Fend.corr.Npix : (2 * mp.Nsbp_bb * mp.Fend.corr.Npix - 2 * mp.Fend.corr.Npix + 1); %start (at 1 because matlab) : step : end (excludes last value)
%     i_e = i_s + 1; %to consider the other Re or Im part of each pixel 
%     i_se = reshape([i_s; i_e], [], 1); % [i_s; i_e] is a 2 row matrix with i_s at 1st row and i_e at the 2nd ; reshape flattens it in a column vector stacking each previous columns  
%                                        %reshape makes it 1 vector column
%                                        %with i_s stacked over i_e
%                                      
% 
%     G_factor = []; %a matrix where the rows are of in the following order: all of pixel 1 for each wavelength, all of pixel 2 for each wavelength, etc
%     for factor = 0:2:(2 * mp.Fend.corr.Npix - 2) %loops through each pixel (step = 2 for Im and Re parts)
%         G_factor = [G_factor; G_reordered(i_se + factor, :)]; %concatenate the rows vertically 
%     end
% 
%     G = reshape(G_factor, [2 * mp.Fend.corr.Npix * mp.Nsbp_bb, Nele_tot]); %re-collapses x (pixels) and y (actuators) into a single dimension

end

