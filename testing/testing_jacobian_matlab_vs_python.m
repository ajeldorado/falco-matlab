a = cat(3, ...
    [1 4 7 10; 1 4 7 10; 13 16 19 22; 13 16 19 22; 25 28 31 34; 25 28 31 34; 37 40 43 46; 37 40 43 46; 49 52 55 58; 49 52 55 58], ...
    [2 5 8 11; 2 5 8 11; 14 17 20 23; 14 17 20 23; 26 29 32 35; 26 29 32 35; 38 41 44 47; 38 41 44 47; 50 53 56 59; 50 53 56 59], ...
    [3 6 9 12; 3 6 9 12; 15 18 21 24; 15 18 21 24; 27 30 33 36; 27 30 33 36; 39 42 45 48; 39 42 45 48; 51 54 57 60; 51 54 57 60]);
%%
npixel = size(a,1)/2; %5
nact_tot = size(a,2);
nsbp = size(a,3);
%%
%G_reordered = reshape(a, [2 * npixel * nsbp, nact_tot]) %old way used in bb_jacobian that is possibly not doing the right thing

G_reordered = reshape(permute(a, [1 3 2]), 2*npixel*nsbp, nact_tot);

i_s = 1:2 * npixel : (2 * nsbp * npixel - 2 * npixel + 1); %start (at 1 because matlab) : step : end (excludes last value)
i_e = i_s + 1; %to consider the other Re or Im part of each pixel 
i_se = reshape([i_s; i_e], [], 1); % [i_s; i_e] is a 2 row matrix with i_s at 1st row and i_e at the 2nd ; reshape flattens it in a column vector stacking each previous columns  
                                       %reshape makes it 1 vector column
                                       %with i_s stacked over i_e
                                     

G_factor = []; %a matrix where the rows are of in the following order: all of pixel 1 for each wavelength, all of pixel 2 for each wavelength, etc
for factor = 0:2:(2 * npixel - 2) %loops through each pixel (step = 2 for Im and Re parts)
    G_factor = [G_factor; G_reordered(i_se + factor, :)]; %concatenate the rows vertically 
end
% disp('G_factor')
% disp(G_factor)

G = reshape(G_factor, [2*npixel*nsbp, nact_tot])

%%
nprobes = 8;
%elements = 1:(2*npixel * nsbp *nprobes );

% Reshape the vector into the desired 3D matrix
%m = reshape(elements, [2*npixel* nsbp, nprobes])
%size(m)
% 
% m_transpose = reshape(m.', nprobes, 2*nsbp, npixel)
% size(m_transpose)
% m2 = permute(m_transpose, [2, 1, 3])
% size(m2)

% Create a vector from 1 to 32
array = 1:(nact_tot*nprobes);

% Reshape the vector into an 8x4 matrix, then transpose it to get 8 rows and 4 columns filled row-wise
m = reshape(array, [nact_tot, nprobes]); % = probes in python
size(G)
size(m)
Gdotm1 = 4*G*m; 
Gdotm2 = reshape(Gdotm1.', nprobes, 2*nsbp, npixel) %same as gdus in python
Gdotm3 = permute(Gdotm2, [2, 1, 3]) %2*Nsbp_bb x Npairs x Npix %same as gdus.transpose((0,2,1))
size(Gdotm3)   
%%

m_pinv = zeros(size(Gdotm3));
for i = 1:size(Gdotm3, 3) %for each page (corresponds to a single pixel)
    m_pinv(:, :, i) = (pinv(Gdotm3(:, :, i))).';
end
disp('m_pinv')
size(m_pinv)
disp(m_pinv)

%%
% zAll_rand = rand(nprobes, npixel)
% zAll4_t_rand = reshape(4*zAll_rand.', [nprobes, 1,npixel])
% size(zAll4_t_rand)


array_zall = 100:(100 + nprobes*npixel - 1);

zall1 = reshape(array_zall, [nprobes, npixel]);

zall4 = reshape(zall1, [nprobes, 1, npixel])

epix = []
for n = 1:2:2*nsbp 
    m_pinv_per_sbp = m_pinv(n:n+1,:, :);
    epix_per_sbp = pagemtimes(m_pinv_per_sbp, zall4);
    epix = [epix, epix_per_sbp]
end
Eest1 = epix(1,:,:) + 1i * epix(2,:,:)
Eest = permute(Eest1, [3,2,1])

%%
Eest_notComplex_rand = pagemtimes(m_pinv, zall4) %E_hat %shape is 2*Nsbp x 1 x Npix 
%%
%reshape to have Npix x Nsbp complex matrix
Eest_complex1_rand = reshape(Eest_notComplex_rand, [2, nsbp, npixel]) %separates Nwpsbp in Re and Im
Eest_complex2_rand = permute(Eest_complex1_rand, [3, 2, 1]) % Npix x Nsbp x 2
Eest = Eest_complex2_rand(:, :, 1) + 1i * Eest_complex2_rand(:, :, 2) % Npix x Nsbp_bb complex matrix


