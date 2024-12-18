



function Eest_bb = pairwise_bb_estimation(mp, jacStruct, DM1Vplus, DM2Vplus, zAll)
    
    DM1Vnom = mp.dm1.V;
    DM2Vnom = mp.dm2.V;
    whichDM = mp.est.probe.whichDM;
    disp('here in pairwise_bb_estimation')

    %Npix = mp.Fend.corr.Npix;
    %Npairs = mp.est.probe.Npairs;

    %mp.Nwpsbp  = mp.Nsbp; %erase this if run another EXAMPLE_try_running_falco with pm.Nwpsbp = 3


    %Where does the 96x1 comes from?? otherwise can't do G_bb*dV
    %mp.dm1.act_ele = (1:mp.dm1.NactTotal).';
    %mp.dm2.act_ele = (1:mp.dm2.NactTotal).';

    gdus = [];%Initialize delta E-fields for each probe image
    
    G_bb = bb_jacobian(mp, jacStruct, mp.dm_ind); % 2*Nsbp*Npix x Nele_tot

    dV1 = []; 
    dV2 = [];
    dV = [];
    zAll4_t = reshape(4*zAll, [mp.est.probe.Npairs, 1, mp.Fend.corr.Npix]);
         
    for iProbe = 1:mp.est.probe.Npairs

        %zAll4_t = reshape(4*zAll.', [mp.est.probe.Npairs, 1, mp.Fend.corr.Npix]);
         

        %dV1minus = [];
        %dV2minus = [];

        if whichDM == 1
            dV1 = DM1Vplus(:, :, iProbe) - DM1Vnom; %dVplus size is nb of pixels (Npix) x nb of pairs of probes (Npairs) (right?)
            dV1 = dV1(mp.dm1.act_ele); %take the active elements (actuators) of dV1 (changes shape to vector column)
            %dE1plus = (4*squeeze(jacStruct.G1(:, :, modeIndex))*dV1plus(mp.dm1.act_ele)).*eye(2*Npix); %is eye the right format for the multiplication?
        
           
    
            %dV1minus = DM1Vminus(:, :, iProbe) - DM1Vnom;
            %dV1minus = dV1minus(mp.dm1.act_ele);
            %dE1minus = 4*squeeze(jacStruct.G1(:, :, modeIndex))*dV1minus(mp.dm1.act_ele).*eye(2*Npix);
    
        elseif whichDM == 2
            dV2 = DM2Vplus(:, :, iProbe) - DM2Vnom; %dVplus size is nb of pixels (Npix) x nb of pairs of probes (Npairs) (right?)
            dV2 = dV2(mp.dm2.act_ele);
            %dE2plus = (4*squeeze(jacStruct.G1(:, :, modeIndex))*dV2plus(mp.dm1.act_ele)).*eye(2*Npix); %is eye the right format for the multiplication?
        
    
    
            %dV2minus = DM2Vminus(:, :, iProbe) - DM1Vnom;
            %dV2minus = dV2minus(mp.dm2.act_ele);
            %dE2minus = 4*squeeze(jacStruct.G1(:, :, modeIndex))*dV2minus(mp.dm1.act_ele).*eye(2*Npix);
    
        end
    
        dV = [dV, dV1; dV2];
        
        %dVminus = [dV1minus; dV2minus];
    
        %gdus = 4*G_bb*dV; % = dE % .*eye(2*Npix); %
        %dEminus = 4*G_bb*dVminus; %.*eye(2*Npix);
    
    
        %gdus = [gdus, dE]; % maybe cat better?
    
    end
       
    gdus1 = 4*G_bb*dV; %2*Nsbp*Npix x Npairs 
     
    % concatenating gdus so that its shape is (number of pixels)x(number of iterations or probes)x(number of wavelengths times two)
    %gdus_cat = cat(1, gdus); %does this makes gdus a 2d array?
    
    
    %gdus_transpose = permute(gdus, [2, 1, 3]);
    
    gdus_transpose = reshape(gdus1.', mp.est.probe.Npairs, 2*mp.Nsbp_bb, mp.Fend.corr.Npix); %each page is one pixel, one page is Npairs x 2*Nsbp
    gdus = permute(gdus_transpose, [2, 1, 3]); %each page is one pixel, one page is 2*Nsbp x Npairs
    
    % computing the pseudo-inverses of each matrix in gdus in 3 steps. First, Z is an array of square matrices and its shape is  (number of pixels)x(number of wavelengths times two)x(number of wavelengths times two)
    %gdus_transpose = permute(gdus, [2, 1, 3]);

%TEST
%     Z = pagemtimes(gdus_transpose, gdus);  %validate the order of the transpose..
%     
%     % then computing the inverse of the square matrices in Z
%     Z_inv = pageinv(Z);
%     
%     %then the pseudo-inverses
%     gdus_pinv = permute(pagemtimes(Z_inv, gdus_transpose), [2, 1, 3]);
%     
%     
%     Eest_notComplex = pagemtimes(gdus_pinv, zAll4_t); %E_hat %shape is 2*Nsbp x 1 x Npix 
%     
%     %reshape to have Npix x Nwpsbp complex matrix
%     Eest_complex1 = reshape(Eest_notComplex, [2, mp.Nsbp_bb, mp.Fend.corr.Npix]); %separates Nwpsbp in Re and Im
%     Eest_complex2 = permute(Eest_complex1, [3, 2, 1]); % Npix x Nsbp_bb x 2
%     Eest = Eest_complex2(:, :, 1) + 1i * Eest_complex2(:, :, 2); % Npix x Nsbp complex matrix
% 
    
    gdus_pinv = zeros(size(gdus));
    for i = 1:2:(2*mp.Nsbp_bb - 1) %loop over the wavelengths (step of 2 because of Re and Im)
        n = i + 1; % to take Im and Re part
        
        gdus_transpose_i = gdus_transpose(:,i:n,:);
        gdus_i = gdus(i:n,:,:);
        Z = pagemtimes(gdus_transpose_i, gdus_i);
        Z_inv = pageinv(Z);
        Zinv_gdustranspose = pagemtimes(Z_inv, gdus_transpose_i);
        gdus_pinv(i:n, :, :) = permute(Zinv_gdustranspose, [2,1,3]);
    
    end
    
    Eest_notComplex = pagemtimes(gdus_pinv, zAll4_t); %E_hat %shape is 2*Nsbp x 1 x Npix 
    
    %reshape to have Npix x Nsbp complex matrix
    Eest_complex1 = reshape(Eest_notComplex, [2, mp.Nsbp_bb, mp.Fend.corr.Npix]); %separates Nsbp in Re and Im: 2 x Nsbp x Npix
    Eest_complex2 = permute(Eest_complex1, [3, 2, 1]); % Npix x Nsbp x 2
    Eest_bb = Eest_complex2(:, :, 1) + 1i * Eest_complex2(:, :, 2); % Npix x Nsbp complex matrix


end