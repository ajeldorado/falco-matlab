



function Eest = pairwise_bb_estimation(mp, jacStruct, DM1Vplus, DM2Vplus, zAll)
    
    DM1Vnom = mp.dm1.V;
    DM2Vnom = mp.dm2.V;
    whichDM = mp.est.probe.whichDM;

    %Npix = mp.Fend.corr.Npix;
    %Npairs = mp.est.probe.Npairs;

    %mp.Nwpsbp  = mp.Nsbp; %erase this if run another EXAMPLE_try_running_falco with pm.Nwpsbp = 3


    %Where does the 96x1 comes from?? otherwise can't do G_bb*dV
    %mp.dm1.act_ele = (1:mp.dm1.NactTotal).';
    %mp.dm2.act_ele = (1:mp.dm2.NactTotal).';

    gdus = [];%Initialize delta E-fields for each probe image
    
    G_bb = bb_jacobian(mp, jacStruct, mp.dm_ind);
   
    dV1 = []; 
    dV2 = [];
    dV = [];
    for iProbe = 1:mp.est.probe.Npairs

        zAll4_t = reshape(4*zAll.', [mp.est.probe.Npairs, 1, mp.Fend.corr.Npix]);
         

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
       
    gdus = 4*G_bb*dV; %outside of the loop right?
     
    % concatenating gdus so that its shape is (number of pixels)x(number of iterations or probes)x(number of wavelengths times two)
    %gdus_cat = cat(1, gdus); %does this makes gdus a 2d array?
    
    
    %gdus_transpose = permute(gdus, [2, 1, 3]);
    
    gdus_transpose = reshape(gdus.', mp.est.probe.Npairs, 2*mp.Nsbp_bb, mp.Fend.corr.Npix);
    gdus = permute(gdus_transpose, [2, 1, 3]);
    
    % computing the pseudo-inverses of each matrix in gdus in 3 steps. First, Z is an array of square matrices and its shape is  (number of pixels)x(number of wavelengths times two)x(number of wavelengths times two)
    %gdus_transpose = permute(gdus, [2, 1, 3]);
    Z = pagemtimes(gdus_transpose, gdus);  %validate the order of the transpose..
    
    % then computing the inverse of the square matrices in Z
    Z_inv = pageinv(Z);
    
    %then the pseudo-inverses
    gdus_pinv = permute(pagemtimes(Z_inv, gdus_transpose), [2, 1, 3]);
    
    
    Eest_notComplex = pagemtimes(gdus_pinv, zAll4_t); %E_hat %shape is 2*Nsbp x 1 x Npix 
    
    %reshape to have Npix x Nwpsbp complex matrix
    Eest_complex1 = reshape(Eest_notComplex, [2, mp.Nsbp_bb, mp.Fend.corr.Npix]); %separates Nwpsbp in Re and Im
    Eest_complex2 = permute(Eest_complex1, [3, 2, 1]); % Npix x Nwpsbp x 2
    Eest = Eest_complex2(:, :, 1) + 1i * Eest_complex2(:, :, 2); % Npix x Nwpsbp complex matrix


%still need to change the shape of Eest
%2712 x 1 complex (Npix x 1)
%

end