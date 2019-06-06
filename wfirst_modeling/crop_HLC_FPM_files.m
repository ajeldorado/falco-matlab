
clear all;

mp.full.phaseb_dir = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/'; % mask design data path


cd([mp.full.phaseb_dir 'hlc_20190210'])


prefix = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/hlc_20190210/run461_nro_';


mp.lambda0 = 575e-9;
mp.fracBW = 0.10;
mp.Nsbp = 19;

lam_occ = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp)*mp.lambda0



fpm_axis = 'p';
% fpm_axis = 's';

transVec = zeros(mp.Nsbp,1);

for wlam = 1:mp.Nsbp

    fn_p_r = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'real.fits'];
    fn_p_i = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'imag.fits'];

    fn_p_r_rot = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'real_rotated.fits'];
    fn_p_i_rot = [prefix  'occ_lam' num2str(lam_occ(wlam),12) 'theta6.69pol'   fpm_axis   '_' 'imag_rotated.fits'];


    Ncrop = 40;

    a = fitsread(fn_p_r);
    b = fitsread(fn_p_i);
    c = fitsread(fn_p_r_rot);
    d = fitsread(fn_p_i_rot);
    
    transVec(wlam) = a(1,1) + 1j*b(1,1);
    
    %--Undo the arbitrary phase offsets
    ab = exp(-1j*angle(transVec(wlam)))*(a+1j*b);
    a = real(ab);
    b = imag(ab);
    cd = exp(-1j*angle(transVec(wlam)))*(c+1j*d);
    c = real(cd);
    d = imag(cd);
    
    %--Crop
    acrop = padOrCropEven(a,Ncrop);
    bcrop = padOrCropEven(b,Ncrop);
    ccrop = padOrCropEven(c,Ncrop);
    dcrop = padOrCropEven(d,Ncrop);

    fprintf('Phase = %.4g\n',angle(acrop(1,1)+1i*bcrop(1,1)));
    
    fitswrite(acrop,[fn_p_r(1:end-5),'_crop.fits']);
    fitswrite(bcrop,[fn_p_i(1:end-5),'_crop.fits']);
    fitswrite(ccrop,[fn_p_r_rot(1:end-5),'_crop.fits']);
    fitswrite(dcrop,[fn_p_i_rot(1:end-5),'_crop.fits']);

    figure(1); imagesc(a); axis xy equal tight; colorbar; drawnow;
    figure(2); imagesc(angle(acrop+1i*bcrop)); axis xy equal tight; colorbar;drawnow;
    
    pause(0.5);
    
end

%--Plot the arbitrary phase shift for each HLC occulter
figure(11); plot(1e9*lam_occ,angle(transVec)/(2*pi),'-bo','Linewidth',2,'Markersize',10);  
title('Phase Discrepancies among HLC Occulter Files','Fontsize',16,'Interpreter','Latex');
xlabel('Wavelength (nm)','Fontsize',20,'Interpreter','Latex');
ylabel('Occ. Substrate Phase (waves)','Fontsize',20,'Interpreter','Latex');
set(gca,'Fontsize',20);
set(gcf,'Color','w');
grid on;
xlim([min(1e9*lam_occ),max(1e9*lam_occ)])
drawnow;
% export_fig('/Users/ajriggs/Downloads/fig_HLC_occ_phase.png','-dpng','-r200');


%% Compare the HLC's pupil and LS compared to the ones generated in FALCO

clear all;

%--Library locations. FALCO and PROPER are required. CVX is optional.
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
%%--Add to the MATLAB Path
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path


mp.full.phaseb_dir = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/'; % mask design data path

cd([mp.full.phaseb_dir 'hlc_20190210'])

prefix = '/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/hlc_20190210/run461_nro_';

pupil = fitsread('run461_nro_pupil.fits');
pupil = padOrCropEven(pupil,310);
LS = fitsread('run461_nro_lyot.fits');
LS = padOrCropEven(LS,310);

figure(11); imagesc(pupil); axis xy equal tight; colorbar; drawnow;
figure(12); imagesc(LS); axis xy equal tight; colorbar; drawnow;

Nbeam = 309;
pupil2 = falco_gen_pupil_WFIRST_CGI_180718(Nbeam,'pixel');
% pupil2(2:end,2:end) = fliplr(pupil2(2:end,2:end));

changes.ID = 0.50;
changes.OD = 0.80;
changes.wStrut = 3.6/100;
changes.flagRot180 = true;
LS2 = falco_gen_pupil_WFIRST_CGI_180718(Nbeam,'pixel',changes);


figure(21); imagesc(pupil2); axis xy equal tight; colorbar; drawnow;
figure(22); imagesc(LS2); axis xy equal tight; colorbar; drawnow;

figure(31); imagesc(pupil-pupil2); axis xy equal tight; colorbar; drawnow;
figure(32); imagesc(LS-LS2); axis xy equal tight; colorbar; drawnow;
