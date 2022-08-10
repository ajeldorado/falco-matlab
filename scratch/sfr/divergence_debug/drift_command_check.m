data_path = ['C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial54'];

drift_it1 = fitsread(fullfile(data_path,'EKF_Series8_Trial54EKF_Series8_Trial54drift_command_it1.fits'));

drift_it231 = fitsread(fullfile(data_path,'EKF_Series8_Trial54EKF_Series8_Trial54drift_command_it231.fits'));


drift_it1 = drift_it1(:,:,1);

drift_it231 = drift_it231(:,:,1);



figure
subplot(1,3,1)
imagesc(drift_it1)
colorbar
title('iteration1')

subplot(1,3,2)
imagesc(drift_it231)
colorbar
title('iteration231')

subplot(1,3,3)
imagesc(drift_it231 - drift_it1)
colorbar
title('iteration231 - iteration1')
%%
data_path = ['C:\Users\sfr\OneDrive - Princeton University\Documents\JPL\testbed_data\dzm_runs\EKF_Series8_Trial55'];

drift_it1 = fitsread(fullfile(data_path,'drift_command_it1.fits'));

drift_it5 = fitsread(fullfile(data_path,'drift_command_it5.fits'));


drift_it1 = drift_it1(:,:,1);

drift_it5 = drift_it5(:,:,1);



figure
subplot(1,3,1)
imagesc(drift_it1)
colorbar
title('iteration1')

subplot(1,3,2)
imagesc(drift_it5)
colorbar
title('iteration5')

subplot(1,3,3)
imagesc(drift_it5 - drift_it1)
colorbar
title('iteration5 - iteration1')