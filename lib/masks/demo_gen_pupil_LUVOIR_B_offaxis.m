% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
clear; 
% addpath('segMirrorFunctions');
%%--Add to the MATLAB Path
mp.path.falco = '~/Repos/falco-matlab/';  %--Location of FALCO
mp.path.proper = '~/Documents/MATLAB/PROPER/'; %--Location of the MATLAB PROPER library
addpath(genpath(mp.path.falco)) %--Add FALCO library to MATLAB path
addpath(genpath(mp.path.proper)) %--Add PROPER library to MATLAB path

apDia = 500;
inputs.Nbeam = apDia/0.925; % number of points across the pupil diameter
inputs.wGap = 6e-3/7.989*apDia; % samples
inputs.numRings = 4;% Number of rings in hexagonally segmented mirror 
inputs.Npad = 2^(nextpow2(apDia));
inputs.ID = 0; % central obscuration radius 
inputs.OD = 1; % pupil outer diameter, can be < 1
inputs.Nstrut = 0;% Number of struts 
inputs.angStrut = [];%Angles of the struts (deg)
inputs.wStrut = []; % Width of the struts (fraction of pupil diam.)

missingSegments = ones(1,hexSegMirror_numSegments( inputs.numRings ));
for index = 0:5
    missingSegments(38+index*4) = 0;
end
inputs.missingSegments = missingSegments;

pupilHG = falco_gen_pupil_customHex(inputs);

[rows, cols] = size(pupilHG);
[X, Y] = meshgrid(-rows/2:rows/2-1);
[~, RHO] = cart2pol(X, Y);

figure(1); imagesc(pupilHG+(RHO<=apDia/2)); axis xy equal tight;
figure(2); imagesc(pupilHG); axis xy equal tight; title('Hypergaussian', 'Fontsize', 20);


%% Generate with PROPER
clear inputs
inputs.Nbeam = apDia*0.96075;  % factor makes the beam size match the hypergaussian approach
% inputs.Nbeam = apDia*0.9607;
pupilPROPER = falco_gen_pupil_LUVOIR_B(inputs);
pupilPROPER = pad_crop(pupilPROPER, size(pupilHG));

figure(3); imagesc(pupilPROPER); axis xy equal tight; title('PROPER', 'Fontsize', 20);

diff = pupilPROPER - pupilHG;
fprintf('Summed squared difference = %.4g\n', sum(diff(:).^2));
figure(4); imagesc(diff); axis xy equal tight; colorbar;


%% Tune the Hypergaussian exponent by comparing to the PROPER-generated aperture
NpupVec = 200:100:1000;
hgExpVec = [44, 64, 84, 110, 130, 160, 194, 224, 250]; % found empirically
figure(5); plot(NpupVec, hgExpVec, '-bd', 'Linewidth', 3);

x = NpupVec(1:end);
y = hgExpVec(1:end);

f0 = x;
f1 = ones(size(x));

A = [f0; f1];

coefVec = A.'\y.';
a0 = coefVec(1);
a1 = coefVec(2);

figure(5); plot(NpupVec, hgExpVec, '-bd', x, a0*x + a1, '--ro', 'Linewidth', 3);
xlabel('Nbeam', 'Fontsize', 20);
ylabel('Hypergaussian Exponent', 'Fontsize', 20);
legend('Measured', 'Best Linear Fit', 'Location', 'best');
set(gca, 'Fontsize', 20);
drawnow;

