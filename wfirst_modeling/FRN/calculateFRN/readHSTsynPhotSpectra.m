function [spMag0Table, spMag0] = readHSTsynPhotSpectra(file_dir, filename, centerLambda, bandWidth, uc)
% Read HSTsynPhotSpectra and extract data

h_planck = uc.h_planck; % meter^2 * kg / second;
c_light  = uc.c_light;  % meter / second;

% Initialize variables.
delimiter = ',';

startRow = 2;
endRow = inf;

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(fullfile(file_dir, filename),'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Create output variable
spMag0Table = table(dataArray{1:end-1}, 'VariableNames', {'lambda_m','Eph_J','a0v','a5v','f5v','g0v','g5v','k0v','k5v','m0v','m5v'});

spMag0 = table2struct(spMag0Table,'ToScalar',true);

bandRange = find(abs(spMag0.lambda_m-centerLambda)<=0.5*bandWidth*centerLambda);

inBandFlux0.lam = spMag0.lambda_m(bandRange);
inBandFlux0.Ephot = h_planck * c_light / centerLambda;

inBandFlux0.spec.a0v = spMag0.a0v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.a5v = spMag0.a5v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.f5v = spMag0.f5v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.g0v = spMag0.g0v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.g5v = spMag0.g5v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.k0v = spMag0.k0v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.k5v = spMag0.k5v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.m0v = spMag0.m0v(bandRange) / inBandFlux0.Ephot;
inBandFlux0.spec.m5v = spMag0.m5v(bandRange) / inBandFlux0.Ephot;

deltaLambda = spMag0.lambda_m(2)-spMag0.lambda_m(1);

inBandFlux0.sum.a0v = sum(inBandFlux0.spec.a0v) * deltaLambda;
inBandFlux0.sum.a5v = sum(inBandFlux0.spec.a5v) * deltaLambda;
inBandFlux0.sum.f5v = sum(inBandFlux0.spec.f5v) * deltaLambda;
inBandFlux0.sum.g0v = sum(inBandFlux0.spec.g0v) * deltaLambda;
inBandFlux0.sum.g5v = sum(inBandFlux0.spec.g5v) * deltaLambda;
inBandFlux0.sum.k0v = sum(inBandFlux0.spec.k0v) * deltaLambda;
inBandFlux0.sum.k5v = sum(inBandFlux0.spec.k5v) * deltaLambda;
inBandFlux0.sum.m0v = sum(inBandFlux0.spec.m0v) * deltaLambda;
inBandFlux0.sum.m5v = sum(inBandFlux0.spec.m5v) * deltaLambda;

spMag0.inBand = inBandFlux0;

return
