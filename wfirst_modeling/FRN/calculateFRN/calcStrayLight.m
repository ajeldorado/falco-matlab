function strayLight = calcStrayLight(scenarioName, Kappa, det, detmpix, dQE, integTime, thruput, CGdata, CGintmpix, IWO, OWO, f_sr, colArea, spMag0, planet)
% Calculate stray light. For EB Fiducial Planet, use stray light table. For
% others, calculate.
scenario = strrep(scenarioName, ' ', '');
scenario = strrep(scenario, '-', '_');

switch scenario
    case 'EBIMGNF_Band1'
        strayLight = 0.122152202814894 * 10e-10;
    case 'EBSPECIFS_Band3'
        strayLight = 0.413420348675147 * 10e-10;
    case 'EBIMGWF_Band4'
        strayLight = 0.137846533553296 * 10e-10;
    otherwise
        lumin = calcLuminescenceBackground(detmpix, det);
        compStar = calcCompanionStarBackground(det, detmpix, CGdata, CGintmpix, integTime, thruput, dQE, IWO, OWO, f_sr, colArea, spMag0, planet);
        strayLightReq = lumin + compStar + 1.0*(lumin + compStar) + 0.1*(lumin + compStar);
    
        strayLight = Kappa * sqrt(det.ENF^2 * strayLightReq * 1000000 * det.pixSize^2 * detmpix * dQE * integTime);
end

return

function luminBackground = calcLuminescenceBackground(detmpix, det)

GCRfluxL2 = 5;
photonsPerRelEvent = 250;
lumRatePerAngle = photonsPerRelEvent / (2*pi);
luminOpticArea = pi/4 * 1^2;
luminOpticThickness = 4;
luminOpticDistance = 0.1;
sBaffling = 0.001;

omegaSignal = detmpix * det.pixSize^2 / luminOpticDistance^2;
omegaIndirect = 2*pi * sBaffling * detmpix/det.npixDet;
phRateDirect = GCRfluxL2 * lumRatePerAngle * luminOpticArea * luminOpticThickness * omegaSignal;
phRateIndirect = GCRfluxL2 * lumRatePerAngle * luminOpticArea * luminOpticThickness * omegaIndirect;

photonRate = phRateDirect + phRateIndirect;
photonRatePerPix = photonRate / detmpix;
photonRateFlux = photonRatePerPix / det.pixSize^2;

MUF = 1.5;
luminBackground = photonRateFlux * MUF/1000000;

return


function companionStar = calcCompanionStarBackground(det, detmpix, CGdata, CGintmpix, integTime, thruput, dQE, IWO, OWO, f_sr, colArea, spMag0, planet)

darkHoleMidRadius = (IWO + OWO) / 2;
PSFpeakMidRadius = CGdata.PSFpeak(find(CGdata.rlamD <= darkHoleMidRadius, 1, 'last'));

% companion star
compStarLeak = 0.0000000006;
compStarMag = 5;
compStarSpecType = 'm5v'; % hardcoded to match spreadsheet
compStarIntegFlux = eval(['spMag0.inBand.sum.', compStarSpecType]); % hardcoded for now, should be referencing a table
compStarFluxInBand = 10^(-0.4 * compStarMag) * compStarIntegFlux;

strayLightBackground = f_sr * compStarFluxInBand * compStarLeak * PSFpeakMidRadius * CGintmpix * colArea * thruput.speck * dQE * integTime;
leakageRate = (strayLightBackground / integTime) / detmpix;
leakageRateSL = leakageRate / det.pixSize^2;

MUF = 1.5;
companionStar = leakageRateSL * MUF/1000000;

return
