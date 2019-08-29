function FRN = calculateFRN(CG_Directory, MUFcaseInput, NItoContrast, mode, scenarioName, detType, CGtauPol, centerLambda, bandWidth, yearsL2, SNR, frameTime, k_pp, maxIntegTime, dutyFactor, suppressOutput, planetNo, bandNo)
% Compute flux ratio noise
%
% B. Nemati and S. Miller 07-Jan-2019

CG_Directory = fullfile(CG_Directory);
uc = unitsConstants;

% Integration Time
integTime = maxIntegTime * dutyFactor;

% Telescope Diameter
Diam = 2.37 * uc.meter;

% Planet
if exist('bandNo', 'var')
    planet = readPlanetTable(planetNo, bandNo);
else
    planet = readPlanetTable(planetNo);
end
        
lam_D = centerLambda / Diam;
planetWA = (planet.Sep/lam_D) + 10^-10;

%% Read all data tables
[inputTables, inputData, QEcurves, ~] = readDataTables(CG_Directory, MUFcaseInput, mode, centerLambda, bandWidth, planetWA, uc);

%% Coronagraph Parameters
CGdata = inputTables.CG.KristTable; % coronagraph design file (Krist table)

CGdata.coreThruput = CGdata.coreThruput * (1/CGtauPol);
CGdata.PSFpeak     = CGdata.PSFpeak     * (1/CGtauPol);
CGdata.occTrans    = CGdata.occTrans    * (1/CGtauPol);
CGintSamp  = CGdata.rlamD(2)-CGdata.rlamD(1); % the intrinsic sampling in the data
CGdesignWL = interp1(CGdata.rlamD, CGdata.r_arcsec, planetWA, 'previous') * uc.arcsec * Diam ./ ...
             interp1(CGdata.rlamD, CGdata.rlamD, planetWA, 'previous');
indCG      = find(CGdata.rlamD <= (planetWA), 1, 'last'); % index of last radial slice < planet's radius
IWO = min(CGdata.rlamD);
OWO = max(CGdata.rlamD);

% Contrast Stability
CS = CStability(NItoContrast, inputData);
rawContrast = CS.rawContrast;
differentialContrast = CS.differentialContrast;

% planet core (area with PSF contour at half-max) in CG analysis intrinsic (not detector) pixels
CGintmpix = CGdata.area_sqarcsec(indCG)*uc.arcsec^2 / (CGintSamp * CGdesignWL/Diam)^2;
speckleFluxRatio = rawContrast * CGdata.PSFpeak(indCG) * CGintmpix ;

%% Instrument Parameters
% Detector
det = setDetector(detType, QEcurves, centerLambda, yearsL2, frameTime, uc);

% Throughput
t_refl = inputData.t_refl;
thruput.refltran = t_refl;
thruput.pol      = 1.0;
% Compute the throughputs, keeping in mind that:
%   occ = fpm x lyot
%   core = psf x pup x occ = psf x pup x fpm x lyot
%
%   planet  = psf pup occ ref fil pol
%   zodi    = 1   pup occ ref fil pol
%   speckle = 1    1   1  ref fil pol(or 1)
MUF_thput = det.MUF_thput; % this is in the spreadsheet as an additional MUF on the throughput
thruput.occulter = CGdata.occTrans(indCG); % FPM x Lyot (non-POLX)
thruput.core     = CGdata.coreThruput(indCG) * MUF_thput; % thruput.pupil * thruput.occulter * thruput.PSF;
thruput.PSF      = thruput.core / thruput.occulter;
% note that the planet thp has no tau_PSF because this is the per pixel thp
% not the thp into the core. The PSF shape itself inherently contains that
thruput.planet   = thruput.core     * thruput.refltran * thruput.pol;
thruput.zodi     = thruput.occulter * thruput.refltran * thruput.pol;
thruput.speck    = 1                * thruput.refltran * thruput.pol;

% Focal Plane
switch mode
    case 'IMG'
        f_sr       = 1.0;
        NyquistWL  = 500 * 1e-9;
        detmpix    = CGdata.area_sqarcsec(indCG) * uc.arcsec^2 * (centerLambda/CGdesignWL)^2 * ((2*Diam)/NyquistWL)^2;
    case 'IFS'
        specResol  = 50; % Spectral resolution is defined as R = lam / dlam
        nLenslets  = 5; % effective no. of lenslets subtended by the PSF core in the IFS
        Nspec      = bandWidth * specResol; % No. of spectral elements within the current band
        f_sr       = 1/Nspec;
        NyquistWL  = 660 * 1e-9; % Wavelength at which detector pixels Nyquist sample the PSF
        detmpix    = nLenslets * (centerLambda/NyquistWL)^2 * 2 * 2; % each lenslet spectrum covers a 2 x 2 pixel area on the IFS detector
end
pixelOnSky = NyquistWL/(2*Diam);
sampling = pixelOnSky/lam_D;
darkCurrentRate = det.darkCurrent * detmpix;
CICrate         = det.CIC * det.frameRate * detmpix;

% Telescope
obsc      = 0.32; % diameter of obscuration relative to primary
strutObsc = 0.07;
colArea   = (pi/4) * Diam^2 * (1-obsc^2) * (1-strutObsc);

%% Observational Parameters
% Star
% Load Spectra for different stellar type, normalized to Vmag = 0
spMag0 = inputData.spMag0; %#ok<NASGU>

inBandZeroMagFlux = eval(['spMag0.inBand.sum.', planet.UseSpec]); % ph/s/m^2
starFlux = inBandZeroMagFlux * 10^(-0.4*planet.V); % ph/s/m^2
absMag = planet.V-5 * log10(planet.PcDist/10);

% Reference Star
% Observing scenario is RDI with one observation of a brighter reference
% and one observation of the target. The noise is computed effective after the
% speckle subtraction. To the observation 0.
RefStarDeltaMag = 3;
BrightnessRatio = 10^(0.4*RefStarDeltaMag);
timeFractionReference = 0.20;
betaRDI = 1 / (BrightnessRatio*timeFractionReference);
v_spl = sqrt(betaRDI);
v_det = betaRDI * sqrt(timeFractionReference);
v_lzo = betaRDI * sqrt(timeFractionReference);
v_ezo = sqrt(betaRDI);

% Zodi
% Zodi brightness - assume local zodi = 22.1 mag /arcsec^2 in V band
% assume  skynoise = (local-zodi + exo-zodi)
magLocalZodi   = 23;
magExoZodi_1AU = 22.17; % exoZodi at 1 AU assumed for a star with the same absolute mag as Sun

% zodi solid angular flux in units of ph/s/m^2/arcsec^2
loZodiAngFlux = inBandZeroMagFlux * 10^(-0.4 * magLocalZodi);
exZodiAngFlux = inBandZeroMagFlux * 10^(-0.4*(absMag-uc.sunAbsMag+magExoZodi_1AU)) / (planet.A^2);
ZodiFlux      = (loZodiAngFlux + exZodiAngFlux) * CGdata.area_sqarcsec(indCG);
loZodiFlux    = loZodiAngFlux * CGdata.area_sqarcsec(indCG);
exZodiFlux    = exZodiAngFlux * CGdata.area_sqarcsec(indCG);

%% Results
% rate at which photo-electrons are generated in the signal core region pixels
plConvertRate  = f_sr * starFlux   * planet.FluxRatio   * colArea * thruput.planet * det.QE;
speckleBkgRate = f_sr * starFlux   * speckleFluxRatio   * colArea * thruput.speck  * det.QE;
zodiBkgRate    = f_sr * ZodiFlux   *                      colArea * thruput.zodi   * det.QE;
loZodiBkgRate  = f_sr * loZodiFlux *                      colArea * thruput.zodi   * det.QE;
exZodiBkgRate  = f_sr * exZodiFlux *                      colArea * thruput.zodi   * det.QE;
photoConvertedRate = (plConvertRate + zodiBkgRate + speckleBkgRate) / detmpix; % per detec pixel
phePixFrame   = photoConvertedRate / det.frameRate;

% delivered QE
chargeTransEff = trapCTE(detType, phePixFrame, yearsL2);
dQE = det.QE * det.effPC * det.effHotPix * det.effCosmic * chargeTransEff;
 
% multipliers for converting noise to flux ratio noise (FRN)
Kappa_cs = CGdata.PSFpeak(indCG) * CGintmpix / thruput.core;
Kappa = 1 / (f_sr * starFlux * colArea * thruput.planet * dQE * integTime);

% stray light
strayLight = calcStrayLight(scenarioName, Kappa, det, detmpix, dQE, integTime, thruput, CGdata, CGintmpix, IWO, OWO, f_sr, colArea, spMag0, planet);

% flux ratio noise contributions in parts per billion
nppb.planetShot  = Kappa * sqrt(det.ENF^2 * plConvertRate   * integTime ) / uc.ppb;
nppb.zodiShot    = Kappa * sqrt(det.ENF^2 * zodiBkgRate     * integTime ) / uc.ppb;
nppb.loZodiShot  = Kappa * sqrt(det.ENF^2 * loZodiBkgRate   * integTime ) / uc.ppb;
nppb.exZodiShot  = Kappa * sqrt(det.ENF^2 * exZodiBkgRate   * integTime ) / uc.ppb;
nppb.specklShot  = Kappa * sqrt(det.ENF^2 * speckleBkgRate  * integTime ) / uc.ppb;
nppb.Det_Dark    = Kappa * sqrt(det.ENF^2 * darkCurrentRate * integTime ) / uc.ppb;
nppb.Det_CIC     = Kappa * sqrt(det.ENF^2 * CICrate         * integTime ) / uc.ppb;
nppb.strayLight  = strayLight / uc.ppb;
nppb.RDI         = sqrt((v_det*nppb.Det_Dark)^2 + (v_det*nppb.Det_CIC)^2 + (v_spl*nppb.specklShot)^2 ...
                   + (v_lzo*nppb.loZodiShot)^2 + (v_ezo*nppb.exZodiShot)^2 + (v_det*nppb.strayLight)^2);
nppb.detecNoise  = sqrt(nppb.Det_Dark^2   + nppb.Det_CIC^2);
nppb.exoSysShot  = sqrt(nppb.planetShot^2 + nppb.specklShot^2 + nppb.zodiShot^2 + nppb.strayLight^2);
nppb.photometric = sqrt(nppb.exoSysShot^2 + nppb.detecNoise^2 + nppb.RDI^2);
nppb.CStab_chop  = Kappa_cs * differentialContrast / uc.ppb;
nppb.CStability  = nppb.CStab_chop / k_pp;
nppb.total       = sqrt(nppb.photometric^2 + nppb.CStability^2);

% scenario parameters
scenario.Name         = scenarioName;
scenario.CGdirectory  = CG_Directory;
scenario.mode         = mode;
scenario.detType      = detType;
scenario.CGtauPol     = CGtauPol;
scenario.centerLambda = centerLambda;
scenario.lam_D        = lam_D;
scenario.pixelOnSky   = pixelOnSky;
scenario.sampling     = sampling;
scenario.bandWidth    = bandWidth;
scenario.SNR          = SNR;
scenario.yearsL2      = yearsL2;
scenario.frameTime    = frameTime;
scenario.k_pp         = k_pp;
scenario.integTime    = integTime;
scenario.MUFcase      = MUFcaseInput;
scenario.planetWA     = planetWA;
scenario.planet       = planet;
scenario.inputTables  = inputTables;
scenario.inputData    = inputData;

if ~suppressOutput
    fprintf('CG Directory : %s\n', CG_Directory);

    fprintf('Mode         : %s\n', mode);
    fprintf('Scenario     : %s\n', scenario.Name);
    fprintf('Planet       : %s\n', planet.Name);
    fprintf('MUF Case     : %s\n', MUFcaseInput);
    fprintf('Annular Zone : %d (for PWA of %.2f)\n', inputData.annzone, planetWA);
end

FRN.scenario     = scenario;
FRN.photomDetail = nppb;
FRN.cstabDetail  = CS;
FRN.photometry   = nppb.photometric;
FRN.cstability   = nppb.CStability;
FRN.total        = nppb.total;

return

function out = setDetector(detCase, QEcurves, centerLambda, yrsL2, frameTime, uc)
% set detector parameters
switch detCase
    case 'EM_PC_SRR'
        amplifierReadNoise = 75; % e- amplifier inherent read noise
        EMgain  = 3600;
        CRtail  = 0.0323*EMgain + 133.5;
        CRflux  = 5 * 1e4; % hits per m2 per sec
        npixDet = 1024^2;
        detArea = 1.7e-4;
        CRhits  = CRflux * detArea *frameTime;

        det.npixDet     = npixDet;
        det.type        = detCase;
        det.maxYrsL2    = 63/12;
        det.doseType    = 'RDF2 95%CL 5yrs';
        det.QE          = interp1(QEcurves.lambda, QEcurves.e2v_Spec, centerLambda, 'previous');
        det.effPC       = exp(-5*amplifierReadNoise / EMgain);
        det.effHotPix   = 1 - (0.05*yrsL2) / det.maxYrsL2; %fraction of pixels that are hot by end of mission level of dosage
        det.effCosmic   = 1 - (CRtail*CRhits) / npixDet;
        det.frameRate   = 1/frameTime;
        det.ENF         = 1;
        det.dark_BE     = [1.5, 2] / uc.hour; % e-/pix/hr
        derate_dark     = 1.0; % calculating for CBE case of EB -- no derate
        derate_CIC      = 1.0; % calculating for CBE case of EB -- no derate
        det.readNoise   = 0;
        det.darkCurrent = derate_dark * (det.dark_BE(1) + (yrsL2/det.maxYrsL2) * (det.dark_BE(2) - det.dark_BE(1)));
        det.CIC         = derate_CIC * (0.000004337*EMgain + 0.0076);
        det.MUF_thput   = 1/1.1;
        det.pixSize     = 1.3e-5;
    case 'EM_PC_REQ'
        %%%%% THIS IS NOT CORRECT YET %%%%%
        error('EM_PC_REQ detector not yet supported')
%         amplifierReadNoise = 75; % e- amplifier inherent read noise
%         EMgain  = 3600;
%         CRtail  = 0.0323*EMgain + 133.5;
%         CRflux  = 5 * 1e4; % hits per m2 per sec
%         npixDet = 1024^2;
%         detArea = 1.7e-4;
%         CRhits  = CRflux * detArea *frameTime;
% 
%         det.npixDet     = npixDet;
%         det.type        = detCase;
%         det.maxYrsL2    = 63/12;
%         det.doseType    = 'RDF2 95%CL 5yrs';
%         det.QE          = interp1(QEcurves.lambda, QEcurves.e2v_Spec, centerLambda, 'previous');
%         det.effPC       = exp(-5*amplifierReadNoise / EMgain);
%         det.effHotPix   = 1 - (0.05*yrsL2) / det.maxYrsL2; %fraction of pixels that are hot by end of mission level of dosage
%         det.effCosmic   = 1 - (CRtail*CRhits) / npixDet;
%         det.frameRate   = 1/frameTime;
%         det.ENF         = 1;
%         det.dark_BE     = [1.5, 2] / uc.hour; % e-/pix/hr
%         derate_dark     = 1.0; % calculating for CBE case of EB -- no derate
%         derate_CIC      = 1.0; % calculating for CBE case of EB -- no derate
%         det.readNoise   = 0;
%         det.darkCurrent = derate_dark * (det.dark_BE(1) + (yrsL2/det.maxYrsL2) * (det.dark_BE(2) - det.dark_BE(1)));
%         det.CIC         = derate_CIC * (0.000004337*EMgain + 0.0076);
%         det.MUF_thput   = 1/1.1;
%         det.pixSize     = 1.3e-5;
    case 'EM_PC_MC_CBE'
        amplifierReadNoise = 40; % e- amplifier inherent read noise
        EMgain  = 1900;
        CRtail  = 0.0323*EMgain + 133.5;
        CRflux  = 5 * 1e4; % hits per m2 per sec
        npixDet = 1024^2;
        detArea = 1.7e-4;
        CRhits  = CRflux * detArea *frameTime;

        det.npixDet     = npixDet;
        det.type        = detCase;
        det.maxYrsL2    = 63/12;
        det.doseType    = 'RDF2 95%CL 5yrs';
        det.QE          = interp1(QEcurves.lambda, QEcurves.e2v_Spec, centerLambda, 'previous');
        det.effPC       = exp(-5*amplifierReadNoise / EMgain);
        det.effHotPix   = 1 - (0.05*yrsL2) / det.maxYrsL2; %fraction of pixels that are hot by end of mission level of dosage
        det.effCosmic   = 1 - (CRtail*CRhits) / npixDet;
        det.frameRate   = 1/frameTime;
        det.ENF         = 1;
        det.dark_BE     = [1.5, 2] / uc.hour; % e-/pix/hr 
        derate_dark     = 1.0; % calculating for CBE case of EB -- no derate
        derate_CIC      = 1.0; % calculating for CBE case of EB -- no derate
        det.readNoise   = 0;
        det.darkCurrent = derate_dark * (det.dark_BE(1) + (yrsL2/det.maxYrsL2) * (det.dark_BE(2) - det.dark_BE(1)));
        det.CIC         = derate_CIC * (0.000004337*EMgain + 0.0076);
        det.MUF_thput   = 1;
        det.pixSize     = 1.3e-5;
end
out = det;

return

function out = trapCTE(detCase, pheOneFrame, yearsL2)
% Based on data from Patrick Morrissey for 95%CL, RDF=2, 5yrs
switch detCase
    case 'EM_PC_SRR'
        testYrsL2 = 5.25;      
        derate_CTE = 0.93;

        radDoseFraction = yearsL2 / testYrsL2;
        dqeFluxSlope =	3.24;
        dqeKnee	= 0.858;
        dqeKneeFlux	= 0.089;

        out = derate_CTE * max(0,min(1+radDoseFraction*(dqeKnee-1),...
            1+radDoseFraction*(dqeKnee-1)+radDoseFraction*dqeFluxSlope*(pheOneFrame-dqeKneeFlux)));
    case 'EM_PC_REQ'
        %%%%% THIS IS NOT CORRECT YET %%%%%
        error('EM_PC_REQ detector not yet supported')
%         testYrsL2 = 5.25;
%         derate_CTE = 0.83;
% 
%         radDoseFraction = yearsL2 / testYrsL2;
%         dqeFluxSlope =	3.24;
%         dqeKnee	= 0.858;
%         dqeKneeFlux	= 0.089;
% 
%         out = derate_CTE * max(0,min(1+radDoseFraction*(dqeKnee-1),...
%             1+radDoseFraction*(dqeKnee-1)+radDoseFraction*dqeFluxSlope*(pheOneFrame-dqeKneeFlux)));
    case 'EM_PC_MC_CBE'
        testYrsL2 = 5.25;      
        derate_CTE = 0.5;

        radDoseFraction = yearsL2 / testYrsL2;
        dqeFluxSlope =	3.24;
        dqeKnee	= 0.858;
        dqeKneeFlux	= 0.089;

        out = derate_CTE * (1+max(0,min(1+radDoseFraction*(dqeKnee-1),...
            1+radDoseFraction*(dqeKnee-1)+radDoseFraction*dqeFluxSlope*(pheOneFrame-dqeKneeFlux))));
end

return
