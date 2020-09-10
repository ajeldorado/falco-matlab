% Produce both the 2-D representation and x and y points for a vertical bowtie mask.
%
% Parameters
% ----------
% Nbeam : float 
%    number of pixel widths across the diameter of the beam
% ID : float
%    inner diameter of bowtie. Units of beam diameters.
% OD : float
%    outer diameter of bowtie. Units of beam diameters.
% rocFillet : float
%    radius of curvature of fillets. Units of beam diameters.
% angDeg : float
%    Opening angle of each half of the bowtie. Units of degrees.
% clockDeg : float
%    Clocking angle of the bowtie mask from the x-axis going
%    counter-clockwise. At zero, the bowtie openings are centered on the
%    x-axis. Units of degrees.
% upsampleFactor : int
%    Number of sub-samples across each edge pixel to use when generating
%    the 2-D representation of the mask with anti-aliased edges. 
% DbeamUM : float
%    Beam diameter. Units of microns.
% stepSize : float
%    step size between points along the arcs to approximate a circle. Units
%    of microns.
% xc : float
%    x center coordinate of mask
% yc : float
%    y center coordinate of mask
% centering : {'pixel', 'interpixel'}
%    centering of the beam in the array. They are equivalent for an
%    odd-sized array.
%
% Notes
% -----
%
% Written by A.J. Riggs at JPL on June 26, 2020.


function [mask, xyAllCell] = falco_gen_rounded_bowtie_vertical_vector(Nbeam, ID, OD, rocFillet, angDeg, clockDeg, upsampleFactor, DbeamUM, stepSize, xc, yc, centering)

flagPlot = false;
flagPlotFinal = false;

rocFilletUM = rocFillet*DbeamUM;%170.; %--width of the dropout region [microns]

% Number of samples across in square output array
maxOffset = Nbeam*max(abs([xc/DbeamUM, yc/DbeamUM])); % max offset [pixels]
switch centering
    case 'pixel'
        Narray = 2*ceil(0.5*(Nbeam + 2*maxOffset + 1));
    case 'interpixel'
        Narray = 2*ceil(0.5*(Nbeam + 2*maxOffset));
    otherwise
        error("Centering must be 'pixel' or 'interpixel'")
end

nPointsCircle = round(2*pi*DbeamUM/2./stepSize); % Number of points to make a full circle at the OD
nPointsFillet = round(2*pi*rocFilletUM/stepSize);
% nPointsCircle = 200; % for CSV files to make CAD model
% nPointsFillet = 30;

pupilCube = zeros(Narray, Narray, 4);
clockRad = pi/180*clockDeg;
angRad = pi/180*angDeg;

xShear = xc/DbeamUM; % pupil diameters
yShear = yc/DbeamUM; % pupil diameters
magFac = 1.;

switch centering
    case 'interpixel'
        x = (-(Narray-1)/2:(Narray-1)/2)/Nbeam;
    otherwise
        x = (-Narray/2:(Narray/2-1))/Nbeam;
end
y = x;
x = x - xShear;
y = y - yShear;
[X, Y] = meshgrid(x,y);
dx = x(2) - x(1);
radius = 0.5;

[THETAS, RHOS] = cart2pol(X, Y);
% THETAS = atan2((X-xShear), (Y-yShear));
% RHOS = sqrt((X-xShear).^2 + (Y-yShear).^2);

slopeA = tan(angRad/2. - clockRad);
slopeB = tan(angRad/2. + clockRad);
%--Vertical openings
mask = double(RHOS>=ID/2. & RHOS <= OD/2. & ((X <= slopeA*Y & X >= -slopeB*Y & Y > 0) | (X >= slopeA*Y & X <= -slopeB*Y & Y < 0)));

if(flagPlot)
    figure(11); imagesc(x, x, mask); axis xy equal tight; colorbar; drawnow;
end

kernel = ones(3);
kernel = kernel/sum(kernel(:));
mask2 = conv2(mask, kernel);
mask2 = pad_crop(mask2, size(mask));
grayInds = find(mask2 > 0 & mask2 < 1);


dxUp = dx/upsampleFactor;
xUp = (-(upsampleFactor-1)/2:(upsampleFactor-1)/2)*dxUp;
[Xup0, Yup0] = meshgrid(xUp);

subpixel = zeros(upsampleFactor, upsampleFactor);

pupil = mask;%zeros(size(mask));

for iInterior = 1:length(grayInds)

    subpixel = 0*subpixel;

    xCenter = X(grayInds(iInterior));
    yCenter = Y(grayInds(iInterior));
    Xup = Xup0 + xCenter;
    Yup = Yup0 + yCenter;
    RHOSup = sqrt((Xup).^2 + (Yup).^2);

    subpixel(RHOSup>=ID/2. & RHOSup <= OD/2. & ((Xup <= slopeA*Yup & Xup >= -slopeB*Yup & Yup > 0) | (Xup >= slopeA*Yup & Xup <= -slopeB*Yup & Yup < 0))) = 1;
    pixelValue = sum(subpixel(:))/upsampleFactor^2;
    pupil(grayInds(iInterior)) = pixelValue;

end
if(flagPlot)
    hold off
end

grayMap = zeros(size(mask));
grayMap(grayInds) = 1;

if(flagPlot)
    figure(1); imagesc(x, x, grayMap); axis xy equal tight; colorbar; drawnow;

    figure(2); imagesc(x, x, mask); axis xy equal tight; colorbar; drawnow;
    figure(3); imagesc(x, x, mask2); axis xy equal tight; colorbar; drawnow;

    figure(4); imagesc(x, x, pupil); axis xy equal tight; colorbar; colormap parula; drawnow;
    figure(5); imagesc(x, x, pupil-fliplr(pupil)); axis xy equal tight; colorbar; drawnow;
end

%% Generate normalized coordinates of the zone to etch including fillets
% Coordinates are normalized to the unmasked pupil diameter along its major axis.

nStrut = 4;

strutWidthVec = zeros(1, nStrut); %wStrut*ones(1, nStrut);
yOffsetVec = zeros(nStrut, 1);
strutAngleVec = 90 + [angDeg/2.+clockDeg, -angDeg/2.+clockDeg, 180+angDeg/2.+clockDeg, 180-angDeg/2+clockDeg]; %--Easier to start -90 degrees off, then rotate later.
strutSlopeVec = tand(strutAngleVec);

Rin = ID/2.0;
Rout = OD/2.0;

RinBias = Rin + rocFillet;
RoutBias = Rout - rocFillet;

th = linspace(0, 360, 1000);%10000);



thROC = linspace(0, 360, 361);
xROC = rocFillet*cosd(thROC);
yROC = rocFillet*sind(thROC);

for istrut = [1, 3] %nStrut:-1:1 % 1:nStrut % 
    
    if(flagPlot)
        figure(33); 
        plot(ID/2*cosd(th), ID/2*sind(th), '--b', OD/2*cosd(th), OD/2*sind(th), '--b', 'Linewidth', 0.5);
        axis xy equal tight;
        xlim(0.55*[-1, 1]); ylim(0.55*[-1, 1])
        hold on;
    end
    xs = -0.5:0.01:0.5; %--for plotting only
    
    %--Strut A, right edge (when going from ID to OD)
    m1 = strutSlopeVec(istrut);
    b1 = yOffsetVec(istrut) - strutWidthVec(istrut)/2.0/cosd(strutAngleVec(istrut));
    %--Strut A, inner circle intersection
    x1ia = (-b1*m1 - sqrt(Rin*Rin*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    x1ib = (-b1*m1 + sqrt(Rin*Rin*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    if istrut == 3
        x1i = max([x1ia, x1ib]);
    else
        x1i = min([x1ia, x1ib]);
    end
    y1i = m1*x1i + b1;
    %--Strut A, outer circle intersection
    x1oa = (-b1*m1 - sqrt(Rout*Rout*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    x1ob = (-b1*m1 + sqrt(Rout*Rout*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    if istrut == 3
        x1o = max([x1oa, x1ob]);
    else
        x1o = min([x1oa, x1ob]);
    end
    y1o = m1*x1o + b1;
    if(flagPlot)
        plot([x1i, x1o], [y1i, y1o], 'rx', 'MarkerSize', 20)
    end
    %--Strut B, left edge (when going from ID to OD)
    if istrut == 2
        istrutB = 1;
    else
        istrutB = istrut+1;
    end
    m2 = strutSlopeVec(istrutB);
    bB = yOffsetVec(istrutB) + strutWidthVec(istrutB)/2.0/cosd(strutAngleVec(istrutB));
    %--Strut B, inner circle intersection
    xm = (-bB*m2 - sqrt(Rin*Rin*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    xp = (-bB*m2 + sqrt(Rin*Rin*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    if istrut == 1 || istrut == 2 || istrut == 6
        x2i = max([xm, xp]);
    else
        x2i = min([xm, xp]);
    end
    y2i = m2*x2i + bB;
    %--Strut B, outer circle intersection
    xm = (-bB*m2 - sqrt(Rout*Rout*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    xp = (-bB*m2 + sqrt(Rout*Rout*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    if istrut == 1 || istrut == 2 || istrut == 6
        x2o = max([xm, xp]);
    else
        x2o = min([xm, xp]);
    end
    y2o = m2*x2o + bB;
    if(flagPlot)
        plot([x2i, x2o], [y2i, y2o], 'rx', 'MarkerSize', 8)
    end
    %--Dropout biased from Strut B, left edge (when going from ID to OD)
    if istrut == 2
        istrutB = 1;
    else
        istrutB = istrut+1;
    end
    m2 = strutSlopeVec(istrutB);
    bB = yOffsetVec(istrutB) + (strutWidthVec(istrutB)/2.0 + rocFillet)/cosd(strutAngleVec(istrutB));
    %--Dropout biased from Strut B, inner circle intersection
    xm = (-bB*m2 - sqrt(RinBias*RinBias*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    xp = (-bB*m2 + sqrt(RinBias*RinBias*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    if istrut == 1
        xfBi = max([xm, xp]);
    else
        xfBi = min([xm, xp]);
    end
    yfBi = m2*xfBi + bB;
    %--Dropout biased from Strut B, outer circle intersection
    xm = (-bB*m2 - sqrt(RoutBias*RoutBias*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    xp = (-bB*m2 + sqrt(RoutBias*RoutBias*(m2*m2 + 1) -bB*bB))/(m2*m2+1);
    if istrut == 1 
        xfBo = max([xm, xp]);
    else
        xfBo = min([xm, xp]);
    end
    yfBo = m2*xfBo + bB;
    if(flagPlot)
        plot([xfBi, xfBo], [yfBi, yfBo], 'gx', 'MarkerSize', 8);
    end
    
    %--Fillet center biased from Strut A, right edge (when going from ID to OD)
    b1 = yOffsetVec(istrut) - (strutWidthVec(istrut)/2.0 + rocFillet)/cosd(strutAngleVec(istrut));
    %--Fillet center biased from Strut A, inner circle intersection
    x1ia = (-b1*m1 - sqrt(RinBias*RinBias*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    x1ib = (-b1*m1 + sqrt(RinBias*RinBias*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    if istrut == 3
        xfAi = max([x1ia, x1ib]);
    else
        xfAi = min([x1ia, x1ib]);
    end
    yfAi = m1*xfAi + b1; 
    %--Fillet center biased from Strut A, outer circle intersection
    x1oa = (-b1*m1 - sqrt(RoutBias*RoutBias*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    x1ob = (-b1*m1 + sqrt(RoutBias*RoutBias*(m1*m1 + 1) -b1*b1))/(m1*m1+1);
    if istrut == 3
        xfAo = max([x1oa, x1ob]);
    else
        xfAo = min([x1oa, x1ob]);
    end
    yfAo = m1*xfAo + b1;
    if(flagPlot)
        plot([xfAi, xfAo], [yfAi, yfAo], 'gx', 'MarkerSize', 8)
    end
%     plot(xROC+x1iBias, yROC+y1iBias, '--k')
    
    %  Intersection Points go clockwise with a before b.
    
    %--Strut A, Outer Fillet 1st Intersection Point (along circle)
    xfAoa = Rout*cos(atan2(yfAo, xfAo));
    yfAoa = Rout*sin(atan2(yfAo, xfAo));
    %--Strut A, Outer Fillet 2nd Intersection Point (along strut)
    b1 = yOffsetVec(istrut) - (strutWidthVec(istrut)/2.0)/cosd(strutAngleVec(istrut)); % along strut
    bf1o = yfAo + xfAo/m1; % line perpendicular to strut and intersecting fillet center
    xfAob = (bf1o-b1) / (m1 + 1/m1);
    yfAob = m1*xfAob + b1;
    [ang1, ang2] = dualatan2ccw(xfAoa-xfAo, yfAoa-yfAo, xfAob-xfAo, yfAob-yfAo); % [radians]
    thfAo = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletAo = rocFillet*cos(thfAo) + xfAo;
    yFilletAo = rocFillet*sin(thfAo) + yfAo;
    
    if(flagPlot)
        plot([xfAoa, xfAob], [yfAoa, yfAob], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Gold'), 'MarkerSize', 8);
        plot(xFilletAo, yFilletAo, '--k');
    end
    
    %--Strut A, Inner Fillet 1st Intersection Point (along strut)
    b1 = yOffsetVec(istrut) - (strutWidthVec(istrut)/2.0)/cosd(strutAngleVec(istrut)); % along strut
    bf1i = yfAi + xfAi/m1; % line perpendicular to strut and intersecting fillet center
    xfAia = (bf1i-b1) / (m1 + 1/m1);
    yfAia = m1*xfAia + b1;
    %--Strut A, Inner Fillet 2nd Intersection Point (along circle)
    xfAib = Rin*cos(atan2(yfAi, xfAi));
    yfAib = Rin*sin(atan2(yfAi, xfAi));
    [ang1, ang2] = dualatan2ccw(xfAia-xfAi, yfAia-yfAi, xfAib-xfAi, yfAib-yfAi); % [radians]
    thfAi = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletAi = rocFillet*cos(thfAi) + xfAi;
    yFilletAi = rocFillet*sin(thfAi) + yfAi;
    
    if(flagPlot)
        plot([xfAia, xfAib], [yfAia, yfAib], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Magenta'), 'MarkerSize', 8);
        plot(xFilletAi, yFilletAi, '--k');
    end
    
    %--Strut B, Inner Fillet 1st Intersection Point (along circle)
    xfBia = Rin*cos(atan2(yfBi, xfBi));
    yfBia = Rin*sin(atan2(yfBi, xfBi));
    %--Strut B, Inner Fillet 2nd Intersection Point (along strut)
    bB = yOffsetVec(istrutB) + (strutWidthVec(istrutB)/2.0)/cosd(strutAngleVec(istrutB)); % along strut
    bfB = yfBi + xfBi/m2; % line perpendicular to strut and intersecting fillet center
    xfBib = (bfB-bB) / (m2 + 1/m2);
    yfBib = m2*xfBib + bB;

    [ang1, ang2] = dualatan2ccw(xfBia-xfBi, yfBia-yfBi, xfBib-xfBi, yfBib-yfBi); % [radians]
    thfBi = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletBi = rocFillet*cos(thfBi) + xfBi;
    yFilletBi = rocFillet*sin(thfBi) + yfBi;
    
    if(flagPlot)
        plot([xfBia, xfBib], [yfBia, yfBib], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Magenta'), 'MarkerSize', 8);
        plot(xFilletBi, yFilletBi, '-k'); 
    end
   
    %--Strut B, Outer Fillet 1st Intersection Point (along strut)
    bB = yOffsetVec(istrutB) + (strutWidthVec(istrutB)/2.0)/cosd(strutAngleVec(istrutB)); % along strut
    bfB = yfBo + xfBo/m2; % line perpendicular to strut and intersecting fillet center
    xfBoa = (bfB-bB) / (m2 + 1/m2);
    yfBoa = m2*xfBoa + bB;
    %--Strut B, Outer Fillet 2nd Intersection Point (along circle)
    xfBob = Rout*cos(atan2(yfBo, xfBo));
    yfBob = Rout*sin(atan2(yfBo, xfBo));
    
    [ang1, ang2] = dualatan2ccw(xfBoa-xfBo, yfBoa-yfBo, xfBob-xfBo, yfBob-yfBo); % [radians]
    thfBo = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletBo = rocFillet*cos(thfBo) + xfBo;
    yFilletBo = rocFillet*sin(thfBo) + yfBo;
    
    if(flagPlot)
        plot([xfBoa, xfBob], [yfBoa, yfBob], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Indigo'), 'MarkerSize', 8);
        plot(xFilletBo, yFilletBo, '-k');     
    end
    
    %--Inner circle angles (CW)
    [thia, thib] = dualatan2cw(xfAib, yfAib, xfBia, yfBia); % [radians]
    thi = linspace(thia, thib, round(nPointsCircle*abs(thib-thia)/(2*pi))); % [radians]

    %--Outer circle angles (CCW)
    [thoa, thob] = dualatan2ccw(xfBob, yfBob, xfAoa, yfAoa); % [radians]
    tho = linspace(thoa, thob, round(nPointsCircle*abs(thob-thoa)/(2*pi))); % [radians]
    
    xArcInner = Rin*cos(thi);
    yArcInner = Rin*sin(thi);
    
    xArcOuter = Rout*cos(tho);
    yArcOuter = Rout*sin(tho);
    
    xAll = [xFilletAo, xFilletAi, xArcInner, xFilletBi, xFilletBo, xArcOuter];
    yAll = [yFilletAo, yFilletAi, yArcInner, yFilletBi, yFilletBo, yArcOuter];
    
    if istrut == 3
        xyAllCell{2} = [xAll; yAll];
    else
        xyAllCell{istrut} = [xAll; yAll];
    end
    if(flagPlot)
        plot(xAll, yAll, '-b','Linewidth', 1); drawnow;
    end
    

    %% Define the fillet regions as binary masks
    %--Sideways openings

    xFilletCenterVec = [xfAo, xfAi, xfBi, xfBo];
    yFilletCenterVec = [yfAo, yfAi, yfBi, yfBo];
    mSecantVec = [(yfAoa-yfAob)/(xfAoa-xfAob), (yfAia-yfAib)/(xfAia-xfAib), (yfBib-yfBia)/(xfBib-xfBia), (yfBob-yfBoa)/(xfBob-xfBoa)];
    bSecantVec = [yfAoa - xfAoa*mSecantVec(1), yfAia - xfAia*mSecantVec(2), yfBia - xfBia*mSecantVec(3), yfBoa - xfBoa*mSecantVec(4)]; 
    xDeltaVec = [(xfAoa-xfAob), (xfAia-xfAib), (xfBia-xfBib), (xfBoa-xfBob)];
    
    nFillet = length(xFilletCenterVec);

    maskFillet = zeros(size(X));
    
    for iFillet = 1:nFillet
        Xf = X - xFilletCenterVec(iFillet);
        Yf = Y - yFilletCenterVec(iFillet);
        [THETASf, RHOSf] = cart2pol(Xf, Yf);

        if xDeltaVec(iFillet) > 0
            pmf = 1;
        else
            pmf = -1;
        end
        
        if istrut == 1
            pms = 1;
        elseif istrut == 3
            pms = -1;
        end
        
        maskFillet = maskFillet | ( (pms*Y > 0) & (RHOSf.^2>=rocFillet.^2) & (pmf*Y >= pmf*(mSecantVec(iFillet)*X + bSecantVec(iFillet))) &...
            (RHOS>=ID/2. & RHOS <= OD/2. & ((X <= slopeA*Y & X >= -slopeB*Y & Y > 0) | (X >= slopeA*Y & X <= -slopeB*Y & Y < 0))));        
    end

    maskFillet = double(maskFillet);
    maskFillet2 = conv2(maskFillet, kernel);
    maskFillet2 = pad_crop(maskFillet2, size(maskFillet));
    grayInds = find(maskFillet2 > 0);

    if(flagPlot)
        figure(6); imagesc(x, x, maskFillet); axis xy equal tight; colorbar; colormap parula; drawnow;
        figure(7); imagesc(x, x, maskFillet2); axis xy equal tight; colorbar; colormap parula; drawnow;
    end


    %% Make grayscale representation of the fillet region

    subpixel = zeros(upsampleFactor, upsampleFactor);

    pupil2 = zeros(size(pupil));

    for iFillet = 1:nFillet

        Xf = X - xFilletCenterVec(iFillet);
        Yf = Y - yFilletCenterVec(iFillet);
        [THETASf, RHOSf] = cart2pol(Xf, Yf);

        if xDeltaVec(iFillet) > 0
            pmf = 1;
        else
            pmf = -1;
        end
        if istrut == 1
            pms = 1;
        elseif istrut == 3
            pms = -1;
        end
        maskFillet = ( (pms*Y > 0) & (RHOSf.^2>=rocFillet.^2) & (pmf*Y >= pmf*(mSecantVec(iFillet)*X + bSecantVec(iFillet))) &...
            (RHOS>=ID/2. & RHOS <= OD/2. & ((X <= slopeA*Y & X >= -slopeB*Y & Y > 0) | (X >= slopeA*Y & X <= -slopeB*Y & Y < 0))));
        maskFillet = double(maskFillet);
        
        maskFillet = double(maskFillet);
        maskFillet2 = conv2(maskFillet, kernel);
        maskFillet2 = pad_crop(maskFillet2, size(maskFillet));
        grayInds = find(maskFillet2 > 0);
        
        if(flagPlot)
            figure(6); imagesc(x, x, maskFillet); axis xy equal tight; colorbar; colormap parula; drawnow;
        end
        
        for iInterior = 1:length(grayInds)

            subpixel = 0*subpixel;

            xCenter = X(grayInds(iInterior));
            yCenter = Y(grayInds(iInterior));
            Xup = Xup0 + xCenter;
            Yup = Yup0 + yCenter;
            RHOSup = sqrt((Xup).^2 + (Yup).^2);

            Xupf = Xup - xFilletCenterVec(iFillet);
            Yupf = Yup - yFilletCenterVec(iFillet);
            [THETASupf, RHOSupf] = cart2pol(Xupf, Yupf);

            if xDeltaVec(iFillet) > 0
                pmf = 1;
            else
                pmf = -1;
            end
            if istrut == 1
                pms = 1;
            elseif istrut == 3
                pms = -1;
            end
            subpixel(  (RHOSupf.^2>=rocFillet.^2) & (pmf*Yup >= pmf*(mSecantVec(iFillet)*Xup + bSecantVec(iFillet))) &...
                (RHOSup>=ID/2. & RHOSup <= OD/2. & ((Xup <= slopeA*Yup & Xup >= -slopeB*Yup & Yup > 0) | (Xup >= slopeA*Yup & Xup <= -slopeB*Yup & Yup < 0)))) = 1;
            
%             subpixel(  (RHOSupf.^2>=rocFillet.^2) & (pmf*Yup >= pmf*(mSecantVec(iFillet)*Xup + bSecantVec(iFillet))) &...
%                 (RHOSup>=ID/2. & RHOSup <= OD/2. & ((Yup <= slopeB*Xup & Yup >= -slopeA*Xup & Xup > 0) | (Yup >= slopeB*Xup & Yup <= -slopeA*Xup & Xup < 0)))) = 1;
%             
            pixelValue = sum(subpixel(:))/upsampleFactor^2;
            pupil2(grayInds(iInterior)) = pixelValue;
            

        end
        if(flagPlot)
            figure(14); imagesc(x, x, pupil2); axis xy equal tight; colorbar; colormap parula; drawnow;
        end
    end
    if(flagPlot)
        hold off
    end

    grayMap = zeros(size(mask));
    grayMap(grayInds) = 1;

    if(flagPlot)
        figure(21); imagesc(x, x, grayMap); axis xy equal tight; colorbar; drawnow;
        figure(22); imagesc(x, x, pupil2); axis xy equal tight; colorbar; drawnow;
        figure(23); imagesc(x, x, pupil-pupil2); axis xy equal tight; colorbar; drawnow;
    end
    pupilCube(:, :, istrut) = pupil2;


end
if(flagPlot)
    hold off
end
%%

mask = pupil - sum(pupilCube, 3);

if(flagPlotFinal)
    figure(20); imagesc(x, x, mask); axis xy equal tight; colorbar; drawnow;
end


end % EOF


%% Custom Functions
function [ang1, ang2] = dualatan2ccw(x1, y1, x2, y2)
    % Function to keep angles increasing when being checked going CCW along
    % a circle.
    
    ang1 = atan2(y1, x1);
    ang2 = atan2(y2, x2);
    
    if (ang1 > 0) && (ang2 < 0)
        ang2 = ang2 + 2*pi;
    end

end


function [ang1, ang2] = dualatan2cw(x1, y1, x2, y2)
    % Function to keep angles increasing when being checked going CW along
    % a circle.
    
    ang1 = atan2(y1, x1);
    ang2 = atan2(y2, x2);
    
    if (ang1 < 0) && (ang2 > 0)
        ang1 = ang1 + 2*pi;
    end

end