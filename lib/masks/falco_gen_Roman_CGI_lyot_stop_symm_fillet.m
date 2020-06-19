% Generate, plot, and write etching coordinates and bitmap for the Roman CGI SPC-WFOV Lyot Stop.
% Use pupil CGI-20200513 as the starting point.
%
% Note: The strut width specified by wStrut is the width before
% symmetrizing the pupil. The actual strut width is larger because of the
% symmetrization.
%
% Written by A.J. Riggs on June 2, 2020.


function lyot = falco_gen_Roman_CGI_lyot_stop_symm_fillet(Nbeam, ID, OD, wStrut, rocFillet, upsampleFactor, centering)


flagPlot = false;
flagPlotFinal = false;

% wHexUM = 21e3; % flat-to-flat width of hexagon [microns]
% wHexOutUM = 23e3; % flat-to-flat width of outer hexagon for making a dropout inner region [microns]
% RoCHexUM = 3.175e3; % radius of curvature at vertices [microns]
% 
% aoiDeg = 0;%7.6; % [degrees]
% aoiAxis = 'x'; % 'x' or 'y'
% flipx = false;
% rotDeg = 10; % counterclockwise rotation of the Lyot stop [degrees]

DbeamUM = 17.000e3; % pupil diameter (major axis) arriving at Lyot plane [microns]
RfilletUM = rocFillet*DbeamUM;%170;%500;%340; %85; %170.; %--width of the dropout region [microns]
% 
% ID = 0.36; % LS inner diameter [pupil diameteters]
% OD = 0.91; % LS inouterner diameter [pupil diameteters]
% wStrut = 3.2/100.; % strut width [pupil diameteters]
% 
stepSizeUM = 1.0; % Distance between polygon points. [microns]
% 
% %--Required Inputs
% upsampleFactor = 100;
% centering = 'pixel';%'interpixel';
% Nbeam = 1000;%120;%1000; % max aperture radius in samples
% % % Narray = 1024;%120;%1024; % Number of samples across in square output array
switch centering
    case 'interpixel'
        Narray = ceil_even(Nbeam);
    case 'pixel'
%         Narray = ceil_odd(Nbeam); 
        Narray = ceil_even(Nbeam+1);
    otherwise
        error('Value of centering not recognized. Must be pixel or interpixel.');
end
xShear = 0.;
yShear = 0.;
mag = 1.0;



Rfillet = RfilletUM/(DbeamUM); %--Normalized dropout width [pupil diameteters]
nPointsCircle = round(2*pi*DbeamUM/2./stepSizeUM); % Number of points to make a full circle at the OD
nPointsFillet = round(2*pi*RfilletUM/stepSizeUM);

kernel = ones(3);
kernel = kernel/sum(kernel(:));

%%


% x = (-(Narray-1)/2:(Narray-1)/2)/Nbeam;
switch centering
    case 'interpixel'
        x = (-(Narray-1)/2:(Narray-1)/2)/Nbeam;
    case 'pixel'
        x = (-Narray/2:(Narray/2-1))/Nbeam;
    otherwise
        error('Value of centering not recognized. Must be pixel or interpixel.');
end
y = x;

x = x/mag;
y = y/mag;

x = x - xShear;
y = y - yShear;

[X, Y] = meshgrid(x,y);
dx = x(2) - x(1);
% radius = 0.5;
[~, RHOS] = cart2pol(X, Y);

%% Pupil Parameters

nStrut = 6;

strutWidthVec = wStrut*ones(1, nStrut);

primaryRadiusYpixels = 4027.25;
ODpixels = 2.*primaryRadiusYpixels;

% primaryRadiusX = 3990.0/ODpixels;
% primaryRadiusY = primaryRadiusYpixels/ODpixels;
% primaryCenterX = 0.;
% primaryCenterY = 0.;
% 
% secondaryRadiusX = 1209.65/ODpixels;
% secondaryRadiusY = 1220.0/ODpixels;
% secondaryCenterX = 0.0/ODpixels;
% secondaryCenterY = -2.95/ODpixels;


strutEndVecX1 = ([843.9, 728.0, 47.5, -192.0, -676.0, -816.65])/ODpixels;
strutEndVecY1 = ([550.85, 580.35, -970.65, -1097.15, 605.55, 458.85])/ODpixels;

strutEndVecX2 = ([1579.9, 3988.0, 2430.0-0.3, -2484.0, -3988.0, -1572.65])/ODpixels;
strutEndVecY2 = ([3866.85, -511.65, -3214.65, -3256.15, -504.45, 3866.85])/ODpixels;

% strutEndVecX1 = ([843.9, 728.0, 47.0+0.5, -192.0, -676.0, -816.65])/ODpixels;
% strutEndVecY1 = ([550.85, 580.35, -970.65, -1097.15, 605.55, 458.85])/ODpixels;
% 
% strutEndVecX2 = ([1579.9, 3988.0, 2430.0-0.3, -2484.0, -3988.0, -1572.65])/ODpixels;
% strutEndVecY2 = ([3866.85, -511.65, -3214.65, -3256.15, -504.45, 3866.85])/ODpixels;

strutCenterVecX = (strutEndVecX1 + strutEndVecX2)/2.0;
strutCenterVecY = (strutEndVecY1 + strutEndVecY2)/2.0;

strutSlopeVec = (strutEndVecY2 - strutEndVecY1)./(strutEndVecX2 - strutEndVecX1);
strutSlopeVecB = -fliplr(strutSlopeVec);

%--Compute y-offsets

yOffsetVec = zeros(nStrut, 1);
for istrut = 1:nStrut
%     yOffsetVec(istrut) = strutCenterVecY(istrut) - strutSlopeVec(istrut)*strutCenterVecX(istrut);

    if istrut == 1 || istrut == 2 || istrut == 3
        yOffsetVec(istrut) = strutCenterVecY(istrut) - strutSlopeVec(istrut)*strutCenterVecX(istrut);
    elseif istrut == 4 || istrut == 5 || istrut == 6
        yOffsetVec(istrut) = strutCenterVecY(istrut) - strutSlopeVec(istrut)*strutCenterVecX(istrut);
    end
end

yOffsetVecB = flipud(yOffsetVec); %--y-offsets going in reverse order

% strutWidthVec = [257.0, 259.0, 258.0, 258.0, 258.5+0.5, 257.0]/ODpixels;

strutAngleVec = atan2(strutEndVecY2-strutEndVecY1, strutEndVecX2-strutEndVecX1)*(180/pi); % degrees

strutAngleVecB = zeros(size(strutAngleVec));
counter = 1;
for istrut = 6:-1:1
    strutAngleVecB(counter) = atan2(-yOffsetVec(istrut), yOffsetVec(istrut)/strutSlopeVec(istrut))*180/pi;
    if counter == 3 || counter == 4
       strutAngleVecB(counter) = strutAngleVecB(counter) - 180; 
    end
    counter = counter + 1;
end

% % strutAngleVec = strutAngleVec + rotDeg;
% % strutAngleVecB = strutAngleVecB + rotDeg;

% tabRadiusVecX = [1342.0+1.0, 1342.0+1.0, 1364.0]/ODpixels;
% tabRadiusVecY = [1352.0+1.0, 1352.0+1.0, 1374.0]/ODpixels;
% tabCenterVecX = [0.0, 0.0, 0.0]/ODpixels;
% tabCenterVecY = [53.85, 53.85, 67.6]/ODpixels;
% 
% lStrut = 0.55;
% 
% deltaAngle = 2.5*pi/16;
% angTabStart = [0.616 - deltaAngle/2.0; 2.54 - deltaAngle/2.0; -1.57 - deltaAngle/2.0];
% angTabEnd   = [0.616 + deltaAngle/2.0; 2.54 + deltaAngle/2.0; -1.57 + deltaAngle/2.0];


%% Determine the ordering of the struts in the data above
if flagPlot
    figure(1); 
    hold on;
    for istrut = 1:nStrut

        plot([strutEndVecX1(istrut), strutEndVecX2(istrut)], [strutEndVecY1(istrut), strutEndVecY2(istrut)], '-bo', 'Linewidth', 0.5*istrut);
        axis xy equal tight;
        xlim(0.5*[-1, 1]); ylim(0.5*[-1, 1])
        drawnow;
        hold on;

        xs = -0.5:0.01:0.5;
        plot(xs, strutSlopeVec(istrut)*xs + yOffsetVec(istrut), ':k')
        xlim(0.5*[-1, 1]); ylim(0.5*[-1, 1])
        drawnow;
        hold on;

    %     pause(1/20);

    end

    for istrut = 4:6
        xs = -0.5:0.01:0.5;
        plot(xs, -strutSlopeVec(istrut)*xs + yOffsetVec(istrut), '--r', 'Linewidth', 1.5)
        xlim(0.5*[-1, 1]); ylim(0.5*[-1, 1])
        drawnow;
        hold on;
    end


    hold off;
end
%% Plot all 6 zones to etch out (but don't include the dropout cuts)

% Rin = ID/2.0;
% Rout = OD/2.0;
% 
% th = linspace(0, 360, 10000);
% 
% figure(2); 
% % plot(ID/2*cosd(th), ID/2*sind(th), '--b', OD/2*cosd(th), OD/2*sind(th), '--b', 'Linewidth', 3);
% axis xy equal tight;
% xlim(0.55*[-1, 1]); ylim(0.55*[-1, 1])
% hold on;
% 
% for istrut = nStrut:-1:1 %1:nStrut
%     xs = -0.5:0.01:0.5;
%     
%     %--Strut A, right edge (when going from ID to OD)
%     mA = strutSlopeVec(istrut);
%     bA = yOffsetVec(istrut) - strutWidthVec(istrut)/2.0/cosd(strutAngleVec(istrut));
%     
%     %--Strut A, inner circle intersection
%     x1ia = (-bA*mA - sqrt(Rin*Rin*(mA*mA + 1) -bA*bA))/(mA*mA+1);
%     x1ib = (-bA*mA + sqrt(Rin*Rin*(mA*mA + 1) -bA*bA))/(mA*mA+1);
%     if istrut == 1 || istrut == 2 || istrut == 3
%         x1i = max([x1ia, x1ib]);
%     else
%         x1i = min([x1ia, x1ib]);
%     end
%     y1i = mA*x1i + bA;
%     
%     %--Strut A, outer circle intersection
%     x1oa = (-bA*mA - sqrt(Rout*Rout*(mA*mA + 1) -bA*bA))/(mA*mA+1);
%     x1ob = (-bA*mA + sqrt(Rout*Rout*(mA*mA + 1) -bA*bA))/(mA*mA+1);
%     if istrut == 1 || istrut == 2 || istrut == 3
%         x1o = max([x1oa, x1ob]);
%     else
%         x1o = min([x1oa, x1ob]);
%     end
%     y1o = mA*x1o + bA;
% %     plot([x1i, x1o], [y1i, y1o], 'rx', 'MarkerSize', 8)
% 
%     
%     %--Strut B, left edge (when going from ID to OD)
%     if istrut == 6
%         istrutB = 1;
%     else
%         istrutB = istrut+1;
%     end
%     mB = strutSlopeVec(istrutB);
%     bB = yOffsetVec(istrutB) + strutWidthVec(istrutB)/2.0/cosd(strutAngleVec(istrutB));
%     
%     %--Strut B, inner circle intersection
%     xm = (-bB*mB - sqrt(Rin*Rin*(mB*mB + 1) -bB*bB))/(mB*mB+1);
%     xp = (-bB*mB + sqrt(Rin*Rin*(mB*mB + 1) -bB*bB))/(mB*mB+1);
%     if istrut == 1 || istrut == 2 || istrut == 6
%         x2i = max([xm, xp]);
%     else
%         x2i = min([xm, xp]);
%     end
%     y2i = mB*x2i + bB;
%     
%     %--Strut B, outer circle intersection
%     xm = (-bB*mB - sqrt(Rout*Rout*(mB*mB + 1) -bB*bB))/(mB*mB+1);
%     xp = (-bB*mB + sqrt(Rout*Rout*(mB*mB + 1) -bB*bB))/(mB*mB+1);
%     if istrut == 1 || istrut == 2 || istrut == 6
%         x2o = max([xm, xp]);
%     else
%         x2o = min([xm, xp]);
%     end
%     y2o = mB*x2o + bB;
% %     plot([x2i, x2o], [y2i, y2o], 'rx', 'MarkerSize', 8)
%     
%     
%     
%     
%     %--Inner circle angles (CW)
%     tia = atan2(y1i, x1i);
%     tib = atan2(y2i, x2i);
% 
%     nPointsInner = round(nPointsCircle*abs(tib-tia)/(2*pi));
%     if istrut == 4
%         tia = 2*pi + tia;
%     end
%     ti = linspace(tia, tib, nPointsInner);
% 
%     %--Outer circle angles (CCW)
%     toa = atan2(y2o, x2o);
%     tob = atan2(y1o, x1o);
%     if istrut == 5
%         tob = tob + 2*pi;
%     end
%     nPointsOuter = round(nPointsCircle*abs(tob-toa)/(2*pi));
%     to = linspace(toa, tob, nPointsOuter);
%     
%     xArcInner = Rin*cos(ti);
%     yArcInner = Rin*sin(ti);
%     
%     xArcOuter = Rout*cos(to);
%     yArcOuter = Rout*sin(to);
%     
%     xAll = [xArcOuter, xArcInner, xArcOuter(1)];
%     yAll = [yArcOuter, yArcInner, yArcOuter(1)];
% 
%     %--Store coordinates for all the mask openings
%     xyMask{istrut} = [xAll; yAll];
%     
%     
% %     plot(xs, strutSlopeVec(istrut)*xs + yOffsetVec(istrut), ':k')
% %     plot(xs, strutSlopeVec(istrut)*xs + yOffsetVec(istrut) - strutWidthVec(istrut)/2.0/cosd(strutAngleVec(istrut)), '-b')
% %     if istrut == 6
% %         plot(xs, strutSlopeVec(1)*xs + yOffsetVec(1), ':k')
% %         plot(xs, strutSlopeVec(1)*xs + yOffsetVec(1) + strutWidthVec(1)/2.0/cosd(strutAngleVec(1)), '-b')
% %     else
% %         plot(xs, strutSlopeVec(istrut+1)*xs + yOffsetVec(istrut+1), ':k')
% %         plot(xs, strutSlopeVec(istrut+1)*xs + yOffsetVec(istrut+1) + strutWidthVec(istrut+1)/2.0/cosd(strutAngleVec(istrut+1)), '-b')
% %     end
% %     plot(xAll, yAll, 'Color',rgb('DeepPink'),'Linewidth', 2);
% 
%     
%     plot(xAll, yAll, '-k','Linewidth', 2); drawnow;
% %     fill(xAll, yAll, 'b'); drawnow;
% %     fill(xAll, yAll, rgb('LightSlateGray')); drawnow;
%     
%   
% end
% hold off



%% Generate normalized coordinates of the dropout cuts for all six regions
% Coordinates are normalized to the unmasked pupil diameter along its major axis.


Rin = ID/2.0;
Rout = OD/2.0;

RinBias = Rin + Rfillet;
RoutBias = Rout - Rfillet;

th = linspace(0, 360, 10000);

if flagPlot
    figure(3); 
    plot(ID/2*cosd(th), ID/2*sind(th), '--b', OD/2*cosd(th), OD/2*sind(th), '--b', 'Linewidth', 0.5);
    axis xy equal tight;
    xlim(0.55*[-1, 1]); ylim(0.55*[-1, 1])
    hold on;
end

% thROC = linspace(0, 360, 361);
% xROC = Rfillet*cosd(thROC);
% yROC = Rfillet*sind(thROC);

% maskNoFillet = zeros(Narray, Narray);
% justFillets = zeros(Narray, Narray);
pupilCube = zeros(Narray, Narray, nStrut);

for istrut = nStrut:-1:1 % 1:nStrut % 
    
%     xs = -0.5:0.01:0.5; %--for plotting only
    
    %--Strut A, right edge (when going from ID to OD)
    mA = strutSlopeVec(istrut);
    bA = yOffsetVec(istrut) - strutWidthVec(istrut)/2.0/cosd(strutAngleVec(istrut));
    %--Strut A, inner circle intersection
    x1ia = (-bA*mA - sqrt(Rin*Rin*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    x1ib = (-bA*mA + sqrt(Rin*Rin*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    if istrut == 1 || istrut == 2 || istrut == 3
        x1i = max([x1ia, x1ib]);
    else
        x1i = min([x1ia, x1ib]);
    end
    y1i = mA*x1i + bA;
    %--Strut A, outer circle intersection
    x1oa = (-bA*mA - sqrt(Rout*Rout*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    x1ob = (-bA*mA + sqrt(Rout*Rout*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    if istrut == 1 || istrut == 2 || istrut == 3
        x1o = max([x1oa, x1ob]);
    else
        x1o = min([x1oa, x1ob]);
    end
    y1o = mA*x1o + bA;
    if flagPlot
        plot([x1i, x1o], [y1i, y1o], 'rx', 'MarkerSize', 8)
    end
    
    %--Strut B, left edge (when going from ID to OD)
    if istrut == 6
        istrutB = 1;
    else
        istrutB = istrut+1;
    end
    mB = strutSlopeVecB(istrutB);
    bB = yOffsetVecB(istrutB) + strutWidthVec(istrutB)/2.0/cosd(strutAngleVecB(istrutB));
    %--Strut B, inner circle intersection
    xm = (-bB*mB - sqrt(Rin*Rin*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    xp = (-bB*mB + sqrt(Rin*Rin*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    if istrut == 1 || istrut == 2 || istrut == 6
        x2i = max([xm, xp]);
    else
        x2i = min([xm, xp]);
    end
    y2i = mB*x2i + bB;
    %--Strut B, outer circle intersection
    xm = (-bB*mB - sqrt(Rout*Rout*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    xp = (-bB*mB + sqrt(Rout*Rout*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    if istrut == 1 || istrut == 2 || istrut == 6
        x2o = max([xm, xp]);
    else
        x2o = min([xm, xp]);
    end
    y2o = mB*x2o + bB;
    if flagPlot
        plot([x2i, x2o], [y2i, y2o], 'rx', 'MarkerSize', 8)
    end
    
    %--Dropout biased from Strut B, left edge (when going from ID to OD)
    if istrut == 6
        istrutB = 1;
    else
        istrutB = istrut+1;
    end
    mB = strutSlopeVecB(istrutB);
    bB = yOffsetVecB(istrutB) + (strutWidthVec(istrutB)/2.0 + Rfillet)/cosd(strutAngleVecB(istrutB));
    %--Dropout biased from Strut B, inner circle intersection
    xm = (-bB*mB - sqrt(RinBias*RinBias*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    xp = (-bB*mB + sqrt(RinBias*RinBias*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    if istrut == 1 || istrut == 2 || istrut == 6
        xfBi = max([xm, xp]);
    else
        xfBi = min([xm, xp]);
    end
    yfBi = mB*xfBi + bB;
    %--Dropout biased from Strut B, outer circle intersection
    xm = (-bB*mB - sqrt(RoutBias*RoutBias*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    xp = (-bB*mB + sqrt(RoutBias*RoutBias*(mB*mB + 1) -bB*bB))/(mB*mB+1);
    if istrut == 1 || istrut == 2 || istrut == 6
        xfBo = max([xm, xp]);
    else
        xfBo = min([xm, xp]);
    end
    yfBo = mB*xfBo + bB;
    if flagPlot
        plot([xfBi, xfBo], [yfBi, yfBo], 'gx', 'MarkerSize', 8);
    end
    
    %--Fillet center biased from Strut A, right edge (when going from ID to OD)
    bA = yOffsetVec(istrut) - (strutWidthVec(istrut)/2.0 + Rfillet)/cosd(strutAngleVec(istrut));
    %--Fillet center biased from Strut A, inner circle intersection
    x1ia = (-bA*mA - sqrt(RinBias*RinBias*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    x1ib = (-bA*mA + sqrt(RinBias*RinBias*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    if istrut == 1 || istrut == 2 || istrut == 3
        xfAi = max([x1ia, x1ib]);
    else
        xfAi = min([x1ia, x1ib]);
    end
    yfAi = mA*xfAi + bA; 
    %--Fillet center biased from Strut A, outer circle intersection
    x1oa = (-bA*mA - sqrt(RoutBias*RoutBias*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    x1ob = (-bA*mA + sqrt(RoutBias*RoutBias*(mA*mA + 1) -bA*bA))/(mA*mA+1);
    if istrut == 1 || istrut == 2 || istrut == 3
        xfAo = max([x1oa, x1ob]);
    else
        xfAo = min([x1oa, x1ob]);
    end
    yfAo = mA*xfAo + bA;
    if flagPlot
        plot([xfAi, xfAo], [yfAi, yfAo], 'gx', 'MarkerSize', 8)
    end    

    %  Intersection Points go clockwise with a before b.
    
    %--Strut A, Outer Fillet 1st Intersection Point (along circle)
    xfAoa = Rout*cos(atan2(yfAo, xfAo));
    yfAoa = Rout*sin(atan2(yfAo, xfAo));
    %--Strut A, Outer Fillet 1st Intersection Point (along strut)
    bA = yOffsetVec(istrut) - (strutWidthVec(istrut)/2.0)/cosd(strutAngleVec(istrut)); % along strut
    bf1o = yfAo + xfAo/mA; % line perpendicular to strut and intersecting fillet center
    xfAob = (bf1o-bA) / (mA + 1/mA);
    yfAob = mA*xfAob + bA;
    [ang1, ang2] = dualatan2ccw(xfAoa-xfAo, yfAoa-yfAo, xfAob-xfAo, yfAob-yfAo); % [radians]
    thfAo = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletAo = Rfillet*cos(thfAo) + xfAo;
    yFilletAo = Rfillet*sin(thfAo) + yfAo;
    
    if flagPlot
        plot([xfAoa, xfAob], [yfAoa, yfAob], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Gold'), 'MarkerSize', 8);
        plot(xFilletAo, yFilletAo, '--k');
    end
    
    %--Strut A, Inner Fillet 1st Intersection Point (along strut)
    bA = yOffsetVec(istrut) - (strutWidthVec(istrut)/2.0)/cosd(strutAngleVec(istrut)); % along strut
    bf1i = yfAi + xfAi/mA; % line perpendicular to strut and intersecting fillet center
    xfAia = (bf1i-bA) / (mA + 1/mA);
    yfAia = mA*xfAia + bA;
    %--Strut A, Inner Fillet 2nd Intersection Point (along circle)
    xfAib = Rin*cos(atan2(yfAi, xfAi));
    yfAib = Rin*sin(atan2(yfAi, xfAi));
    [ang1, ang2] = dualatan2ccw(xfAia-xfAi, yfAia-yfAi, xfAib-xfAi, yfAib-yfAi); % [radians]
    thfAi = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletAi = Rfillet*cos(thfAi) + xfAi;
    yFilletAi = Rfillet*sin(thfAi) + yfAi;
    
    if flagPlot
        plot([xfAia, xfAib], [yfAia, yfAib], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Magenta'), 'MarkerSize', 8);
        plot(xFilletAi, yFilletAi, '--k');
    end
    
    %--Strut B, Inner Fillet 1st Intersection Point (along circle)
    xfBia = Rin*cos(atan2(yfBi, xfBi));
    yfBia = Rin*sin(atan2(yfBi, xfBi));
    %--Strut B, Inner Fillet 2nd Intersection Point (along strut)
    bB = yOffsetVecB(istrutB) + (strutWidthVec(istrutB)/2.0)/cosd(strutAngleVecB(istrutB)); % along strut
    bfB = yfBi + xfBi/mB; % line perpendicular to strut and intersecting fillet center
    xfBib = (bfB-bB) / (mB + 1/mB);
    yfBib = mB*xfBib + bB;

    [ang1, ang2] = dualatan2ccw(xfBia-xfBi, yfBia-yfBi, xfBib-xfBi, yfBib-yfBi); % [radians]
    thfBi = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletBi = Rfillet*cos(thfBi) + xfBi;
    yFilletBi = Rfillet*sin(thfBi) + yfBi;
    
    if flagPlot
        plot([xfBia, xfBib], [yfBia, yfBib], 'o','Color', rgb('Magenta'),'MarkerFaceColor', rgb('Magenta'), 'MarkerSize', 8);
        plot(xFilletBi, yFilletBi, '-k'); 
    end
   
    %--Strut B, Outer Fillet 1st Intersection Point (along strut)
    bB = yOffsetVecB(istrutB) + (strutWidthVec(istrutB)/2.0)/cosd(strutAngleVecB(istrutB)); % along strut
    bfB = yfBo + xfBo/mB; % line perpendicular to strut and intersecting fillet center
    xfBoa = (bfB-bB) / (mB + 1/mB);
    yfBoa = mB*xfBoa + bB;
    %--Strut B, Outer Fillet 2nd Intersection Point (along circle)
    xfBob = Rout*cos(atan2(yfBo, xfBo));
    yfBob = Rout*sin(atan2(yfBo, xfBo));
    
    [ang1, ang2] = dualatan2ccw(xfBoa-xfBo, yfBoa-yfBo, xfBob-xfBo, yfBob-yfBo); % [radians]
    thfBo = linspace(ang1, ang2, round(nPointsFillet*abs(ang2-ang1)/(2*pi))); % [radians]
    xFilletBo = Rfillet*cos(thfBo) + xfBo;
    yFilletBo = Rfillet*sin(thfBo) + yfBo;
    
    if flagPlot
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
    
    xyMask{istrut} = [xAll; yAll];
    
    if flagPlot
        plot(xAll, yAll, '-b','Linewidth', 1); drawnow;
    %     fill(xAll, yAll, 'b'); drawnow;
    %     fill(xAll, yAll, rgb('LightSlateGray')); drawnow;
    end

    %% Define the fillet regions as binary masks

    %--Line cutting off everything not in the fillet region
    xFilletCenterVec = [xfAo, xfAi, xfBi, xfBo];
    yFilletCenterVec = [yfAo, yfAi, yfBi, yfBo];
    mSecantVec = [(yfAoa-yfAob)/(xfAoa-xfAob), (yfAia-yfAib)/(xfAia-xfAib), (yfBib-yfBia)/(xfBib-xfBia), (yfBob-yfBoa)/(xfBob-xfBoa)];
    bSecantVec = [yfAoa - xfAoa*mSecantVec(1), yfAia - xfAia*mSecantVec(2), yfBia - xfBia*mSecantVec(3), yfBoa - xfBoa*mSecantVec(4)]; 
    xDeltaVec = [(xfAoa-xfAob), (xfAia-xfAib), (xfBia-xfBib), (xfBoa-xfBob)];
    
    % Line connecting the inner corners. Need to use this as a boundary to
    % avoid pixels being included on the other side of the pupil.
    mID = (y2i-y1i)/(x2i-x1i);
    bID = y1i - mID*x1i;
    
    nFillet = length(xFilletCenterVec);
%     maskFillet = zeros(size(X));
    subpixel = zeros(upsampleFactor, upsampleFactor);
    allFillets = zeros(Narray, Narray);
    
    for iFillet = 1:nFillet
        Xf = X - xFilletCenterVec(iFillet);
        Yf = Y - yFilletCenterVec(iFillet);
        [~, RHOSf] = cart2pol(Xf, Yf);

        if istrut <= 3
            pmA = -1;
        else
            pmA = 1;
        end
        
        if any(istrut == [1, 2, 6])
            pmB = 1;
        else
            pmB = -1;
        end
        
        if any(istrut == [1, 5, 6])
            pmid = 1;
        else
            pmid = -1;
        end
        
        %--Fillet-specific        
        if istrut == 1
            pmf = -1;
            if iFillet == 1
                pmf = 1;
            end
        elseif istrut == 2
            pmf = 1;
            if iFillet == 4
                pmf = -1;
            end
        elseif istrut == 3
            pmf = -1;
            if any(iFillet == [2, 3])
                pmf = 1;
            end
        elseif istrut == 4
            pmf = 1;
            if any(iFillet == 1)
                pmf = -1;
            end
        elseif istrut == 5
            pmf = -1;
            if any(iFillet == 4)
                pmf = 1;
            end
        elseif istrut == 6
            pmf = 1;
            if any(iFillet == [2, 3])
                pmf = -1;
            end
        end
        
%         maskFillet = maskFillet | ((RHOSf.^2 >= Rfillet.^2) & (pmf*Y >= pmf*(mSecantVec(iFillet)*X + bSecantVec(iFillet)))...
%             & (RHOS>=ID/2. & RHOS <= OD/2. & (pmA*Y >= pmA*(mA*X + bA)) & (pmB*Y >= pmB*(mB*X + bB)) & (pmid*Y >= pmid*(mID*X + bID))  )) ;

        maskFillet = ((RHOSf.^2 >= Rfillet.^2) & (pmf*Y >= pmf*(mSecantVec(iFillet)*X + bSecantVec(iFillet)))...
            & (RHOS>=ID/2. & RHOS <= OD/2. & (pmA*Y >= pmA*(mA*X + bA)) & (pmB*Y >= pmB*(mB*X + bB)) & (pmid*Y >= pmid*(mID*X + bID))  )) ;

        maskFillet = double(maskFillet);
        maskFillet2 = conv2(maskFillet, kernel);
        maskFillet2 = pad_crop(maskFillet2, size(maskFillet));
        grayIndsFillet = find(maskFillet2 > 0);
        
        if(flagPlot)
            figure(9); imagesc(x, x, maskFillet); axis xy equal tight; colorbar; colormap parula; drawnow;
        end
        
        dxUp = dx/upsampleFactor;
        xUp = (-(upsampleFactor-1)/2:(upsampleFactor-1)/2)*dxUp;
        [Xup0, Yup0] = meshgrid(xUp);
        
        %--Make the grayscale version
        for iInterior = 1:length(grayIndsFillet)

            subpixel = 0*subpixel;

            xCenter = X(grayIndsFillet(iInterior));
            yCenter = Y(grayIndsFillet(iInterior));
            Xup = Xup0 + xCenter;
            Yup = Yup0 + yCenter;
            RHOSup = sqrt((Xup).^2 + (Yup).^2);

            Xupf = Xup - xFilletCenterVec(iFillet);
            Yupf = Yup - yFilletCenterVec(iFillet);
            [~, RHOSupf] = cart2pol(Xupf, Yupf);

            subpixel((RHOSupf.^2 >= Rfillet.^2) & (pmf*Yup >= pmf*(mSecantVec(iFillet)*Xup + bSecantVec(iFillet)))...
                & (RHOSup>=ID/2. & RHOSup <= OD/2. & (pmA*Yup >= pmA*(mA*Xup + bA)) & (pmB*Yup >= pmB*(mB*Xup + bB)) & (pmid*Yup >= pmid*(mID*Xup + bID))  )) = 1;
        
            pixelValue = sum(subpixel(:))/upsampleFactor^2;
            allFillets(grayIndsFillet(iInterior)) = pixelValue;
            

        end
        if(flagPlot)
            figure(14); imagesc(x, x, allFillets); axis xy equal tight; colorbar; colormap parula; drawnow;
        end        
    end

%         %--Just the annular segment
%         maskFillet = maskFillet | (RHOS>=ID/2. & RHOS <= OD/2. & (pmA*Y >= pmA*(mA*X + bA)) & (pmB*Y >= pmB*(mB*X + bB)) & (pmid*Y >= pmid*(mID*X + bID)) );
    maskAnnZone = double(RHOS>=ID/2. & RHOS <= OD/2. & (pmA*Y >= pmA*(mA*X + bA)) & (pmB*Y >= pmB*(mB*X + bB)) & (pmid*Y >= pmid*(mID*X + bID)));
    maskAnnZoneBlur = conv2(maskAnnZone, kernel);
    maskAnnZoneBlur = pad_crop(maskAnnZoneBlur, size(maskAnnZone));
    grayIndsAnnZone = find(maskAnnZoneBlur > 0 & maskAnnZoneBlur < 1);

    if flagPlot
        grayMap = zeros(Narray);
        grayMap(grayIndsAnnZone) = 1;
        figure(15); imagesc(x, x, grayMap); axis xy equal tight; colorbar; colormap parula; drawnow;
    end
    
    %--Make the grayscale version
    annZone = maskAnnZone; %zeros(Narray, Narray);
    subpixel = zeros(upsampleFactor, upsampleFactor);
    for iInterior = 1:length(grayIndsAnnZone)

        subpixel = 0*subpixel;

        xCenter = X(grayIndsAnnZone(iInterior));
        yCenter = Y(grayIndsAnnZone(iInterior));
        Xup = Xup0 + xCenter;
        Yup = Yup0 + yCenter;
        RHOSup = sqrt((Xup).^2 + (Yup).^2);

        subpixel(RHOSup>=ID/2. & RHOSup <= OD/2. & (pmA*Yup >= pmA*(mA*Xup + bA)) & (pmB*Yup >= pmB*(mB*Xup + bB)) & (pmid*Yup >= pmid*(mID*Xup + bID))) = 1;

        pixelValue = sum(subpixel(:))/upsampleFactor^2;
        annZone(grayIndsAnnZone(iInterior)) = pixelValue;

    end
    
    pupilCube(:, :, istrut) = annZone - allFillets;
    
    if(flagPlot)
        figure(17); imagesc(x, x, pupilCube(:, :, istrut)); axis xy equal tight; colorbar; drawnow;
    end

end
if flagPlot
    hold off
end

%%
lyot = sum(pupilCube, 3);

if flagPlotFinal
    figure(123); imagesc(x, x, lyot); axis xy equal tight; colorbar; colormap gray; drawnow;
    set(gcf, 'Color', 'w');
    hold on;

    for istrut = 1:nStrut
        plot(xyMask{istrut}(1,:), xyMask{istrut}(2,:), '-.r', 'Linewidth', 2);
        hold on;
    end
    hold off
end

%% Verify symmetry of mask

% if flagPlot
%     xyMaskFliplr = xyMask;
%     for istrut = 1:nStrut
%         xyMaskFliplr{istrut}(1, :) = -xyMaskFliplr{istrut}(1, :);
%     end
% 
%     figure(21);
%     set(gcf, 'Color', 'w');
% 
%     for istrut = 1:nStrut
%         plot(1e-3*xyMask{istrut}(1,:), 1e-3*xyMask{istrut}(2,:), '-k', 'Linewidth', 1);
%         hold on;
%     end
%     for istrut = 1:nStrut
%         plot(1e-3*xyMaskFliplr{istrut}(1,:), 1e-3*xyMaskFliplr{istrut}(2,:), '--r', 'Linewidth', 1);
%         axis xy equal tight;
%         title('Visual Symmetry Inspection', 'Fontsize', 16);
%         xlabel('mm', 'Fontsize', 20);
%         ylabel('mm', 'Fontsize', 20);
%         set(gca, 'Fontsize', 20);
%         hold on;
% 
%     %     xlim(1e-3*0.5*wHexOut*[-1, 1]); ylim(1e-3*0.5*wHexOut*[-1, 1])
%         drawnow;
%     end
%     hold off;
% end

%% Scale, rotate, stretch, flip coordinates as necessary

% for istrut = 1:nStrut
%     
%     %--Scale
%     xyMask{istrut} = xyMask{istrut}*DbeamUM; % for plotting only
%     
%     %--Rotate
%     if rotDeg ~= 0
%         rotMat = [cosd(rotDeg), -sind(rotDeg);...
%                   sind(rotDeg), cosd(rotDeg)];
%         for ii = 1:size(xyMask{istrut}, 2)
%             xyMask{istrut}(:, ii) = rotMat*xyMask{istrut}(:, ii);
%         end        
%     end
%     
%     %--Stretch by 1/cosd(aoiDeg)
%     if lower(aoiAxis) == 'x'
%         xyMask{istrut}(1, :) = xyMask{istrut}(1, :)/cosd(aoiDeg); % for plotting only
%     elseif lower(aoiAxis) == 'y'
%         xyMask{istrut}(2, :) = xyMask{istrut}(2, :)/cosd(aoiDeg); % for plotting only      
%     else
%         error('Axis of AOI must be x or y')
%     end
% 
%     %--Flip in x (i.e., across the y-axis)
%     if flipx
%         xyMask{istrut}(1, :) = -xyMask{istrut}(1, :); % for plotting only       
%     end
%     
% end

%% 
% [xyHexCubeEtch, xyHexCubePlot] = gen_ls_outer_hexagon_for_etching(wHexUM, wHexOutUM, RoCHexUM);



%% Plot how the mask will look
% if flagPlot
%     figure(6);
%     plot(1e-3*xyHexCubePlot(1,:), 1e-3*xyHexCubePlot(2,:), '-k', 'Linewidth', 2); 
%     % fill(1e-3*xyHexCubePlot(1,:), 1e-3*xyHexCubePlot(2,:), rgb('LightGray')); 
%     axis xy equal tight;
%     xlabel('mm', 'Fontsize', 20);
%     ylabel('mm', 'Fontsize', 20);
%     set(gca, 'Fontsize', 20);
%     set(gcf, 'Color', 'w');
%     hold on;
%     for istrut = 1:nStrut
%         plot(1e-3*xyMask{istrut}(1,:), 1e-3*xyMask{istrut}(2,:), '-k', 'Linewidth', 1); axis xy equal tight;
%     %     fill(1e-3*xyMask{istrut}(1,:), 1e-3*xyMask{istrut}(2,:), 'w')
%         xlim(1e-3*0.6*wHexOutUM*[-1, 1]); ylim(1e-3*0.6*wHexOutUM*[-1, 1])
%         drawnow;
%     end
%     hold off;
% 
%     % export_fig('fig_flight_hlc_ls_plotting.png', '-dpng', '-r400');
% end

%% Plot the dropouts that will be etched out.
% if flagPlot
%     figure(4); 
%     xlabel('microns', 'Fontsize', 20);
%     ylabel('microns', 'Fontsize', 20);
%     set(gca, 'Fontsize', 20);
%     set(gcf, 'Color', 'w');
%     hold on;
%     for istrut = 1:nStrut
% 
%         plot(xyMask{istrut}(1,:), xyMask{istrut}(2,:), '-k', 'Linewidth', 2); axis xy equal tight;
%     %     fill(xyMask{istrut}(1,:), xyMask{istrut}(2,:), rgb('Silver'))
%         drawnow;
% 
%         plot(xyHexCubeEtch(1,:,istrut), xyHexCubeEtch(2,:,istrut), '-k', 'Linewidth', 2);
%     %     fill(xyHexCubeEtch(1,:,istrut), xyHexCubeEtch(2,:,istrut), rgb('Silver'))
%         axis xy equal tight;
%         xlim(0.6*wHexOutUM*[-1, 1]); ylim(0.6*wHexOutUM*[-1, 1])
%         drawnow;
%     end
%     hold off;
% 
%     % export_fig('fig_flight_hlc_ls_etching.png', '-dpng', '-r400');
% end

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

