% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function [tCoef] = falco_thin_film_material_def(lam, aoi, t_Ni, t_PMGI, pol)
%
% Calculates the thin-film complex transmission for the provided
% combinations of metal and dielectric thicknesses and list of wavelengths.
%
% INPUTS:
%   lam: Wavelength [m]
%   aoi:    Angle of incidense [deg]
%   t_Ni:   Nickel layer thickness [m]
%   t_PMGI: PMGI layer thickness [m]
%   pol: = 0 for TE(s) polarization, = 1 for TM(p) polarization
%
% OUTPUTS:
%   cMask(t_PMGI,t_ni): complex field transmission coeffient. Scalar,
%   complex value.
%
% REVISION HISTORY:
% Modified on 2018-05-01 by A.J. Riggs.
% Created on 2017-12-11 by Erkin Sidick.
% -------------------------------------------------------------------------

function [tCoef, rCoef] = falco_thin_film_material_def(lam, aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, d0, pol)

%--New logic: Titanium layer goes beneath Nickel only. Always include them
%together. Subtract off the thickness of the Ti layer from the intended Ni
%layer thickness.

Nmetal = length(t_Ni_vec);
t_Ti_vec = zeros(Nmetal,1);

for ii = 1:Nmetal
    if(t_Ni_vec(ii) > t_Ti_base) %--For thicker layers
        t_Ni_vec(ii) = t_Ni_vec(ii) - t_Ti_base;
    else %--For very thin layers.
        t_Ti_vec(ii) = t_Ni_vec(ii);
        t_Ni_vec(ii) = 0;
    end
end
% if(t_Ni > t_Ti) %--For thicker layers
%     t_Ni = t_Ni - t_Ti;
% else %--For very thin layers.
%     t_Ti = t_Ni;
%     t_Ni = 0;
% end


lam_nm = lam * 1.0e9;    % m --> nm
theta  = aoi*pi/180;     % deg --> rad

lam_u = lam*1.0e6;
npmgi = 1.524 + 5.176e-03./lam_u.^2 + 2.105e-4./lam_u.^4;

vsilica = [400	1.470127387
405	1.469591804
410	1.469076614
415	1.468580736
420	1.468103157
425	1.467642933
430	1.467199179
435	1.466771066
440	1.466357818
445	1.465958706
450	1.465573044
455	1.46520019
460	1.464839538
465	1.464490518
470	1.464152593
475	1.463825257
480	1.463508031
485	1.463200465
490	1.462902132
495	1.462612628
500	1.462331573
505	1.462058604
510	1.46179338
515	1.461535575
520	1.461284882
525	1.461041008
530	1.460803676
535	1.460572623
540	1.460347597
545	1.460128359
550	1.459914683
555	1.459706352
560	1.45950316
565	1.459304911
570	1.459111418
575	1.458922501
580	1.45873799
585	1.458557722
590	1.458381542
595	1.4582093
600	1.458040855
605	1.457876069
610	1.457714813
615	1.457556962
620	1.457402396
625	1.457251001
630	1.457102667
635	1.456957289
640	1.456814766
645	1.456675001
650	1.4565379
655	1.456403375
660	1.45627134
665	1.456141712
670	1.456014412
675	1.455889364
680	1.455766494
685	1.455645731
690	1.455527009
695	1.455410261
700	1.455295425
705	1.455182439
710	1.455071246
715	1.454961789
720	1.454854014
725	1.454747868
730	1.4546433
735	1.454540263
740	1.454438708
745	1.454338591
750	1.454239867
755	1.454142494
760	1.454046432
765	1.453951639
770	1.453858079
775	1.453765714
780	1.453674508
785	1.453584427
790	1.453495437
795	1.453407505
800	1.453320601
805	1.453234694
810	1.453149754
815	1.453065752
820	1.452982662
825	1.452900457
830	1.452819109
835	1.452738595
840	1.45265889
845	1.45257997
850	1.452501812
855	1.452424394
860	1.452347694
865	1.452271691
870	1.452196365
875	1.452121696
880	1.452047665
885	1.451974252
890	1.451901441
895	1.451829213
900	1.451757551
905	1.451686439
910	1.45161586
915	1.451545798
920	1.451476239
925	1.451407167
930	1.451338567
935	1.451270427
940	1.451202731
945	1.451135467
950	1.451068621
955	1.451002181
960	1.450936135
965	1.45087047
970	1.450805174
975	1.450740237
980	1.450675648
985	1.450611394
990	1.450547467
995	1.450483855
1000	1.450420548];

lam_silica = vsilica(:,1);  % nm
n_silica   = vsilica(:,2);
nsilica    = interp1(lam_silica, n_silica, lam_nm, 'linear');

% ---------------------------------------------
vnickel = [387.5	1.61	2.3
400	1.61	2.36
413.3	1.61	2.44
427.5	1.62	2.52
442.8	1.62	2.61
459.2	1.64	2.71
476.9	1.66	2.81
495.9	1.67	2.93
516.6	1.71	3.06
539.1	1.75	3.19
563.6	1.8	3.33
590.4	1.85	3.48
619.9	1.93	3.65
636	1.98	3.74
653	2.02	3.82
670	2.08	3.91
689	2.14	4
709	2.21	4.09
729	2.28	4.18
751	2.36	4.25
775	2.43	4.31
800	2.48	4.38
827	2.53	4.47
855	2.59	4.55
886	2.65	4.63
918	2.69	4.73
954	2.74	4.85
992	2.8	4.97
1033	2.85	5.1]; 

lam_nickel = vnickel(:,1);  % nm
n_nickel   = vnickel(:,2);
k_nickel   = vnickel(:,3);
nnickel    = interp1(lam_nickel, n_nickel, lam_nm, 'linear');
knickel    = interp1(lam_nickel, k_nickel, lam_nm, 'linear');

% ---------------------------------------------
titanium = [413.0000    2.1400    2.9800
  431.0000    2.2100    3.0100
  451.0000    2.2700    3.0400
  471.0000    2.3200    3.1000
  496.0000    2.3600    3.1900
  521.0000    2.4400    3.3000
  549.0000    2.5400    3.4300
  582.0000    2.6000    3.5800
  617.0000    2.6700    3.7200
  659.0000    2.7600    3.8400
  704.0000    2.8600    3.9600
  756.0000    3.0000    4.0100
  821.0000    3.2100    4.0100
  892.0000    3.2900    3.9600
  984.0000    3.3500    3.9700];

lam_ti = vnickel(:,1);  % nm
n_ti   = vnickel(:,2);
k_ti   = vnickel(:,3);
nti    = interp1(lam_ti, n_ti, lam_nm, 'linear');
kti    = interp1(lam_ti, k_ti, lam_nm, 'linear');
% ---------------------------------------------
Nmetal = length(t_Ni_vec);
Ndiel = length(t_PMGI_vec);

tCoef = zeros(Ndiel,Nmetal); %--initialize2
rCoef = zeros(Ndiel,Nmetal); %--initialize

for jj = 1:Ndiel
    dpm = t_PMGI_vec(jj);
    
    for ii = 1:Nmetal
        dni = t_Ni_vec(ii);
        dti = t_Ti_vec(ii);
        
        nvec = [1 1 npmgi nnickel-1i*knickel nti-1i*kti nsilica];
        %dvec = [d0-dpm dpm dni];
        dvec = [d0-dpm-dni-dti dpm dni dti];
        
%         [R, T, rr, tt] = thin_film_filter_2(nvec, dvec, theta, lam, pol);
        %--OUTPUTS:
            % R = normalized reflected intensity coefficient
            % T =        "   transmitted     "
            % rr = complex field reflection coefficient
            % tt =      "        transmission "
        %amp(jj,ii) = tt; %sqrt(T);
        %pha(jj,ii) = atan2(imag(tt), real(tt));
        
        [~, ~, rr, tt] = falco_thin_film_solver(nvec, dvec, theta, lam, pol);
        tCoef(jj,ii) = tt; %--Complex field transmission coeffient
        rCoef(jj,ii) = rr; %--Complex field transmission coeffient
    end
end

%pha = 0;
        
return
% ---------------------------------------------
%
