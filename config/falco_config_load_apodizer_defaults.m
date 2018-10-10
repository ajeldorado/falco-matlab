% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% mp = falco_config_load_apodizer_defaults(mp,DM)
%
% Function to store mask configuration data for apodized coronagraphs.
% -Add more apodizers as desired.
%
% 
% REVISION HISTORY:
% ----------------
% Created on 2018-05-29 by A.J. Riggs.




function mp = falco_config_load_apodizer_defaults(mp)


%% Apodizer (Shaped Pupil Properties (Plane P3)

if(mp.flagApod == false); mp.SPname = 'none'; end
% if(isfield(mp,'SPname')==false); mp.SPname = 'luvoirA5bw10'; end

%--Load the inner and outer radii of the apodizer's concentric rings
switch mp.whichPupil
    case{'WFIRST_onaxis','WFIRST20180103'}
        switch mp.SPname
            case '20170714' %-- 18% Bandwidth, NOT concentric rings
                %--Dummy values (since it is not a concentric-ring mask)
                mp.rEdgesLeft = 1;
                mp.rEdgesRight = 1;
            case '32WA194' %-- 18% Bandwidth, concentric rings
                mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot100_3300Dpup9850_32WA194_45LS79_8FPres4_BW18N9_c90_Nring9_D20mm_10um_left_WS.dat');
                mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot100_3300Dpup9850_32WA194_45LS79_8FPres4_BW18N9_c90_Nring9_D20mm_10um_right_WS.dat');  
            case '31WA220' %-- 10% Bandwidth,  concentric rings
                mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_3300Dpup9850_31WA220_43LS82_8FPres4_BW10N9_c90_Nring10_D20mm_10um_left_WS.dat');
                mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_3300Dpup9850_31WA220_43LS82_8FPres4_BW10N9_c90_Nring10_D20mm_10um_right_WS.dat');  
            case '27WA166' %-- 10% Bandwidth,  concentric rings (worse initial contrast)
                mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot400_3300Dpup9850_27WA166_45LS91_8FPres4_BW10N9_c70_Nring8_D20mm_10um_left_WS.dat');
                mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot400_3300Dpup9850_27WA166_45LS91_8FPres4_BW10N9_c70_Nring8_D20mm_10um_right_WS.dat');
        end

    case{'LUVOIRA0'}
        if(strcmpi(mp.SPname,'luvoirA0bw10')) %-- 10% Bandwidth
            mp.rEdgesLeft  = load('out_RSPLC1D_A_maxTrPH_2848Dpup9960_33WA228_38LS83_10FPres5_BW10N11_c100_Nring10_left.dat');
            mp.rEdgesRight = load('out_RSPLC1D_A_maxTrPH_2848Dpup9960_33WA228_38LS83_10FPres5_BW10N11_c100_Nring10_right.dat');  
        end  
        
    case{'LUVOIRA5'}
        if(strcmpi(mp.SPname,'luvoirA5bw10')) %-- 10% Bandwidth
            mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_1123Dpup9900_34WA264_24LS78_10FPres4_BW10N6_c100_Nring14_D20mm_10um_left_WS.dat');
            mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot200_1123Dpup9900_34WA264_24LS78_10FPres4_BW10N6_c100_Nring14_D20mm_10um_right_WS.dat');  
        end
        
    case{'LUVOIR_B_offaxis'}
        if(strcmpi(mp.SPname,'luvoirBbw20')) %-- 10% Bandwidth
            mp.rEdgesLeft  = load('out_RSPLC1D_nl_maxTrPH_Nlyot400_0Dpup9900_25WA350_6LS74_10FPres4_BW20N9_c103_Nring19_D20mm_10um_left_WS.dat');
            mp.rEdgesRight = load('out_RSPLC1D_nl_maxTrPH_Nlyot400_0Dpup9900_25WA350_6LS74_10FPres4_BW20N9_c103_Nring19_D20mm_10um_right_WS.dat');  
        end
end

%--Scale the ring diameters from the inscribed diameter values to the circumscribed diameter values. 
if(mp.flagApod)
    mp.rEdgesLeft = mp.rEdgesLeft/mp.P1.Dfac;
    mp.rEdgesRight = mp.rEdgesRight/mp.P1.Dfac;
end

%% FPM Properties

mp.F3.ang = 180; % default, per-side opening angle of FPM (degrees)

switch mp.coro
    case{'HLC','APHLC','SPHLC'}
        if(strcmpi(mp.SPname,'20170714'))
            mp.F3.Rin = 2.6; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 9.0; % outer hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.ang = 65; % opening angle of horizontal bowtie (degrees)
        elseif(strcmpi(mp.SPname,'32WA194'))
            mp.F3.Rin = 3.2; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 19.4; % outer hard-edge radius of the focal plane mask, in lambda0/D
        elseif(strcmpi(mp.SPname,'31WA220'))
            mp.F3.Rin = 3.1; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 22.0; % outer hard-edge radius of the focal plane mask, in lambda0/D
        elseif(strcmpi(mp.SPname,'27WA166'))
            if(isfield(mp.F4.corr,'Rin')==false)
                if(isfield(mp.F3,'RinA')); mp.F4.corr.Rin  = mp.F3.RinA;
                else; mp.F3.Rin = 2.7; mp.F4.corr.Rin  = mp.F3.Rin;
                end
            end  %--lambda0/D, inner radius of correction region
            if(isfield(mp.F3,'Rout')==false); mp.F3.Rout = 16.6; end % outer hard-edge radius of the focal plane mask, in lambda0/D    

        %--LUVOIR
        elseif(strcmpi(mp.SPname,'luvoirA0bw10'))
            mp.F3.Rin = 3.3*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
            mp.F3.Rout = 22.8*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
        elseif(strcmpi(mp.SPname,'luvoirA5bw10'))
        %     if( exist('mp.P1.Dfac','var')==false ); mp.P1.Dfac = 15.2/13.7; end %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
            mp.F3.Rin = 3.367*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
            mp.F3.Rout = 26.4*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
        elseif(strcmpi(mp.SPname,'luvoirBbw20'))
        %     if( exist('mp.P1.Dfac','var')==false ); mp.P1.Dfac = 15.2/13.7; end %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
            mp.F3.Rin = 2.5*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
            mp.F3.Rout = 35.0*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
        end 
        
        
        %--Convert to one-step propagation in SPHLC if the inner region is too large and would overlap with the outer iris.
        if(strcmpi(mp.coro,'SPHLC'))
            if( mp.F3.Rin*sqrt(2) >= mp.F3.Rout ) %--If the inner regions corners would start overlapping with the outer FPM iris.
                mp.F3.Rin = mp.F3.Rout; %--Have the inner region be large enough to contain the outer iris
                mp.flagSPHLConeFPM = true
            else
                mp.flagSPHLConeFPM = false
            end
        end
        
        
    otherwise
        
        %--WFIRST
        if(strcmpi(mp.SPname,'20170714'))
            mp.F3.Rin = 2.6; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 9.0; % outer hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.ang = 65; % opening angle of horizontal bowtie (degrees)
        elseif(strcmpi(mp.SPname,'32WA194'))
            mp.F3.Rin = 3.2; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 19.4; % outer hard-edge radius of the focal plane mask, in lambda0/D
        elseif(strcmpi(mp.SPname,'31WA220'))
            mp.F3.Rin = 3.1; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 22.0; % outer hard-edge radius of the focal plane mask, in lambda0/D
        elseif(strcmpi(mp.SPname,'27WA166'))
            mp.F3.Rin = 2.7; % inner hard-edge radius of the focal plane mask, in lambda0/D
            mp.F3.Rout = 16.6; % outer hard-edge radius of the focal plane mask, in lambda0/D

        %--LUVOIR
        elseif(strcmpi(mp.SPname,'luvoirA0bw10'))
            mp.F3.Rin = 3.3*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
            mp.F3.Rout = 22.8*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
        elseif(strcmpi(mp.SPname,'luvoirA5bw10'))
        %     if( exist('mp.P1.Dfac','var')==false ); mp.P1.Dfac = 15.2/13.7; end %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
            mp.F3.Rin = 3.367*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
            mp.F3.Rout = 26.4*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
        elseif(strcmpi(mp.SPname,'luvoirBbw20'))
        %     if( exist('mp.P1.Dfac','var')==false ); mp.P1.Dfac = 15.2/13.7; end %--Ratio of OD_circumscribed to OD_inscribed for the non-circular outer aperture.
            mp.F3.Rin = 2.5*mp.P1.Dfac; % inner hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
            mp.F3.Rout = 35.0*mp.P1.Dfac; % outer hard-edge radius of the focal plane mask, in lambda0/Dcircumscribed
        end
        
end



%% Lyot Stop Properties
switch mp.whichPupil
	case{'Simple','SimplePROPER'}
        %--Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle
        mp.P4.full.Nbeam = mp.P1.full.Nbeam; %--Number of pixels across the re-imaged pupil at the Lyot plane (independent of beam centering)
        mp.P4.full.Nbeam
        mp.P4.compact.Nbeam = mp.P1.compact.Nbeam; %--Number of pixels across the aperture or beam (independent of beam centering)
        %--Make or read in Lyot stop (LS) for the 'full' model
        if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.3; end % Inner radius of Lyot stop
        if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.85; end % Outer radius of Lyot stop 
        mp.LS_num_strut = 4; % Number of struts in Lyot stop 
        mp.LS_strut_angs = [0 90 180 270];%Angles of the struts 
        mp.LS_strut_width = 0.01;% Size of Lyot stop spiders 
        
    case{'WFIRST_onaxis','WFIRST20180103'} % WFIRST is case specific 
        %--Make or read in Lyot stop (LS) for the 'full' model
        if(strcmpi(mp.SPname,'20170714'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.25; mp.P4.ODnorm = 0.82; mp.P4.ang = 115;  end
            mp.LS_strut_width = 0; %--Dummy value
        elseif(strcmpi(mp.SPname,'32WA194'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.45; mp.P4.ODnorm = 0.79; end
        elseif(strcmpi(mp.SPname,'31WA220'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.43; mp.P4.ODnorm = 0.82; end
        elseif(strcmpi(mp.SPname,'27WA166'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.45; mp.P4.ODnorm = 0.91; end
        end

    case{'LUVOIRA0'} % LUVOIR with big central obscuration
        if(strcmpi(mp.SPname,'luvoirA0bw10'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.38/mp.P1.Dfac; end
            if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.83/mp.P1.Dfac; end
        end
        
    case{'LUVOIRA5'} % LUVOIR with small central obscuration
        if(strcmpi(mp.SPname,'luvoirA5bw10'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.24/mp.P1.Dfac; end
            if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.78/mp.P1.Dfac; end
        end
        
        
    case{'LUVOIR_B_offaxis'} % Segmented, off-axis LUVOIR
        if(strcmpi(mp.SPname,'luvoirBbw20'))
            if(isfield(mp.P4,'IDnorm')==false);  mp.P4.IDnorm = 0.06/mp.P1.Dfac; end
            if(isfield(mp.P4,'ODnorm')==false);  mp.P4.ODnorm = 0.74/mp.P1.Dfac; end
        end
end



end %--END OF FUNCTION



