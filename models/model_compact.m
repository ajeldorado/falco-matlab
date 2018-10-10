% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_compact(mp,   modvar)
%--Blind model used by the estimator and controller
%  Does not include unknown aberrations/errors that are in the full model.
%
% REVISION HISTORY:
% --------------
% Modified on 2017-10-17 by A.J. Riggs to have model_compact.m be a wrapper. All the 
%  actual compact models have been moved to sub-routines for clarity.
% Modified on 19 June 2017 by A.J. Riggs to use lower resolution than the
%   full model.
% model_compact.m - 18 August 2016: Modified from hcil_model.m
% hcil_model.m - 18 Feb 2015: Modified from HCIL_model_lab_BB_v3.m
% ---------------
%
% INPUTS:
% -mp = structure of model parameters
% -DM = structure of DM settings
% -modvar = structure of model variables
%
%
% OUTPUTS:
% -Eout
%  -> Note: When computing the control Jacobian, Eout is a structure
%  containing the Jacobians. Otherwise, it is just the E-field at the
%  second focal plane.
%
% modvar structure fields (4):
% -sbpIndex
% -wpsbpIndex
% -whichSource
% -flagGenMat


function Eout = model_compact(mp, modvar,varargin)

modvar.wpsbpIndex = 0; %--Dummy index since not needed in compact model

% Set default values of input parameters
normFac = mp.F4.compact.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
flagEval = false;             % flag to use a different (usually higher) resolution at final focal plane for evaluation
flagNewNorm = false;
%--Enable different arguments values by using varargin
icav = 0;                     % index in cell array varargin
while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
        case {'normoff','unnorm','nonorm'} % Set to 0 when finding the normalization factor
            normFac = 0; 
            flagNewNorm = true;
            %fprintf('model_compact: Not normalized.\n');
        case {'eval'} % Set to 0 when finding the normalization factor
            flagEval = true; 
        otherwise
            error('model_compact: Unknown keyword: %s\n', varargin{icav});
          
    end
end

%--Normalization factor for compact evaluation model
if( (flagNewNorm==false) && (flagEval==true) )
    normFac = mp.F4.eval.I00(modvar.sbpIndex); % Value to normalize the PSF. Set to 0 when finding the normalization factor
end

%--Set the wavelength
if(isfield(modvar,'lambda'))
    lambda = modvar.lambda;
else
    lambda = mp.sbp_centers(modvar.sbpIndex);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input E-fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Include the tip/tilt in the input wavefront
if(isfield(mp,'ttx'))
    %--Scale by lambda/lambda0 because ttx and tty are in lambda0/D
    x_offset = mp.ttx(modvar.ttIndex)*(mp.lambda0/lambda);
    y_offset = mp.tty(modvar.ttIndex)*(mp.lambda0/lambda);

    TTphase = (-1)*(2*pi*(x_offset*mp.P2.compact.XsDL + y_offset*mp.P2.compact.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex);  

elseif strcmpi(modvar.whichSource,'offaxis') %--Use for throughput calculations 
    TTphase = (-1)*(2*pi*(modvar.x_offset*mp.P2.compact.XsDL + modvar.y_offset*mp.P2.compact.YsDL));
    Ett = exp(1i*TTphase*mp.lambda0/lambda);
    Ein = Ett.*mp.P1.compact.E(:,:,modvar.sbpIndex); 
    
else %--Backward compatible with code without tip/tilt offsets in the Jacobian
    Ein = mp.P1.compact.E(:,:,modvar.sbpIndex);  
end



%--Apply a Zernike (in amplitude) at input pupil if specified
if(isfield(modvar,'zernIndex')==false)
    modvar.zernIndex = 1;
end
%--Only used for Zernike sensitivity control, which requires the perfect 
% E-field of the differential Zernike term.
if(modvar.zernIndex~=1)
    indZernVec = find(mp.jac.zerns==modvar.zernIndex); %--Index in vector of RMS values for Zernikes.
    indsZnoll = modvar.zernIndex; %--Just send in 1 Zernike mode
    zernMat = falco_gen_norm_zernike_maps(mp.P1.compact.Nbeam,mp.centering,indsZnoll); %--Cube of normalized (RMS = 1) Zernike modes.
    % figure(1); imagesc(zernMat); axis xy equal tight; colorbar; 
    Ein = Ein.*zernMat*(2*pi*1i/lambda)*mp.jac.Zcoef(indZernVec);
end


%--Select the type of coronagraph
switch mp.coro 
    
    case{'EHLC'} %--DMs, optional apodizer, extended FPM with metal and dielectric modulation and outer stop, and LS. Uses 1-part direct MFTs to/from FPM
        %--Complex transmission map of the FPM.
        FPM = falco_gen_EHLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact');
        Eout = model_compact_EHLC(mp,   lambda, normFac, Ein, FPM, flagEval);
                
    case{'HLC','APHLC'} %--DMs, optional apodizer, FPM with optional metal and dielectric modulation, and LS. Uses Babinet's principle about FPM.
        %--Complex transmission map of the FPM.
        FPM = falco_gen_HLC_FPM_complex_trans_mat( mp,modvar.sbpIndex,modvar.wpsbpIndex,'compact');
        Eout = model_compact_HLC(mp,   lambda, normFac, Ein, FPM, flagEval);
        
    case{'SPHLC','FHLC'} %--DMs, optional apodizer, complex/hybrid FPM with outer diaphragm, LS. Uses 2-part direct MFTs to/from FPM
        %Eout = model_compact_SPHLC(mp,   lambda, Ein, normFac);
        
        
    
    case{'LC','DMLC','APLC'} %--DMs, optional apodizer, FPM with/without phase contribution, and LS.
        Eout = model_compact_LC(mp,   lambda, Ein, normFac, flagEval);  
       
    case{'SPLC','FLC'} %--DMs, optional apodizer, binary-amplitude FPM with outer diaphragm, LS
        Eout = model_compact_SPLC(mp,   lambda, Ein, normFac, flagEval);
            
    case{'vortex','Vortex','VC','AVC'} %--DMs, optional apodizer, vortex FPM, LS
        Eout = model_compact_VC(mp, lambda, Ein, normFac, flagEval);      
        
%     case{'SPC','APP','APC'} %--Pupil-plane mask only
%         Eout = model_compact_APC(mp,   modvar);             

    otherwise
        disp('ERROR: CASE NOT RECOGNIZED IN model_compact.m');        
end
    

end % End of function


    
