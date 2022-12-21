function mp = falco_ctrl_strategy_example(mp, metrics, itr)
    % mp = falco_control_strategy_example(mp, metrics, itr)
    %
    % modify parameters like:
    %    exposure time
    %    regularization
    %    probe amplitude
    % based on metrics like normalized intensity, iteration #, or
    % whatever you like
    %
    % metrics = 'out' struct from falco_store_intensities, etc.
    % examples from Joon's hcim are here:
    % /home/bseo/HCIT/hcim/efc/config/CtrlStrategy_20190701.py
    
    % metrics example:
    % Note: most of the metrics arrays below are size (Nitr, Nsbp)
    %
    %                   Nitr: 150
    %           log10regHist: [150×1 double]
    %                   ctrl: [1×1 struct]
    %                    dm1: [1×1 struct]
    %                    dm2: [1×1 struct]
    %                    dm8: [1×1 struct]
    %                    dm9: [1×1 struct]
    %                  Zsens: [0×0×150 double]
    %      complexProjection: [149×1 double]
    %     complexCorrelation: [149×1 double]
    %              InormHist: [151×1 double]
    %           IrawCorrHist: [151×1 double]
    %          IrawScoreHist: [151×1 double]
    %           IestCorrHist: [150×1 double]
    %          IestScoreHist: [150×1 double]
    %          IincoCorrHist: [150×1 double]
    %         IincoScoreHist: [150×1 double]
    %        normIntMeasCorr: [150×1 double]
    %       normIntMeasScore: [150×1 double]
    %         normIntModCorr: [150×1 double]
    %        normIntModScore: [150×1 double]
    %       normIntUnmodCorr: [150×1 double]
    %      normIntUnmodScore: [150×1 double]
    %                  thput: [151×1 double]
    %                   Fend: [1×1 struct]
    %          serialDateVec: [150×1 double]
    %              sbp_width: 0
    %      tb_report_initial: [1×1 struct]
    %                    Itr: 116
    %          datetimeArray: [150×1 datetime]
    %           InormHist_tb: [1×1 struct]
    %            EforSpectra: {1×116 cell}
    %              smspectra: {1×116 cell}
    %                     sm: {1×116 cell}
    %                 alpha2: {1×116 cell}

    % exposure time
    if mean(metrics.normIntMeasCorr(itr,:)) < 1e-7
        mp.tb.info.sbp_texp = 5.0*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass (seconds)
        mp.tb.info.sbp_texp_probe = 5.0*ones(mp.Nsbp,1); % Exposure time for each sub-bandpass when probing (seconds)
    end
    
    % dark hole weighting, does falco have this?
    
    % regularization
    if mean(metrics.normIntMeasCorr(itr,:)) < 1e-7 && mod(itr, 10) == 0 && itr > 21
        % beta bump
        mp.ctrl.dmfacVec = 0.5;        %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
        mp.ctrl.log10regVec = -5;   %--log10 of the regularization exponents (often called Beta values)
        
    elseif mean(metrics.normIntMeasCorr(itr,:)) < 1e-7
        mp.ctrl.dmfacVec = 0.5;        %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
        mp.ctrl.log10regVec = -4:-3;   %--log10 of the regularization exponents (often called Beta values)
    
    end
    
    % probe intensity
    % -- automatic in falco
    % lower the upper bound
    if mean(metrics.normIntMeasCorr(itr,:)) < 1e-7
        mp.est.InormProbeMax = 1e-6; % upper bound
    end

    
end % main function
