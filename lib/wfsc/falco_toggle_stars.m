% Copyright 2018-2021, by the California Institute of Technology. ALL RIGHTS
  % RESERVED. United States Government Sponsorship acknowledged. Any
  % commercial use must be negotiated with the Office of Technology Transfer
  % at the California Institute of Technology.
  % -------------------------------------------------------------------------
  %
  % Toggle stars on/off for multi-star wavefront control (MSWC).
  %
  % INPUTS
  % ------
  % mp : structure of model parameters
  % starIndex : which star(s) to turn on:
  %             1 = on-axis star only
  %             2 = off-axis star only
  %             0 or 'both' = both stars
  % initStarWeights : saved initial star weights vector (from mp.star.weights)
  % initTbCurrents : saved initial testbed currents struct with fields:
  %                  .onax  = on-axis star current
  %                  .offax = off-axis star current
  %
  % OUTPUTS
  % -------
  % mp : modified structure of model parameters
  %
  % NOTES
  % -----
  % This function works for model (mp.flagSim == true)
  % This function has been tested only for the OMC testbed (mp.flagSim == false)
  % This function modifies:
  %  - mp.star.weights
  %  - mp.tb.star.current
  %  - mp.tb.offaxisstar.current
  %  - mp.tb.info.PSFpeaks
  % and pauses for mp.est.toggledMSWC_waitTime if in testbed mode.

  function mp = falco_toggle_stars(mp, starIndex, initStarWeights, initTbCurrents)

      % Convert 'both' to 0 for convenience
      if ischar(starIndex) && strcmpi(starIndex, 'both')
          starIndex = 0;
      end

      % Set star weights
      if starIndex == 0  % Both stars
          mp.star.weights = initStarWeights;

          if ~mp.flagSim
              % Restore both testbed sources to initial currents
              mp.tb.star.current = initTbCurrents.onax;
              mp.tb.offaxisstar.current = initTbCurrents.offax;
          end

      elseif starIndex == 1  % On-axis star only
          mp.star.weights = zeros(1, mp.star.count);
          mp.star.weights(1) = initStarWeights(1);

          if ~mp.flagSim
              mp.tb.star.current = mp.tb.info.star_power;
              mp.tb.offaxisstar.current = 0;
              pause(1.0);
          end

      elseif starIndex == 2  % Off-axis star only
          mp.star.weights = zeros(1, mp.star.count);
          mp.star.weights(2) = initStarWeights(2);

          if ~mp.flagSim
              mp.tb.star.current = 0;
              mp.tb.offaxisstar.current = mp.tb.info_offaxisstar.star_power;
              pause(1.0);
          end

      else
          error('starIndex must be 0 (both), 1 (on-axis), or 2 (off-axis)');
      end

      % Ensure non-negative currents and pause for testbed stabilization
      if ~mp.flagSim
          if mp.tb.star.current < 0
              mp.tb.star.current = 0;
          end
          if mp.tb.offaxisstar.current < 0
              mp.tb.offaxisstar.current = 0;
          end

          pause(mp.est.toggledMSWC_waitTime);
      end

  end
