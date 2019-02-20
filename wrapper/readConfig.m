
function mp = readConfig(coro_type, default_config, user_config)

switch upper(coro_type)
  case 'HLC'
      %% Read HLC user configuration
      uc = ini2struct(user_config);

      %% Add PROPER and FALCO paths
      addpath(genpath(uc.path.falco));
      addpath(genpath(uc.path.proper));

      %% Read default HLC configuration file
      mp = ini2struct(default_config);

      %% Merge default pupil parameters
      mp = mergestruct(mp, mp.(strcat('whichPupil_', uc.whichPupil)));

      %% Merge default controller parameters
      mp = mergestruct(mp, mp.(strcat('controller_', uc.controller)));

      %% Merge default dm9 parameters
      if (isfield(uc.dm9,'inf0name') == false)
        mp = mergestruct(mp, mp.dm9_inf0name);
      else
        mp = mergestruct(mp, mp.(strcat('dm9_inf0name_', uc.dm9.inf0name)));
      end

      %% Merge system and user config
      mp = mergestruct(mp, uc);

      %% Override controller parameters from user config
      mp = mergestruct(mp, uc.(strcat('controller_', uc.controller)));

      %% Resolution must be the same at pupils P1 and P4 for the Vortex and Lyot coronagraphs.
      switch lower(mp.coro)
          case{'lc','hlc','vortex','vc','avc'}
              mp.P4.full.Nbeam = mp.P1.full.Nbeam;  % P4 must be the same as P1 for Vortex and H/LC.
              mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;  % P4 must be the same as P1 for Vortex and H/LC.
      end

      %% Specs for Correction (Corr) region and the Scoring (Score) region.
      if(isfield(mp.F4.corr,'Rin')==false)
          if(isfield(mp.F3,'RinA'))
              mp.F4.corr.Rin  = mp.F3.RinA;
          else
            body
              mp.F4.corr.Rin  = mp.F3.Rin;
          end
      end  %--lambda0/D, inner radius of correction region

      if (isfield(mp.F4.score,'Rin')==false); mp.F4.score.Rin = mp.F4.corr.Rin; end

      if (mp.controller == 'plannedEFC')
          [mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);
      end

  otherwise
      error('Error: Coronagraph %s not supported. Stopping.', coro_type)
end

return
