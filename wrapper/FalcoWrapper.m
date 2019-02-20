function FalcoWrapper(coro, default_config, user_config, n_threads, wdir, prefix)

% Wrapper for running FALCO Matlab code in parallel.
%
% Outputs:
%
% Required Inputs:
%   coro = Type of Coronagraph
%   user_config = User configuration file
%   n_threads = Number of parallel threads
%   wdir = working directory
%   prefix = out file prefix
%
% 2018 Nov 1 Navtej Created FALCO wrapper
%

  % Start time for the routine
  start_time = datetime;

  %% Read system and user configuration
  mp = readConfig(coro, default_config, user_config);

  %% Set working directory
  mp.path.config = fullfile(wdir, 'brief/')
  mp.path.ws = fullfile(wdir, 'ws/')

  %% Set out prefix
  mp.runLabel = prefix

  %% Speed up parallel pool launch
  distcomp.feature('LocalUseMpiexec', false);

  %% Set maximum number of threads
  LASTN = maxNumCompThreads(n_threads);

  %% Check if parfor is getting used
  if mp.flagParfor == true
      %% Start Parallel Pool
      t = tempname();
      mkdir(t);
      try
         pc = parcluster();
         pc.JobStorageLocation=t; 
         pause(1+60*rand());    	
         parpool(pc, n_threads);
      catch ME
         disp(ME);
      end
  end
 
  %% Run FALCO
  out = falco_wfsc_loop(mp);

  %% Remove the parpool temporary directory 
  if exist(t, 'dir')
     rmdir(t, 's'); 
  end

  disp(['Elapsed Time : ', datestr(datetime - start_time, 'HH:MM:SS')]);

end
