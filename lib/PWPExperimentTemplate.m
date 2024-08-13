%% PWP


%% Load Defaults
EXAMPLE_defaults_DST_LC_design

mp.SeriesNum = 11;
mp.TrialNum = 1;
mp.relinItrVec = 1;

mp.Nitr = 2;
mp.Nsbp = 1;  

mp.dm1.V = out_pre.out.dm1.Vall(:, :, end);
mp.dm2.V = out_pre.out.dm2.Vall(:, :, end);

 %% Flesh out workspace 
[mp, out] = falco_flesh_out_workspace(mp);
