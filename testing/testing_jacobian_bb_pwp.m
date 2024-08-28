Npairs = 6;
Nact = 36;
Npix =mp.Fend.corr.Npix; %2717 here

DM1Vnom = load('/home/hcst/falco-matlab/data/DM1Vnom_pwp.mat');
DM2Vnom = mp.dm2.V;
whichDM = mp.est.probe.whichDM;
DM1Vnom = DM1Vnom.DM1Vnom;


%DM1Vplus  = rand([Nact, Nact, Npairs]) + DM1Vnom; %rand part is dDM1Vprobe
%DM1Vminus = zeros([Nact, Nact, Npairs]);
%DM2Vplus  = rand([Nact, Nact, Npairs]) + DM2Vnom;
%DM2Vminus = zeros([Nact, Nact, Npairs]);

DM1Vplus  =  load('/home/hcst/falco-matlab/data/DM1Vplus_pwp.mat');
DM2Vplus  =  1;
%jacStruct = load('/home/hcst/falco-matlab/data/jac/jacStruct_pwp.mat');
DM1Vplus = DM1Vplus.DM1Vplus;


dm_inds = [1];

zAll = load('/home/hcst/falco-matlab/data/zAll_pwp.mat');
zAll = zAll.zAll;


Eest_bb = pairwise_bb_estimation(mp, jacStruct, DM1Vplus, DM2Vplus, zAll);
%%
save('/home/hcst/falco-matlab/data/Eest_bb.mat', 'Eest_bb');




