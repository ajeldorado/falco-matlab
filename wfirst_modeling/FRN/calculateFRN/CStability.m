function CS = CStability(NItoContrast, inputData)
% Compute contrast stability
% 
% S. Miller 11-Jan-2019

% Unpack data
dist = inputData.dist;
distUnits = inputData.distUnits;
sens = inputData.sens;
initRC = inputData.initRC;

% Calculate errors
M.ThermHP  = sum(sens .* (dist((0*21)+(1:21),:) .* distUnits).^2) * 10^-9;
V.ThermHP  = sum(sens .* (dist((1*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dM.ThermHP = sum(sens .* (dist((2*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.ThermHP = sum(sens .* (dist((3*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.ThermLP  = sum(sens .* (dist((4*21)+(1:21),:) .* distUnits).^2) * 10^-9;
V.ThermLP  = sum(sens .* (dist((5*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dM.ThermLP = sum(sens .* (dist((6*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.ThermLP = sum(sens .* (dist((7*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.PupShear  = sum(sens .* (dist((8*21)+(1:21),:)  .* distUnits).^2) * 10^-9;
V.PupShear  = sum(sens .* (dist((9*21)+(1:21),:)  .* distUnits).^2) * 10^-9;
dM.PupShear = sum(sens .* (dist((10*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.PupShear = sum(sens .* (dist((11*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.RWA_TT  = sum(sens .* (dist((12*21)+(1:21),:) .* distUnits).^2) * 10^-9;
V.RWA_TT  = sum(sens .* (dist((13*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dM.RWA_TT = sum(sens .* (dist((14*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.RWA_TT = sum(sens .* (dist((15*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.RWA_WFE  = sum(sens .* (dist((16*21)+(1:21),:) .* distUnits).^2) * 10^-9;
V.RWA_WFE  = sum(sens .* (dist((17*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dM.RWA_WFE = sum(sens .* (dist((18*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.RWA_WFE = sum(sens .* (dist((19*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.DM_setting  = sum(sens .* (dist((20*21)+(1:21),:) .* distUnits).^2) * 10^-9;
V.DM_setting  = sum(sens .* (dist((21*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dM.DM_setting = sum(sens .* (dist((22*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.DM_setting = sum(sens .* (dist((23*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.DM_Therm  = sum(sens .* (dist((24*21)+(1:21),:) .* distUnits).^2) * 10^-9;
V.DM_Therm  = sum(sens .* (dist((25*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dM.DM_Therm = sum(sens .* (dist((26*21)+(1:21),:) .* distUnits).^2) * 10^-9;
dV.DM_Therm = sum(sens .* (dist((27*21)+(1:21),:) .* distUnits).^2) * 10^-9;

M.CGI_Initial_NI = initRC(1);
V.CGI_Initial_NI = initRC(2);

Mcell = struct2cell(M);
Vcell = struct2cell(V);
dMcell = struct2cell(dM);
dVcell = struct2cell(dV);
M.raw_NI_components = sum(reshape([Mcell{:}], size(Mcell)));
V.raw_NI_components = sum(reshape([Vcell{:}], size(Vcell)));
dM.differential_NI = sqrt(2*M.raw_NI_components * sum(reshape([dMcell{:}], size(dMcell))));
dV.differential_NI = sqrt(sum(reshape([dVcell{:}], size(dVcell)).^2));

sum_raw_NI_components = M.raw_NI_components + V.raw_NI_components;
sum_differential_NI = sqrt(dM.differential_NI^2 + dV.differential_NI^2);

CS.NIcontr.M = M;
CS.NIcontr.V = V;
CS.NIcontr.dM = dM;
CS.NIcontr.dV = dV;

CS.rawContrast = sum_raw_NI_components / NItoContrast;
CS.differentialContrast = sum_differential_NI / NItoContrast;

return
