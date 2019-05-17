function [ThroughputTable, t_refl] = readThroughputTable(file_dir, filename, mode, centerLambda)
% Find reflection throughput based on mode and center lambda

ThroughputTable = readtable(fullfile(file_dir, filename));
lambda = strcat('lam_', num2str(centerLambda*10^9));

switch mode
    case 'IMG'
        t_refl = ThroughputTable{1,lambda};
    case 'IFS'
        t_refl = ThroughputTable{2,lambda};
end

return
