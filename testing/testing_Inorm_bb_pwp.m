% data
out.InormHist_BB_PWP = [
    5.34501202328435e-07,...
    5.63370981319209e-07,...
    5.97951804215536e-07,...
    7.03542659933285e-07,...
    8.12642181066931e-07,...
    1.01832753854799e-06,...
    1.25339179995879e-06,...
    1.60380746354332e-06,...
    2.05872487401094e-06,...
    2.63819703712478e-06,...
    3.97634599534390e-06,...
];

out.InormHist_BB_PWP_02bw = [
    5.92245099532892e-07,...
    6.09904039366215e-07,...
    6.20397177896299e-07,...
    6.59808104964928e-07,...
    6.83539141421749e-07,...
    7.47737724130085e-07,...
    7.88709836292809e-07,...
    8.84803182330166e-07,...
    9.49129044677861e-07,...
    1.08838183051661e-06,...
    1.18598369807467e-06,...
];

out.InormHist_PWP = [
    5.56808355283671e-07, ...
    9.02187841745229e-07, ...
    4.22612929991637e-07, ...
    3.89298877570106e-07, ...
    3.92831444244147e-07, ...
    3.85594325158116e-07, ...
    3.93787311673024e-07, ...
    3.78095548108334e-07, ...
    3.85915969770109e-07, ...
    3.78054885364199e-07, ...
    3.86371065005478e-07
];


% number of iterations
iterations = 1:length(out.InormHist_PWP);

% Plot the first dataset (BB PWP)
figure;
hold on;
plot(iterations, out.InormHist_BB_PWP, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', "#0072BD", 'MarkerFaceColor', "#0072BD");
plot(iterations, out.InormHist_BB_PWP_02bw, '--', 'LineWidth', 2, 'MarkerSize',6,  'Marker', 'o', 'Color', "#0072BD", 'MarkerFaceColor', "#0072BD");
plot(iterations, out.InormHist_PWP, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'auto');

hold off;

% Add labels and title
xlabel('Iteration', 'FontSize', 14);
ylabel('log(Inorm)', 'FontSize', 14);
%title('InormHist as a Function of Iteration');
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
yticks([5e-7 1e-6 5e-6]);
legend({'Broadband Estimator 10% BW ','Broadband estimator 2% BW' , 'Classic Estimator 10% BW'}, 'Location', 'northwest');

% Grid and format
grid on;
ax = gca;
ax.YAxis.Exponent = -6;
ax.XAxis.FontSize = 12; % Increase font size of x-axis ticks
ax.YAxis.FontSize = 12; % Increase font size of y-axis ticks

%ax.YAxis.TickLabelFormat = '';

