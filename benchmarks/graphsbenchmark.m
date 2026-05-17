% Benchmark data
Tmax = log10([1e2, 2e2,5e2,1e3]);

EDM  = [ 2.4416, 3.5238, 9.6201,  30.0018];
DM   = [6.1134, 9.2916,  24.3305, 74.7792];
RSSA_dep = [ 9.6303, 13.8283, 31.8367, 125.0828];   % RSSA with dependency graph
%RSSA = [0.28304, 5.4154, 21.9386, 121.3152, 332.309];
FRM  = [6.1484, 8.9848, 25.0621, 77.4519];

figure;
hold on;
scale=100;
plot(Tmax, EDM,  '-s', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(Tmax, DM,   '-o', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(Tmax, RSSA_dep, '-d', 'LineWidth', 1.2, 'MarkerSize', 6);
%plot(Tmax, RSSA, '-^', 'LineWidth', 1.2, 'MarkerSize', 6);
plot(Tmax, FRM,  '-v', 'LineWidth', 1.2, 'MarkerSize', 6);

xlabel('log$_{10}T_{\max}$', 'Interpreter', 'latex');
ylabel('Execution time [s]', 'Interpreter', 'latex');

legend({'EDM', 'DM', 'RSSA*', 'FRM'}, ...
       'Location', 'northwest');

grid on;
set(gca, 'FontSize', 12);

hold off;