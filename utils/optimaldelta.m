clear;
close all;

% Definizione modello
% BRModelDefinitionNoWater; % se serve
BRModelReducedwithoutIMA;

numSimulations = 20; % Numero di simulazioni per ciascun delta
deltaValues = 0.05:0.01:0.45; % Variazione del delta
averageTimes = zeros(length(deltaValues),1); % Array per tempi medi

for j = 1:length(deltaValues)
    delta = deltaValues(j);
    performanceTimes = zeros(numSimulations, 1); % Preallocazione

    parfor i = 1:numSimulations
        tic;
        [T, Dynamics] = simRSSAOptimized(vMinusBRr, vPlusBRr, cBRr, initialStateBRr, delta, 5e2, 1e-2);
        performanceTimes(i) = toc;
    end

    averageTimes(j) = mean(performanceTimes); % Tempo medio per questo delta
    disp(['Delta = ', num2str(delta), ', Avg Time = ', num2str(averageTimes(j))]);
end

% Trova il delta che minimizza il tempo medio
[~, idxMin] = min(averageTimes);
bestDelta = deltaValues(idxMin);
disp(['Delta ottimale (min tempo medio) = ', num2str(bestDelta)]);

% Salva dati per successivi plot
save('deltaPerformanceData.mat', 'deltaValues', 'averageTimes', 'bestDelta');

% Plot
figure;
plot(deltaValues, averageTimes, '-o', 'LineWidth', 2);
xlabel('\delta');
ylabel('Average elapsed time (s)');
title('Average elapsed time as a function of \delta');
grid on;