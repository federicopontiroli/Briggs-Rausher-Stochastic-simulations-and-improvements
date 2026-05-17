clear;
close all;
%BRModelDefinitionNoWater;
BRModelReducedwithoutIMA;
%[T,Dynamics]=simOptimized_disc(vMinusBRr,vPlusBRr,cBRr,initialStateBRr,2e4,1e-2);
%[T,Dynamics]=simRSSAOptimized(vMinusBR,vPlusBR,cBR,initialStateBR,0.09,5000,1e-3);
%[T,Dynamics]=simRSSA(vMinusBR,vPlusBR,cBR,initialStateBR,0.05,1000,1e-3);
%[T,Dynamics]=simFRM(vMinusBR,vPlusBR,cBR,initialStateBR,1e4,1e-3);
%[T,Dynamics]=simDM_optimized_disc(vMinusBR,vPlusBR,cBR,initialStateBR,1e3,1e-3);
%[T,Dynamics]=simNRM(vMinusBR,vPlusBR,cBR,initialStateBR,1e6,1e-1);

numSimulations = 50; % Number of simulations to run
performanceTimes = zeros(numSimulations, 1); % Preallocate array for performance times

parfor i = 1:numSimulations
    tic; % Start timer
    [T, Dynamics] = simRSSAOptimized(vMinusBRr, vPlusBRr, cBRr, initialStateBRr,0.1, 2e2, 1e-2);
    performanceTimes(i) = toc; % Store elapsed time
end

averagePerformanceTime = mean(performanceTimes); % Calculate average performance time
disp(['Average Performance Time: ', num2str(averagePerformanceTime), ' seconds']);