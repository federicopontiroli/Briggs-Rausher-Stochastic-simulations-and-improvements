
%Selectvariables
%HIO,HIO2,Mn,MnOH
%1e-22 good
%[T,Dynamics]= simDM_optimized(vMinusBR,vPlusBR,cBR,initialStateBR,700);
variableNames = {'I','HOI','HIO_2', 'I_2', 'IO_2', 'Mn', 'HO_2', 'MnOH'};
%varsToPlot = {'I','HOI','HIO_2','I_2','IO_2','HO_2'};
varsToPlot = {'I_2'};
avogadro = 6.02214076e23;
% Trova gli indici corrispondenti
idx = find(ismember(variableNames, varsToPlot));

%colors = turbo(length(idx));
colors = [1 0 0; 0 0 1; 0 1 0; 1 0 1; 0 1 1; 0 0 0]; % Define colors: red and blue
figure; hold on
for k = 1:length(idx)
    i = idx(k);
    plot(T, Dynamics(:,i), ...
         'LineWidth', 1, ...
         'Color', colors(k,:));
end
xlabel('Time');
ylabel('Number of molecules');


% Add legend
legend(varsToPlot, 'Location', 'best');
hold off
%xlim([0 3000])


figure
hold on

scatter(-log10(Dynamics(:,1)/1e10), Dynamics(:,4)/1e10, 5, 'k', ...
        'filled', 'MarkerFaceAlpha', 0.05)
xlabel('-log_{10}(#molecules of I/10^{10})');
ylabel('#molecules of I_2/10^{10}');
grid on

figure; hold on
for k = 1:length(idx)
    i = idx(k);
    plot(T, -log10(Dynamics(:,i)/1e10), ...
         'LineWidth', 1, ...
         'Color', colors(k,:));
end
xlabel('Time');
ylabel('-log_{10} of the empirical concentration');

% Add legend
legend({'[I]', '[I_2]'}, 'Location', 'best');
hold off