%  Briggs–Rauscher reaction with Ferulic Acid Injection (Modified Plotting)
clear; clc;

% Parameters
p.a = 0.0225; p.b = 0.33; p.h = 0.056; p.q = 0.0015; p.c = 55.5; p.v0 = 0.004;

% Kinetic constants k1-k12
p.k = [5*1e3 2*1e10 3.1*1e12 7.35*1e3 1e5 180.2 3.2*1e4 7.5*1e5 12.5 37 0.013 8*1e5];
p.kFA  = 1e4;   % Scavenging rate constant (Antioxidant activity)
p.kdeg = 0;     % Negligible natural degradation

% Phase 1: Establish Baseline Oscillations (0 to 1200s)
tspan1 = [0 2000];
% y0: [I-, HOI, HIO2, I2, IO2, Mn2+, HO2, O2, IMA, FA]
y0_start = [0.0004 0 0 0 1e-5 p.v0 10^(-3) 1e-5 1e-5 0]; 
opts = odeset('RelTol',1e-9, 'AbsTol',1e-12);
[t1, y1] = ode15s(@sysode_FA, tspan1, y0_start, opts, p);

% Phase 2: Inject Ferulic Acid at t=1200s
y0_inj = y1(end, :);
y0_inj(10) = 3e-3; % Injection concentration
[t2, y2] = ode15s(@sysode_FA, [2000 30000], y0_inj, opts, p);

% Combine data
t = [t1; t2]; 
y = [y1; y2];

% Plotting HO2 and FA
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% Subplot 1: HO2 (Hydroperoxyl Radical)
subplot(2,1,1)
% Plotting log10 to capture the orders of magnitude in oscillations
plot(t, -log10(abs(y(:,7)) + 1e-18), 'b', 'LineWidth', 1.5)
hold on; 
xline(2000, '--r', 'FA Injection', 'LabelVerticalAlignment', 'bottom');
ylabel('-log_{10}[HO_2]'); 
title('Concentration of Hydroperoxyl Radicals (HO_2)');
grid on;

% Subplot 2: Ferulic Acid (FA)
subplot(2,1,2)
plot(t, y(:,10), 'g', 'LineWidth', 1.5); 
hold on; 
xline(2000, '--r');
ylabel('[Ferulic Acid] (M)'); 
xlabel('Time (s)'); 
title('Consumption of Ferulic Acid');
grid on;
%       I HOI  HIO2  I2  IO2  Mn   HO2     

% ODE Function 
function dy = sysode_FA(~, y, p)
    dy = zeros(10,1);
    % Main BR Reactions (using p.k indices)
    dy(1) = - p.k(1)*p.a*p.h^2*y(1) - p.k(2)*p.h*y(3)*y(1) - p.k(3)*p.h*y(1)*y(2) + p.k(11)*p.c*y(4) + p.k(9)*p.q*y(4) + p.k(10)*p.b*y(2);
    dy(2) = p.k(1)*p.a*p.h^2*y(1) + 2*p.k(2)*p.h*y(3)*y(1) - p.k(3)*p.h*y(1)*y(2) + p.k(11)*p.c*y(4) + p.k(5)*y(3)^2 - p.k(10)*p.b*y(2);
    dy(3) = p.k(1)*p.a*p.h^2*y(1) - p.k(2)*p.h*y(3)*y(1) - p.k(4)*p.h*p.a*y(3) + p.k(12)*p.c*y(5)^2 - 2*p.k(5)*y(3)^2 + p.k(6)*p.c*y(5)*y(6);
    dy(4) = p.k(3)*p.h*y(1)*y(2) - p.k(11)*p.c*y(4) - p.k(9)*p.q*y(4);
    dy(5) = 2*p.k(4)*p.a*p.h*y(3) - 2*p.k(12)*p.c*y(5)^2 - p.k(6)*p.c*y(6)*y(5);
    dy(6) = -p.k(6)*p.c*y(6)*y(5) + p.k(7)*p.b*(p.v0 - y(6));
    
    % Radical Scavenging Term (Reaction with Antioxidant)
    % dy(7) is HO2. Scavenging has a negative sign.
    dy(7) = p.k(7)*p.b*(p.v0 - y(6)) - 2*p.k(8)*y(7)^2 - p.kFA*y(7)*y(10);
    
    dy(8) = p.k(8)*y(7)^2 + p.k(10)*y(2)*p.b; % Oxygen production
    dy(9) = p.k(9)*y(4)*p.q;                 % IMA production
    dy(10) = - p.kFA*y(7)*y(10);             % FA consumption
end

figure;
labels = {"I", "HOI", "HIO_2", "I_2", "IO_2", "Mn", "HO_2"};
idx = [1 2 3 4 5 7];

tiledlayout(1,6, 'TileSpacing','compact', 'Padding','compact')

for k = 1:length(idx)
    i = idx(k);
    nexttile
    plot(t, 1e6*y(:,i), 'LineWidth', 1)
    xline(2000, '--r', 'FA Injection', 'LabelVerticalAlignment', 'bottom');
    title(labels{i})
    grid on
end

sgtitle('10^{6} \times Concentrations over Time')
xlabel('Time')
