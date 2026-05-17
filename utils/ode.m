a = 0.0225;
b = 0.33;
h = 0.056;
q = 0.0015;
c = 55.5; 
v0 = 0.004;
%k = [1.43*1e3 2*1e10 3.1*1e12 7.3*1e3 6*1e5 180.2 3.2*1e4 5.4*1e5 40 37 0.0396 1.7*1e7];
k = [5*1e3 2*1e10 3.1*1e12 7.35*1e3 1e5 180.2 3.2*1e4 7.5*1e5 12.5 37 0.013 8*1e5];
function dy = sysode7(~,y,a,b,h,q,v0,k,c) 
dy = zeros(9,1);
dy(1) = - k(1)*a*h^2*y(1) - k(2)*h*y(3)*y(1) - k(3)*h*y(1)*y(2) + k(11)*c*y(4) + k(9)*q*y(4) + k(10)*b*y(2);
dy(2) = k(1)*a*h^2*y(1) + 2*k(2)*h*y(3)*y(1) - k(3)*h*y(1)*y(2) + k(11)*c*y(4) + k(5)*y(3)^2 - k(10)*b*y(2);
dy(3) = k(1)*a*h^2*y(1) - k(2)*h*y(3)*y(1) - k(4)*h*a*y(3) + k(12)*c*y(5)^2 - 2*k(5)*y(3)^2 + k(6)*c*y(5)*y(6);
dy(4) = k(3)*h*y(1)*y(2) - k(11)*c*y(4) - k(9)*q*y(4);
dy(5) = 2*k(4)*a*h*y(3) - 2*k(12)*c*y(5)^2 - k(6)*c*y(6)*y(5);
dy(6) = -k(6)*c*y(6)*y(5) + k(7)*b*(v0 - y(6));
dy(7) = k(7)*b*(v0 - y(6)) - 2*k(8)*y(7)^2;
dy(8)=k(8)*y(7)^2+k(10)*y(2)*b;
dy(9)=k(9)*y(4)*q;
end
tspan = [0 2000];
y0=[0.0004 0 0 0 1e-5 v0 10^(-3) 1e-5 1e-5];
opts = odeset('RelTol',1e-10,'AbsTol',1e-16);
[t,y] = ode15s(@(t,y) sysode7(t,y,a,b,h,q,v0,k,c), tspan, y0, opts); 

% Create 7 plots, one for each variable
%We implement the system of 9 odes (including the products O2 and IMA),
%however, for our purposes we plot only the solutions of the first 7
%reactants

loglabels = {"-log_{10}([I^{-}])(M)", "-log_{10}([HOI])(M)", "-log_{10}([HIO_2])(M)", "-log_{10}([I_2])(M)", "-log_{10}([IO_2])(M)", "-log_{10}([Mn^{2+}])(M)", "-log_{10}([HO2])(M)"};
for i = 1:7
    figure
    plot(t, -log10(y(:, i)), 'LineWidth', 1);
    ylabel(loglabels{i})
    xlabel('Time(s)');
    hold on 
    grid on;

end

%%%%%LIMIT CYCLE%%%%%


figure

t_cutoff = 21;
markerlimitb = find(((-log10(y(:,1)))>=5.21 & -log10(y(:,1))<=7.95 & y(:,4)*1e4 >0 &y(:,4)*1e4 <=5.759 & t > t_cutoff)); 
plot(-log10(y(:,1)),y(:,4)*1e4,...
    '- k', ...                    % Line style '-' and color 'k' for black
    'Marker', 'o', ...           % Marker style 'o'
    'Color', 'b', ...            % Marker color (blue)
    'MarkerFaceColor', 'b', ...  % Marker fill color (blue)
    'MarkerEdgeColor', 'b', ...  % Marker edge color (blue)
    'MarkerSize', 2, ...         % Size of the markers
    'MarkerIndices', markerlimitb); 
hold on; % Hold on to add more plots
markerlimitr = find((-log10(y(:,1))>=8.53 & -log10(y(:,1))<=10.04 & y(:,4)*1e4 > 0.046 & y(:,4)*1e4 <=5.93 )); 
plot(-log10(y(:,1)),y(:,4)*1e4,...
    '- k', ...                    % Line style '-' and color 'k' for black
    'Marker', 'o', ...           % Marker style 'o'
    'Color', 'k', ...            % Marker color (red)
    'MarkerFaceColor', 'y', ...  % Marker fill color (red)
    'MarkerEdgeColor', 'y', ...  % Marker edge color (red)
    'MarkerSize', 2, ...         % Size of the markers
    'MarkerIndices', markerlimitr); 
xlabel('-log_{10}([I^-])(M)'), ylabel('[I_2](M)'), grid on;
figure
markerIndices = find((4*t >= 16 & 4*t <= 2743) | (4*t >= 3558 & 4*t <= 6763) | (4*t > 7573)); % Adjusted to include all indices in the range

% Plot blue dots
plot(4*t, -log10(y(:,1)), ...
    '- k', ...                    % Line style '-' and color 'k' for black
    'Marker', 'o', ...           % Marker style 'o'
    'Color', 'b', ...            % Marker color (blue)
    'MarkerFaceColor', 'b', ...  % Marker fill color (blue)
    'MarkerEdgeColor', 'b', ...  % Marker edge color (blue)
    'MarkerSize', 2, ...         % Size of the markers
    'MarkerIndices', markerIndices); % Use the updated marker indices
hold on; % Hold on to add more plots
markerIndices2 = find((4*t >=2760 & 4*t <=3512 ) | (4*t >= 6780 & 4*t <= 7540));
% Plot red dots
plot(4*t, -log10(y(:,1)), ...
    '- k', ...                    % Line style '-' and color 'k' for black
    'Marker', 'o', ...           % Marker style 'o'
    'Color', 'k', ...            % Marker color (red)
    'MarkerFaceColor', 'y', ...  % Marker fill color (red)
    'MarkerEdgeColor', 'y', ...  % Marker edge color (red)
    'MarkerSize', 2, ...         % Size of the markers
    'MarkerIndices', markerIndices2); % Use the same marker indices
xlabel('Time(s)'), ylabel('-log_{10}([I^-])(M)')
grid on

figure
% ====== PLOT ======
hI = plot(t, -log10(y(:,1)), 'LineWidth', 1);

slowColorB = [0.1 0.3 0.7]; %Blue
hold on

yl = ylim;
% Adjust the y-limits to ensure the patches cover the entire range
ylim([yl(1) yl(2)])

% ====== PATCH REALI (sfondo) ======
patch([0 675 675 0],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([675 700 700 675],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([700 876 876 700],[yl(1) yl(1) yl(2) yl(2)], ...
      [1 1 0],'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([876 903 903 876],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([903 1680 1680 903],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1680 1705 1705 1680],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1705 1879 1879 1705],[yl(1) yl(1) yl(2) yl(2)], ...
      [1 1 0],'FaceAlpha',0.5,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1879 1906 1906 1879],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1906 2000 2000 1906],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.5,'EdgeColor','none', ...
      'HandleVisibility','off')

% ====== PATCH FITTIZI PER LEGENDA ======
hSlowB = patch(NaN,NaN,slowColorB,'FaceAlpha',0.5,'EdgeColor','none');
hSlowY = patch(NaN,NaN,[1 1 0],'FaceAlpha',0.5,'EdgeColor','none');
hFast = patch(NaN,NaN,[0 1 0],'FaceAlpha',0.1,'EdgeColor','none');
% ====== LEGENDA ======
legend([hI hSlowB hSlowY hFast], {'-log_{10}([I^{-}])','slow(Radical)', 'slow(Non-Radical)', 'fast'}, 'Location','best')
xlabel('Time(s)');

hold off
grid on


%====== PLOT IO2 ======
figure
hIO2 = plot(t, -log10(y(:,5)), 'LineWidth', 1);

slowColorB = [0.1 0.3 0.7]; %Blue
hold on

yl = ylim;
%Adjust the y-limits to ensure the patches cover the entire range
ylim([yl(1) yl(2)])

%====== PATCH REALI (sfondo) ======
patch([0 675 675 0],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([675 700 700 675],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([700 876 876 700],[yl(1) yl(1) yl(2) yl(2)], ...
      [1 1 0],'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([876 903 903 876],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([903 1680 1680 903],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1680 1705 1705 1680],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1705 1879 1879 1705],[yl(1) yl(1) yl(2) yl(2)], ...
      [1 1 0],'FaceAlpha',0.5,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1879 1906 1906 1879],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1906 2000 2000 1906],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.5,'EdgeColor','none', ...
      'HandleVisibility','off')

%====== PATCH FITTIZI PER LEGENDA ======
hSlowB = patch(NaN,NaN,slowColorB,'FaceAlpha',0.5,'EdgeColor','none');
hSlowY = patch(NaN,NaN,[1 1 0],'FaceAlpha',0.5,'EdgeColor','none');
hFast = patch(NaN,NaN,[0 1 0],'FaceAlpha',0.1,'EdgeColor','none');
%====== LEGENDA ======
legend([hIO2 hSlowB hSlowY hFast], {'-log_{10}([IO_2])','slow(Radical)', 'slow(Non-Radical)', 'fast'}, 'Location','best')
xlabel('Time(s)');
hold off
grid on

%====== PLOT HO2======
figure
hHO2 = plot(t, -log10(y(:,7)), 'LineWidth', 1);

slowColorB = [0.1 0.3 0.7]; %Blue
hold on

yl = ylim;
%Adjust the y-limits to ensure the patches cover the entire range
ylim([yl(1) yl(2)])

%====== PATCH REALI (sfondo) ======
patch([0 675 675 0],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([675 700 700 675],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([700 876 876 700],[yl(1) yl(1) yl(2) yl(2)], ...
      [1 1 0],'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')
patch([876 903 903 876],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([903 1680 1680 903],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.50,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1680 1705 1705 1680],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1705 1879 1879 1705],[yl(1) yl(1) yl(2) yl(2)], ...
      [1 1 0],'FaceAlpha',0.5,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1879 1906 1906 1879],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 1 0],'FaceAlpha',0.1,'EdgeColor','none', ...
      'HandleVisibility','off')

patch([1906 2000 2000 1906],[yl(1) yl(1) yl(2) yl(2)], ...
      slowColorB,'FaceAlpha',0.5,'EdgeColor','none', ...
      'HandleVisibility','off')

%====== PATCH FITTIZI PER LEGENDA ======
hSlowB = patch(NaN,NaN,slowColorB,'FaceAlpha',0.5,'EdgeColor','none');
hSlowY = patch(NaN,NaN,[1 1 0],'FaceAlpha',0.5,'EdgeColor','none');
hFast = patch(NaN,NaN,[0 1 0],'FaceAlpha',0.1,'EdgeColor','none');
%====== LEGENDA ======
legend([hHO2 hSlowB hSlowY hFast], {'-log_{10}([HO_2])','slow(Radical)', 'slow(Non-Radical)', 'fast'}, 'Location','best')

hold off
xlabel('Time(s)');
grid on

