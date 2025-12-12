%% Data analysis
clear; clc; close all

% Definede input parameters
fileNames = {'ellipse_results_KamerTemperatuur.mat', ...
             'ellipse_results_NaturalConvection.mat', ...
             'ellipse_results_Nucleation0.mat', ...
             'ellipse_results_Transition_updated.mat', ...
             'ellipse_results_Film.mat'};
             
dataset_names = ["kamertemp", "natural", "nucleation", "transition", "film"];

true_angles = [90, 90, 90;
               75, 75, 75;
               60, 60, 60;
               30, 30, 30];   

n_repeats = 3;                  % Number of repeats
I_max = size(true_angles,1);    % Number of impact angles
H_max = 3;                      % Number of drop heights
T_max = numel(fileNames);       % Number of datasets (i.e. number of temperatures regimes)

%% Fluid properties and calculate V0 and D0 from mass

% Fluid properties
rho = 1063.56;          % [kg/m^3]
mu = 3.352 * 10^-3;     % [Pa.s]
sigma = 5.241 * 10^-2;  % [N/m]
g = 9.81;               % [m/s^2]

% Mass (in mg)
mass_values_mg = [16, 15.4, 16.4];
mass_mean_mg = mean(mass_values_mg);

% Conversions
mass_mean_kg = mass_mean_mg * 1e-6;
mass_std_kg  = std(mass_values_mg * 1e-6) * 2;

V_error = mass_std_kg / rho;
V_calculated_from_mass = mass_mean_kg / rho;

% Setting pippet minus the error
V_uL = 15 - (V_error * 1e9);                 % volume in microliter
V = V_uL * 1e-9;                             % convert to m^3

% Diameter calculation
D0 = (6*V/pi)^(1/3);                      % diameter in meter
D0_error = (1/3) * (D0 / V) * V_error;          % diameter in meter

% Command window prints
fprintf('\n---- BLOOD DROPLET ANALYSIS ----\n');
fprintf('Measured masses (mg):         %.2f, %.2f, %.2f\n', mass_values_mg);
fprintf('Average mass (mg):            %.3f mg\n', mass_mean_mg);
fprintf('Average mass (kg):            %.3e kg\n', mass_mean_kg);
fprintf('Standard deviation times two (95 percent CI) mass (kg): %.3e kg\n', mass_std_kg);
fprintf('----------------------------------------\n');
fprintf('Volume from mass (m^3):       %.3e m^3\n', V_calculated_from_mass);
fprintf('Volume from mass (µL):        %.3f µL\n', V_calculated_from_mass * 1e9);
fprintf('\n');
fprintf('Volume error margin (m^3):    %.3e m^3\n', V_error);
fprintf('Volume error margin (µL):     %.3f µL\n', V_error * 1e9);
fprintf('\n');
fprintf('Entered volume:               %f µL\n', V_uL);
fprintf('Converted volume:             %.3e m^3\n', V);
fprintf('----------------------------------------\n');

%% LoadAndFilterEllipseData
Results = LoadAndFilterEllipseData(fileNames, dataset_names, I_max, H_max, T_max);

%% Put this data from the loop above into a cell (buildDatasetsWAlphaWL)
[datasets_W_Alpha_WL_matrix, datasets_W_Alpha_WL_vec] = buildDatasetsWAlphaWL(Results, dataset_names);

%% computeHeightAveragedMetrics
[mean_alpha_matrix, mean_WL_matrix, accuracy_matrix, Alpha_SEM_matrix] = computeHeightAveragedMetrics(Results, dataset_names, true_angles, H_max, I_max);

%% t-distribution
alpha = 0.05;
df = n_repeats - 1;
t_kritisch = tinv(1 - alpha/2, df); % N=3, df=2, 95% CI

% W_error is determined inside the errorbar by multiplying with t_critical.
CI_Alpha = Alpha_SEM_matrix * t_kritisch; % Averaged over the different drop heights


%% Colors
% Colors for the matricis
color_90 = [1 0 0];       % Red
color_75 = [1 0.5 0];     % orange
color_60 = [0 0.8 0.8];   % Cyaan
color_30 = [0 0 1];       % Vlue

% Colormatrix
colors_list_matrix = [color_90; color_75; color_60; color_30; ...
                      color_90; color_75; color_60; color_30; ...
                      color_90; color_75; color_60; color_30];

colors_list = [color_90; color_75; color_60; color_30];

% Other colors
colors = [
    0,     0.4470, 0.7410;   % blue
    0.8500, 0.3250, 0.0980;  % orange
    0.9290, 0.6940, 0.1250;  % yellow
    0.4940, 0.1840, 0.5560   % purple
];

%% Plot
%% Drop Impacts at the True Impact Angle
%T_max index values: 1=kamertemp, 2=natural, 3=nucleation, 4=transition, 5=film
%I_max index values: 1=90°, 2=75°, 3=60°, 4=30°

% Plot things
regime_names = ["Room Temp. (23°C)", "Natural Conv. (90°C)", "Nucleation (120°C)", "Transition (160°C)", "Film (230°C)"];
heights = [30, 60, 90];
true_angle_values = true_angles(:,end);
markers = {'s', 'd', 'o'};
colors2  = {'k', 'k', 'k'};

% Plot
for k = 1:I_max

    figure('Position',[100 100 900 700])
    t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for i = 1:T_max
    
        % Select the dataset
        name = dataset_names{i};
        regime_name = regime_names(i);
        raw_alpha_data = Results.(name).raw_alpha;
    
        % Determine the position of the next tile
        tile = nexttile;
        
        hold(tile,'on')
    
        true_angle = true_angle_values(k);
    
        for H = 1:3
            alpha_vals = raw_alpha_data{k, H};
    
            plot(tile,1:n_repeats, alpha_vals, ...
                markers{H}, ...
                'MarkerSize', 10, ...
                'MarkerFaceColor', colors2{H}, ...
                'MarkerEdgeColor', colors2{H}, ...
                'LineStyle', 'none', ...
                'DisplayName', sprintf("H = %d cm", heights(H)));
        end
    
        yline(tile, true_angle, '--k', 'LineWidth', 1.5, ...
              'DisplayName', sprintf("True Angle (%d°)", true_angle))
    
        xlabel(tile,'Repeats')
        ylabel(tile,'Measured Impact Angle (deg)')
        title(tile, sprintf("Spread — %s", regime_name))
        
        xticks(tile, 1:n_repeats)
        xlim(tile, [0.5, n_repeats + 0.5])
        
        grid(tile,'on')
        hold(tile,'off')
    end
    lgd = legend(tile,'Location','northoutside');
    lgd.Layout.Tile = 6;
    sgtitle(sprintf('Drop Impacts at the True Impact Angle - Alpha = %d°', true_angle))
end

%% Accuracy and SEM for Alpha
angle_labels = ["90 deg", "75 deg", "60 deg", "30 deg"];
regime_names = ["Room Temp.(23°C)", "Natural Conv.(90°C)", "Nucleation(120°C)", "Transition(160°C)", "Film(230°C)"];

figure('Position', [100, 100, 800, 800]);

% Accuracy
subplot(2, 1, 1);
bar(accuracy_matrix, 'grouped');

set(gca, 'XTickLabel', regime_names);
xlabel('Temperature Regime');
ylabel('Mean |measured - true| (deg)');
title('Accuracy of Impact Angle Measurement (Averaged over Drop Heights)');
legend(angle_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
grid on;

% Standard Error of the Mean (SEM)
subplot(2, 1, 2);
bar(CI_Alpha, 'grouped');

set(gca, 'XTickLabel', regime_names);
xlabel('Temperature Regime');
ylabel('Mean SEM (deg)');
title('Standard Error of the Mean (Averaged over Drop Heights) — 95% CI, N=3, t_{crit}=4.303');
grid on;

%% Impact hoek met error bars
% Bar plot of the impact angels with 95% CI errorbars
figure('Position',[200 200 900 500]);

% Bar plot
hb = bar(mean_alpha_matrix, 'grouped');
hold on;

% Number of regimes en impact angles
[ngroups, nbars] = size(mean_alpha_matrix);

% X-positions for the bar locations
x = nan(nbars, ngroups);
for b = 1:nbars
    x(b,:) = hb(b).XEndPoints;
end

% Errorbars
for b = 1:nbars
    errorbar(x(b,:), mean_alpha_matrix(:,b), CI_Alpha(:,b), ...
        'k', 'linestyle', 'none', 'LineWidth', 1.3);
end

true_angle_values = [90, 75, 60, 30];

for b = 1:nbars
    yline(true_angle_values(b), '--', 'LineWidth', 1.5, 'Color', colors(b,:));
end

grid on
hold off

%% Width vs temperature per drop height
set(gca, 'XTickLabel', regime_names)
xlabel('Temperature Regime')
ylabel('Measured Impact Angle (deg)')
title('Measured Impact Angles with 95% Confidence Intervals (Averaged over Drop Heights)')
legend(angle_labels, 'Location', 'northoutside', 'Orientation','horizontal')
grid on
hold off

colors_4 = [1 0 0; 1 0.5 0; 0 0.8 0.8; 0 0 1];
angle_labels = {'I = 90°','I = 75°','I = 60°','I = 30°'};
regime_temps = [23 90 120 160 230];
heights = [30, 60, 90];

for H = 1:3
    figure; hold on
    for I = 1:4
        W_vals_T = nan(1, length(regime_temps));
        for T = 1:length(regime_temps)
            idx = (H-1)*4 + I;
            W_vals_T(T) = datasets_W_Alpha_WL_vec{T,1}(idx);
            W_error_T(T) = datasets_W_Alpha_WL_vec{T,4}(idx);
        end
        
        errorbar(regime_temps, W_vals_T * 1000, (W_error_T *t_kritisch) * 1000, '-o', ...
            'Color', colors_4(I,:), ...
            'MarkerFaceColor', colors_4(I,:), ...
            'LineWidth', 1.5, 'MarkerSize', 8)
    end
    xlabel('Temperature (°C)')
    ylabel('Mean Bloodstain Width, $W_{max}$ (mm)', 'Interpreter','latex')
    title(sprintf('W_{max} vs Temperature — H = %d cm', heights(H)))
    legend(angle_labels, 'Location', 'northwest')
    grid on
    hold off
end

%% Model from Laan et al. (2015) plot with exp. data
% Theory curve
A = 1.24; 
X_theory = logspace(-1, 3, 100); 
Y_theory = sqrt(X_theory) ./ (A + sqrt(X_theory));

figure;
h_theory = loglog(X_theory, Y_theory, 'k--', 'LineWidth', 2); 
hold on 

% Heights in meters (H1, H1, H1, H1, H2, H2, H2, H2, H3, H3, H3, H3)
h_exp = [0.30, 0.30, 0.30, 0.30, ...
         0.60, 0.60, 0.60, 0.60, ...
         0.90, 0.90, 0.90, 0.90]; 

% Impact velocity
v_exp = 0.82*sqrt(2 * g * h_exp);

We = (rho .* v_exp.^2 .* D0) ./ sigma;
Re = (rho .* v_exp .* D0) ./ mu;


% Loop the datasets over the measuert width
for i = 1:5
    % Get the width data W
    W_current = datasets_W_Alpha_WL_vec{i, 1};
    Alpha_current = deg2rad(datasets_W_Alpha_WL_vec{i, 2});
    marker_shape = datasets_W_Alpha_WL_vec{i, 6};
    
    % X-axis
    X_base = We .* (Re.^(-0.4)) .* (sin(Alpha_current).^(1.6)); % X-as waarden zijn voor iedereen gelijk

    % Y-axis
    Y_current = (W_current ./ D0) .* (Re.^(-0.2)) .* (sin(Alpha_current).^(-0.2));
    
    % Plot
    scatter(X_base, Y_current, 60, colors_list_matrix, 'filled', ...
            'Marker', marker_shape, 'MarkerEdgeColor', 'k');
end

grid on
xlim([0.5 100])
xlabel('$We Re^{-2/5} \sin^{8/5}(\alpha)$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$W_{max}/D_0 ~ Re^{-1/5} \sin^{-1/5}(\alpha)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Universal Experimental data vs Model', 'FontSize', 12);

% Making the legend
legend_entries = [];
legend_labels  = {};

% Add theory curve
legend_entries(end+1) = h_theory;
legend_labels{end+1}  = 'Theory';

% Add temperature markers (black)
for i = 1:5
    h_temp = plot(nan, nan, datasets_W_Alpha_WL_vec{i,6}, ...
                  'MarkerFaceColor','k', 'MarkerEdgeColor','k');
    legend_entries(end+1) = h_temp;
    legend_labels{end+1}  = datasets_W_Alpha_WL_vec{i,7};   % Temperature label
end

% Add angle markers (coloured)
angle_colors = {color_90, color_75, color_60, color_30};
angle_labels = {'90^{\circ} Impact', '75^{\circ} Impact', ...
                '60^{\circ} Impact', '30^{\circ} Impact'};

for k = 1:4
    h_ang = plot(nan, nan, 'o', 'MarkerFaceColor', angle_colors{k}, ...
                 'MarkerEdgeColor','k');
    legend_entries(end+1) = h_ang;
    legend_labels{end+1}  = angle_labels{k};
end

legend(legend_entries, legend_labels, 'Location','southeast');
