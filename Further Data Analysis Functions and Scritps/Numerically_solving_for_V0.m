%% Numerically solving for the impact velocity
clear; clc; close all

% Definede input parameters
fileNames = {'ellipse_results_KamerTemperatuur.mat', ...
             'ellipse_results_NaturalConvection.mat', ...
             'ellipse_results_Nucleation.mat', ...
             'ellipse_results_Transition.mat', ...
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

%% Loading the data and use self-written functions
% LoadAndFilterEllipseData
Results = LoadAndFilterEllipseData(fileNames, dataset_names, I_max, H_max, T_max);

% Put this data from the loop above into a cell (buildDatasetsWAlphaWL)
[datasets_W_Alpha_WL_matrix, datasets_W_Alpha_WL_vec] = buildDatasetsWAlphaWL(Results, dataset_names);

% computeHeightAveragedMetrics
[mean_alpha_matrix, mean_WL_matrix, accuracy_matrix, Alpha_SEM_matrix] = computeHeightAveragedMetrics(Results, dataset_names, true_angles, H_max, I_max);

%% t-distribution
alpha = 0.05;
df = n_repeats - 1;
t_kritisch = tinv(1 - alpha/2, df); % N=3, df=2, 95% CI

% W_error is determined inside the errorbar by multiplying with t_critical.
CI_Alpha = Alpha_SEM_matrix * t_kritisch; % Averaged over the different drop heights

%% Numerically solving for the Impact Velocity

% Fitting parameter
A = 1.24;
A2 = 1.7277;
%A2 = 1.2401;

C2 = 0.8035; % A = 1.2401

% Empty arrays
W_all = [];
Alpha_all = [];

V0_results = zeros(size(W_all));
V0_guess = 3; 

V0_all = [];
V0_per_dataset = cell(5,1);

% Extract all data for the width and impact angle and put it all in one
% column vector
for i = 1:5
    W_current = datasets_W_Alpha_WL_vec{i, 1};
    Alpha_current = deg2rad(datasets_W_Alpha_WL_vec{i, 2});
    
    V0_results = zeros(numel(W_current),1);

    %W_all = [W_all;W_current(:)];
    %Alpha_all = [Alpha_all;Alpha_current(:)];

    for j = 1:numel(W_current)
        W = W_current(j);
        alpha = Alpha_current(j);
    
        f = @(V0) ( ...
            ( (rho * V0 * D0 / mu)^(1/5) ) * ...
            ( ( ((rho * V0^2 * D0 / sigma) * (rho * V0 * D0 / mu)^(-2/5))^0.5 * sin(alpha) ) / ...
            ( A + ((rho * V0^2 * D0 / sigma) * (rho * V0 * D0 / mu)^(-2/5))^0.5 * (sin(alpha))^(4/5) ) ) ...
            ) - (W / D0);
    
        V0_results(j) = fzero(f, V0_guess);
    end

    V0_all = [V0_all; V0_results];

    V0_per_dataset{i} = V0_results;
end

% Plot a histogram of the results
bins = 0.174;

figure(1)
histogram(V0_all, 'BinWidth', bins, ...
    'FaceColor',[0.2 0.6 0.8], ...
    'EdgeColor','k');

xlabel('Numerically calculated impact velocity $V_0$ (m/s)', ...
       'Interpreter','latex')
ylabel('Counts')
title('Distribution of numerically calculated impact velocity')
grid on
hold on

% Heights in meters
h_exp2 = [0.3, 0.6, 0.9];

% Theoretische impact velocities
v_exp2 = C2*sqrt(2 * g * h_exp2);
%v_exp2 = sqrt(2 * g * h_exp2);

% xlines for the calculated velocities with v0(h) = Cv.sqrt(2gh)
for i = 1:length(v_exp2)
    xline(v_exp2(i), '--r', 'Label', ['H = ' num2str(h_exp2(i)) ' m'], 'LineWidth', 1.5);
end
hold off