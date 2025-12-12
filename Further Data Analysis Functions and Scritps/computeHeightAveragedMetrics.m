function [mean_alpha_matrix, mean_WL_matrix, accuracy_matrix, Alpha_SEM_matrix] = ...
         computeHeightAveragedMetrics(Results, dataset_names, true_angles, H_max, I_max)
% This function computes height averaged ellipse metrics from preprocessed
% data stored in the Results structure.

% Inputs:
% Results: Structure containing ellipse results per dataset
% dataset_names: Cell array with dataset identifiers
% true_angles: True (expected) impact angles per index (matrix in the form of n x m where n is the amount of repeats and m is the amount of impact angles)
% Example:
% true_angles = [90, 90, 90;
%               75, 75, 75;
%               60, 60, 60;
%               30, 30, 30];   

% H_max: Number of drop heights (integer)
% I_max: Number of impact angle conditions (integer)

% Outputs:
% mean_alpha_matrix: Averaged Alpha (impact angle) over all heights
% mean_WL_matrix: Averaged Width/Length ratio over heights
% accuracy_matrix: Deviation between measured Alpha and true angles
% Alpha_SEM_matrix: SEM of Alpha, combined across heights

numDatasets = length(dataset_names);

% Initialise matrices
mean_alpha_matrix = zeros(numDatasets, I_max);
mean_WL_matrix    = zeros(numDatasets, I_max);
accuracy_matrix   = zeros(numDatasets, I_max);
Alpha_SEM_matrix  = zeros(numDatasets, I_max);

    %% Process each dataset
    for d = 1:numDatasets
        name = dataset_names{d};

        % Accuracy alpha
        mean_devi = abs(Results.(name).Alpha_mean - true_angles);
        Results.(name).accuracy_overH = mean(mean_devi, 2);

        % Average over height (Alpha, WL)
        Results.(name).Alpha_mean_overH = mean(Results.(name).Alpha_mean, 2);
        Results.(name).WL_mean_overH    = mean(Results.(name).WL_mean, 2);

        % Aggerated SEM
        SEM_squared = Results.(name).Alpha_SEM.^2;
        sum_SEM_squared = sum(SEM_squared, 2);
        Results.(name).Alpha_SEM_overH = sqrt(sum_SEM_squared) / H_max;

        % Fill matrices
        mean_alpha_matrix(d,:) = Results.(name).Alpha_mean_overH;
        mean_WL_matrix(d,:)    = Results.(name).WL_mean_overH;
        accuracy_matrix(d,:)   = Results.(name).accuracy_overH;
        Alpha_SEM_matrix(d,:)  = Results.(name).Alpha_SEM_overH;
    end
end
