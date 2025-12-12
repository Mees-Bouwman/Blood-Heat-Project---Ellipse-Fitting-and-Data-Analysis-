function Results = LoadAndFilterEllipseData(fileNames, datasetNames, I_max, H_max, T_max)
% All data is structured in a struct for easy calling the data in the
% plotting script.

%% Constants and Index Naming
% I is the impact angle index: 1 = 90°, 2 = 75°, 3 = 60°, 4 = 30°.
% H is the drop height index: 1 = 30 cm, 2 = 60 cm, 3 = 90 cm.

% H_max: Number of drop heights (integer)
% I_max: Number of impact angle conditions (integer)
% T_max: Number of temperature regimes (integer)

%% Example formating fileNames and datasetNames:
% Definede input parameters
%fileNames = {'ellipse_results_KamerTemperatuur.mat', ...
%             'ellipse_results_NaturalConvection.mat', ...
%             'ellipse_results_Nucleation0.mat', ...
%             'ellipse_results_Transition_updated.mat', ...
%             'ellipse_results_Film.mat'};

%dataset_names = ["kamertemp", "natural", "nucleation", "transition", "film"];

%% Load and Prepare Data
data = cell(1, T_max);

    for T = 1:T_max
        % Load the structure 'results' or 'resultsX' from the .mat file X can
        % be any number and is meant to differentiate between different
        % datasets
        loaded_data = load(fileNames{T});
        
        % Find the variable named 'results' or 'resultsX'
        field_names = fieldnames(loaded_data);
        
        data{T} = loaded_data.(field_names{1});
    end

%% Result Structure for loop
% Not all columns (data) in the table are loaded in here, but they can be added in easily.

Results = struct();

    for T = 1:T_max
        res = data{T};
        Rownames = string(res.Properties.RowNames);
        name = datasetNames{T};   % Dataset name
    
        % Empty matrices
        Results.(name).W_mean = zeros(I_max, H_max);
        Results.(name).W_SEM = zeros(I_max, H_max);
        Results.(name).WL_mean = zeros(I_max, H_max);
        Results.(name).Alpha_mean = zeros(I_max, H_max);
        Results.(name).Alpha_SEM = zeros(I_max, H_max);
        
        % Create empty cells for filtered raw data
        Results.(name).raw_alpha = cell(I_max, H_max);
        Results.(name).raw_W = cell(I_max, H_max);
    
        for H = 1:H_max
            for I = 1:I_max
                
                % Mask needed for filtering the data out of the table(s)
                mask = contains(Rownames, "I_" + I + " H_" + H);
    
                % Load filtered data from mask
                W_vals = 2 * res(mask,:).SemiMinor_mm / 1000;
                WL_vals = res(mask,:).WL_Ratio;
                Alpha_vals = res(mask,:).ImpactAngle_deg;
                
                % Store filtered raw data
                % Alpha and Width values
                Results.(name).raw_alpha{I,H} = Alpha_vals;
                Results.(name).raw_W{I,H} = W_vals;
                
                % Means and Standard Error of the Mean
                % Width
                Results.(name).W_mean(I,H) = mean(W_vals);
                Results.(name).W_SEM(I,H) = std(W_vals)/sqrt(numel(W_vals));
                
                % W/L-ratio
                Results.(name).WL_mean(I,H) = mean(WL_vals);
                
                % Alpha impact angle in degrees
                Results.(name).Alpha_mean(I,H) = mean(Alpha_vals);
                Results.(name).Alpha_SEM(I,H) = std(Alpha_vals)/sqrt(numel(Alpha_vals));
            end
        end
    end
end