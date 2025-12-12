function [datasets_W_Alpha_WL_matrix, datasets_W_Alpha_WL_vec] = buildDatasetsWAlphaWL(Results, dataset_names)
% This function needs the Results struct obtained from
% LoadAndFilterEllipseData.m and the dataset_names

% Input:
% Results: Structure containing ellipse results per dataset
% dataset_names: Cell array with dataset identifiers

% Output:
% Both are the same data but in a different format
% datasets_W_Alpha_WL_matrix: Stores matrix data in a cell
% datasets_W_Alpha_WL_vec: Stores vetor data in a cell

% Markers and labels
markers = {'o','s','d','^','v'};
labels  = {'Room temp.','Natural conv.','Nucleation','Transition','Film'};

%% MATRIX variant (4x3)
datasets_W_Alpha_WL_matrix = cell(length(dataset_names), 7);

    for d = 1:length(dataset_names)
        name = dataset_names{d};

        % Mean-values
        datasets_W_Alpha_WL_matrix{d,1} = Results.(name).W_mean;      
        datasets_W_Alpha_WL_matrix{d,2} = Results.(name).Alpha_mean;  
        datasets_W_Alpha_WL_matrix{d,3} = Results.(name).WL_mean;     

        % SEM-values
        datasets_W_Alpha_WL_matrix{d,4} = Results.(name).W_SEM;       
        datasets_W_Alpha_WL_matrix{d,5} = Results.(name).Alpha_SEM;   

        % Markers and labels
        datasets_W_Alpha_WL_matrix{d,6} = markers{d};
        datasets_W_Alpha_WL_matrix{d,7} = labels{d};
    end


%% VECTOR VARIANT (1x12)
datasets_W_Alpha_WL_vec = cell(length(dataset_names), 7);

    for d = 1:length(dataset_names)
        name = dataset_names{d};

        % Vectorise (4x3 -> 1x12)
        Results.(name).W_mean_vec     = reshape(Results.(name).W_mean,     1, []);
        Results.(name).Alpha_mean_vec = reshape(Results.(name).Alpha_mean, 1, []);
        Results.(name).WL_mean_vec    = reshape(Results.(name).WL_mean,    1, []);

        Results.(name).Alpha_SEM_vec  = reshape(Results.(name).Alpha_SEM,  1, []);
        Results.(name).W_SEM_vec      = reshape(Results.(name).W_SEM,      1, []);

        % Fill cell-array
        % Mean-values
        datasets_W_Alpha_WL_vec{d,1} = Results.(name).W_mean_vec;
        datasets_W_Alpha_WL_vec{d,2} = Results.(name).Alpha_mean_vec;
        datasets_W_Alpha_WL_vec{d,3} = Results.(name).WL_mean_vec;
        
        % SEM-values
        datasets_W_Alpha_WL_vec{d,4} = Results.(name).W_SEM_vec;
        datasets_W_Alpha_WL_vec{d,5} = Results.(name).Alpha_SEM_vec;
        
        % Markers and labels
        datasets_W_Alpha_WL_vec{d,6} = markers{d};
        datasets_W_Alpha_WL_vec{d,7} = labels{d};
    end
end
