% Ellipse fitting with manual region selection and parameter analysis
clear; clc; close all

%% User input prompt
prompt = {'Number of repetitions for calibration length (default = 3):',
    'Reference length (mm, default = 10):',
    'Number of photos to do an analysis on:'
    };

dlgtitle = 'User Input';
dims = [1 45];
definput = {'4','10', '3'};   % default values

answer = inputdlg(prompt, dlgtitle, dims, definput);

if ~isempty(answer)
    n_number_repetitions = str2double(answer{1});
    per_hoeveelheid_mm = str2double(answer{2});
    n_number_photos = str2double(answer{3});
else
    disp('User cancelled.');
    return
end

%% Select all photos to analyse (in order)
% Make sure the filenames are correct using the indexing system explained
% in the report
[filelist, pathname] = uigetfile( ...
    {'*.jpg;*.jpeg;*.png;*.tif;*.bmp', ...
     'Images (*.jpg, *.jpeg, *.png, *.tif, *.bmp)'}, ...
    'Select ALL photos at once', ...
    'MultiSelect','on');

if isequal(filelist,0)
    error('No photo selected. Script stopped.');
end

if ischar(filelist)
    filelist = {filelist};   % convert to cell array
end

%% Ellipse selectie en info berekenen/opslaan
for rep = 1:n_number_photos
    filename = filelist{rep};
    
     % Filename without extension used as row name
    [~, rowName, ~] = fileparts(filename);
    
    fprintf('\n===== Photo %d/%d: %s =====\n', rep, n_number_photos, rowName);
    
    img_path = fullfile(pathname, filename);
    img = imread(img_path);
    grayy = im2gray(img);
    gray = imadjust(grayy);
    gray = imgaussfilt(gray, 2);
    
    
    %% Calibration
    figure(1)
    imshow(grayy)
    axis on; axis equal;
    title(['Click ', num2str(n_number_repetitions), ' calibration lengths']);
    
    lengte_kalibratie_pixels = zeros(n_number_repetitions,1);
    
    for i = 1:n_number_repetitions
        d_kali = drawline('LineWidth',0.5);
        pos = d_kali.Position;
        lengte_kalibratie_pixels(i) = (sqrt((pos(2,1) - pos(1,1))^2 + (pos(2,2) - pos(1,2))^2))/per_hoeveelheid_mm;    % number of pixels per mm
        delete(d_kali)
    end
    
    % Compute and print calibration information
    lengte_kalibratie_pixels_mean = mean(lengte_kalibratie_pixels);
    lengte_kalibratie_pixels_std = std(lengte_kalibratie_pixels);
    lengte_kalibratie_pixels_min = min(lengte_kalibratie_pixels);
    lengte_kalibratie_pixels_max = max(lengte_kalibratie_pixels);
    
    alpha = 0.05;
    df = n_number_repetitions - 1;
    t_critical = tinv(1 - alpha/2, df); % Voor N=3, df=2, 95% CI

    fprintf('Calibration\n');
    fprintf('Mean Calibration length over %.0f measurements: %.2f ± %.2f pixels per mm\n', ...
        length(lengte_kalibratie_pixels), lengte_kalibratie_pixels_mean, t_critical*lengte_kalibratie_pixels_std/(sqrt(length(lengte_kalibratie_pixels))));
    fprintf('\n');
    fprintf('Min. Calibration length: %.2f pixels per mm\n', lengte_kalibratie_pixels_min);
    fprintf('Max. Calibration length: %.2f pixels per mm\n', lengte_kalibratie_pixels_max);
    close(1)
    
    %% Manually selecting an area
    % EdgeDetection (Canny)
    imshow(gray);
    title('Select an area to crop the image');
    rect = drawrectangle('Color','y','LineWidth',1);
    wait(rect);
    
    % ROI parameters
    roi = rect.Position;
    rect.Visible = 'off';

    % Crop the image
    gray_cropped = imcrop(gray, roi);
    edges = edge(gray_cropped, 'Canny', [0.01 0.15], 0.8);

    figure(2)
    imshow(gray_cropped)
    hold on
    [B,~] = bwboundaries(edges, 'noholes');
    for k = 1:length(B)
        boundary = B{k};
        h1 = plot(boundary(:,2), boundary(:,1), 'g.', 'MarkerSize', 5);
    end
    title('Select an elliptical region (double-click to confirm)')
    
    % Select the required edges
    rect = drawrectangle('Color', 'y', 'LineWidth',0.5);
    wait(rect);
    rect.Visible = 'off';
    
    title('Select detected edges that are not needed (double-click to confirm)')
    % Region with pixels that should be excluded
    rect1 = drawfreehand('Color', 'r', 'LineWidth',0.5);
    wait(rect1);
    rect1.Visible = 'off';
    
    % Create mask for detected edges
    mask = createMask(rect);     % Main region
    mask1 = createMask(rect1);   % Region with pixels to be excluded
    edges_roi = edges & mask & ~mask1;
    
    [Y, X] = find(edges_roi);
    x_sel = X;
    y_sel = Y;
    
    h2 = plot(x_sel, y_sel, 'b.', 'MarkerSize', 5);
    legend([h1, h2], {'Edges' , 'Used edges'})
    hold off
    
    %% Ellipsfit
    a = fitellip(x_sel, y_sel);
    v = solveellipse(a);
    a1 = v(1); a2 = v(2); rx0 = v(3); ry0 = v(4); theta = v(5);
    
    %% Compute parameters
    semi_major = max(a1, a2);
    semi_minor = min(a1, a2);
    
    eccentricity = sqrt(1 - (semi_minor^2 / semi_major^2));
    alpha = rad2deg(asin((semi_minor) / (semi_major)));
    inclination_theta = 90 - alpha;
    
    Oppervlakte_pix = pi * semi_major * semi_minor;
    W_L_Verhouding = semi_minor/semi_major;
    D_eq = sqrt(4*semi_minor*semi_major);
    
    % Convert to mm
    semi_major_mm = semi_major / lengte_kalibratie_pixels_mean;
    semi_minor_mm = semi_minor / lengte_kalibratie_pixels_mean;
    
    Oppervlakte_mm = pi * semi_major_mm * semi_minor_mm;
    D_eq_mm = sqrt(4*semi_minor_mm*semi_major_mm);
    
    %% Plot results
    [rx, ry] = drawellip(a, x_sel, y_sel);
    
    figure(3)
    imshow(gray_cropped)
    hold on
    h3 = plot(rx, ry, 'g', 'LineWidth', 1);
    h4 = plot(rx0, ry0, 'r+', 'MarkerSize', 9);
    title('Ellipse fit')
    legend([h3, h4], {'EllipseFit' , 'Ellipse-center'})
    hold off

    %% After all calculations: store results in the table
    results(rep,:) = table( ...
        semi_major, semi_minor, ...
        semi_major_mm, semi_minor_mm, ...
        Oppervlakte_mm, D_eq_mm, ...
        W_L_Verhouding, eccentricity, alpha, inclination_theta, ...
        rx0, ry0, ...
        'VariableNames', { ...
            'SemiMajor_pix','SemiMinor_pix', ...
            'SemiMajor_mm','SemiMinor_mm', ...
            'Area_mm2','DiameterEq_mm', ...
            'WL_Ratio','Eccentricity','ImpactAngle_deg','InclinationTheta_deg', ...
            'CenterX_pix','CenterY_pix' ...
        });
    results.Properties.RowNames{rep} = rowName;

end

disp('Results are saved')

save('ellipse_results_Test.mat', 'results');
load('ellipse_results_Test.mat');