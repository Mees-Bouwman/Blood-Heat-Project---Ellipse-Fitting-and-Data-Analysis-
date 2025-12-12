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

%% Ellipse selection and parameter calculation/storage
for rep = 1:n_number_photos
    filename = filelist{rep};
    
    % Filename without extension used as row name
    [~, rowName, ~] = fileparts(filename);
    
    fprintf('\n===== Photo %d/%d: %s =====\n', rep, n_number_photos, rowName);
    
    % Covert image for better contrast
    img_path = fullfile(pathname, filename);
    img = imread(img_path);
    gray = im2gray(img);
    gray = histeq(gray);
    
    %% Calibration
    figure(1)
    imshow(gray)
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
    
   %% Ellipse Selection + obtain coordinates
    figure(2)
    imshow(gray)
    axis on; axis equal;
    title('Draw an ellipse around the droplet or stain')
    hold on
    
    h = drawellipse('Color', 'y', 'LineWidth',0.5, 'FaceAlpha', 0);
    wait(h)
    
    hCenter = h.Center;
    hSemiAxis = h.SemiAxes;
    hRotationAngle = h.RotationAngle;
    
    x0 = hCenter(1);
    y0 = hCenter(2);
    a = hSemiAxis(1);
    b = hSemiAxis(2);
    phi = - deg2rad(hRotationAngle);
    
    R = [cos(phi), -sin(phi);
         sin(phi),  cos(phi)];

    % Four corner ellipse points before rotation
    p = [0,  -b;   % left
         a,  0;   % top
         0,  b;   % right
         -a, 0];  % bottom

    % Rotate and translate points
    P = (R * p')' + [x0, y0];
    
    close(2)
    
    % Four coordinates
    P1 = P(1,:); P2 = P(2,:); P3 = P(3,:); P4 = P(4,:);
    
    %% Call function to compute ellipse parameters from four points
    % The four points P1–P4 are row vectors [X, Y].
    [r0, rx, ry, semi_major, semi_minor] = ellipseFourPoints(P1, P2, P3, P4);
    
    eccentricity = sqrt(1 - ((semi_minor^2) / (semi_major^2)));
    alpha = rad2deg(asin((2 * semi_minor) / (2 * semi_major)));
    inclination_theta = 90 - alpha;

    Oppervlakte_ellips = pi * semi_major * semi_minor;
    W_L_Ratio = semi_minor/semi_major;
    D_eq = sqrt(4*semi_minor*semi_major);
    
    % Convert to millimeters using calibration
    semi_major_mm = semi_major / lengte_kalibratie_pixels_mean;
    semi_minor_mm = semi_minor / lengte_kalibratie_pixels_mean;

    Area_mm = pi * semi_minor_mm * semi_major_mm;
    D_eq_mm = sqrt(4*semi_minor_mm*semi_major_mm);
    
    %% Plot everything
    figure(3)
    imshow(img)
    axis on; axis equal; hold on
    
    h1 = plot(rx, ry, 'g', 'LineWidth', 1);
    h2 = plot(P(:,1), P(:,2), 'ko', 'MarkerSize',5, 'MarkerFaceColor','r');
    h3 = plot(transpose(r0(1)), transpose(r0(2)), 'r+', 'MarkerSize',9);
    
    title('Fitting an ellipse')
    legend([h1, h2, h3], {'EllipseFit', 'Four selected points' , 'Ellipse-center'})
    hold off

    %% After all calculations: store results in the table
    results(rep,:) = table( ...
        semi_major, semi_minor, ...
        semi_major_mm, semi_minor_mm, ...
        Area_mm, D_eq_mm, ...
        W_L_Ratio, eccentricity, alpha, inclination_theta, ...
        r0(1), r0(2), ...
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

save('ellipse_results_test.mat', 'results');
load('ellipse_results_test.mat');