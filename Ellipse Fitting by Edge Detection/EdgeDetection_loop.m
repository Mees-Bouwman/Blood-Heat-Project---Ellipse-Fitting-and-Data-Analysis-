% Ellips-fit met handmatige regioselectie en parameteranalyse
clear; clc; close all

%% Gebruiker prompt
prompt = {'Number of repetitions for calibration length (default = 3):',
    'Reference length (mm, default = 10):',
    'Number of photos to do an analysis on:'
    };

dlgtitle = 'User Input';
dims = [1 45];
definput = {'4','10', '3'};   % standaardwaarden

answer = inputdlg(prompt, dlgtitle, dims, definput);

if ~isempty(answer)
    n_aantal_herhalingen = str2double(answer{1});
    per_hoeveelheid_mm = str2double(answer{2});
    n_aantal_fotos = str2double(answer{3});
else
    disp('User cancelled.');
    return
end

%% Alle foto's selecteren die geanalyseerd moeten worden (op volgorde)
% Zorg ervoor dat de namen van de foto's juist zijn
[filelist, pathname] = uigetfile( ...
    {'*.jpg;*.jpeg;*.png;*.tif;*.bmp', ...
     'Afbeeldingen (*.jpg, *.jpeg, *.png, *.tif, *.bmp)'}, ...
    'Selecteer ALLE fotos tegelijk', ...
    'MultiSelect','on');

if isequal(filelist,0)
    error('Geen foto geselecteerd. Script gestopt.');
end

if ischar(filelist)
    filelist = {filelist};   % maak cel-array
end

%% Ellipse selectie en info berekenen/opslaan
for rep = 1:n_aantal_fotos
    filename = filelist{rep};
    img_path = fullfile(pathname, filename);
    
    % Bestandsnaam zonder extensie als ROWNAME
    [~, rowName, ~] = fileparts(filename);
    
    fprintf('\n===== Foto %d/%d: %s =====\n', rep, n_aantal_fotos, rowName);
    
    img_path = fullfile(pathname, filename);
    img = imread(img_path);
    grayy = im2gray(img);
    gray = imadjust(grayy);
    gray = imgaussfilt(gray, 2);
    
    
    %% Kalibratie
    figure(1)
    imshow(grayy)
    axis on; axis equal;
    title(['Klik ', num2str(n_aantal_herhalingen), ' kalibratie lengtes aan']);
    
    lengte_kalibratie_pixels = zeros(n_aantal_herhalingen,1);
    
    for i = 1:n_aantal_herhalingen
        d_kali = drawline('LineWidth',0.5);
        pos = d_kali.Position;
        lengte_kalibratie_pixels(i) = (sqrt((pos(2,1) - pos(1,1))^2 + (pos(2,2) - pos(1,2))^2))/per_hoeveelheid_mm;    % aantal pixels per mm
        delete(d_kali)
    end
    
    % Kalibratielengte berekenen en printen naar de Command Window
    lengte_kalibratie_pixels_mean = mean(lengte_kalibratie_pixels);
    lengte_kalibratie_pixels_std = std(lengte_kalibratie_pixels);
    lengte_kalibratie_pixels_min = min(lengte_kalibratie_pixels);
    lengte_kalibratie_pixels_max = max(lengte_kalibratie_pixels);
    
    fprintf('Kalibratie\n');
    fprintf('Gemiddelde kalibratie lengte over %.0f metingen: %.2f ± %.2f pixels per mm\n', ...
        length(lengte_kalibratie_pixels), lengte_kalibratie_pixels_mean, lengte_kalibratie_pixels_std/(sqrt(length(lengte_kalibratie_pixels))));
    fprintf('\n');
    fprintf('Min. Kalibratie lengte: %.2f pixels per mm\n', lengte_kalibratie_pixels_min);
    fprintf('Max. Kalibratie lengte: %.2f pixels per mm\n', lengte_kalibratie_pixels_max);
    close(1)
    
    %% Handmatige selectie van regio
    % Randdetectie (Canny)
    imshow(gray);
    title('Selecteer een regio om uit te knippen');
    rect = drawrectangle('Color','y','LineWidth',1);
    wait(rect);
    
    % ROI parameters ophalen
    roi = rect.Position;    % [x y width height]
    rect.Visible = 'off';

    % Knip de afbeelding bij
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
    title('Selecteer een ellips regio (dubbelklik om te bevestigen)')
    
    % Selecteer nodige randen
    rect = drawrectangle('Color', 'y', 'LineWidth',0.5);
    wait(rect);
    rect.Visible = 'off';
    
    % Regio met pixels die niet meegenomen moeten worden
    rect1 = drawfreehand('Color', 'r', 'LineWidth',0.5);
    wait(rect1);
    rect1.Visible = 'off';
    
    % Mask maken voor gedetecteerde randen
    mask = createMask(rect);     % Hoofdregio
    mask1 = createMask(rect1);   % Regio met pixels die niet meegenomen moeten worden
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
    
    %% Parameters berekenen
    semi_major = max(a1, a2);
    semi_minor = min(a1, a2);
    
    eccentricity = sqrt(1 - (semi_minor^2 / semi_major^2));
    alpha = rad2deg(asin((semi_minor) / (semi_major)));
    inclination_theta = 90 - alpha;
    
    Oppervlakte_pix = pi * semi_major * semi_minor;
    W_L_Verhouding = semi_minor/semi_major;
    D_eq = sqrt(4*semi_minor*semi_major);
    
    % Omrekenen naar mm
    semi_major_mm = semi_major / lengte_kalibratie_pixels_mean;
    semi_minor_mm = semi_minor / lengte_kalibratie_pixels_mean;
    
    Oppervlakte_mm = pi * semi_major_mm * semi_minor_mm;
    D_eq_mm = sqrt(4*semi_minor_mm*semi_major_mm);
    
    %% Plot resultaat
    [rx, ry] = drawellip(a, x_sel, y_sel);
    
    figure(3)
    imshow(gray_cropped)
    hold on
    h3 = plot(rx, ry, 'g', 'LineWidth', 1);
    h4 = plot(rx0, ry0, 'r+', 'MarkerSize', 9);
    title('Ellipse fit')
    legend([h3, h4], {'EllipseFit' , 'Ellipse-center'})
    hold off

    %% Na alle berekeningen: opslaan in de table
    results41(rep,:) = table( ...
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
    results41.Properties.RowNames{rep} = rowName;

end

disp('Results are saved')

save('ellipse_results_Film_randet.mat', 'results41');
load('ellipse_results_Film_randet.mat');   % laadt de table 'results'
