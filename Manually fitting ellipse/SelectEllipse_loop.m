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
    gray = im2gray(img);
    gray = histeq(gray);
    
    %% Kalibratie
    figure(1)
    imshow(gray)
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
    
   %% Ellipse Selection + obtain coordinates
    figure(2)
    imshow(gray)
    axis on; axis equal;
    title('Teken een ellips rond de druppel of vlek')
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
    
    p = [0,  -b;   % links
         a,  0;   % boven
         0,  b;   % rechts
         -a, 0];  % onder
    
    P = (R * p')' + [x0, y0];
    
    close(2)
    
    % Coördinaten vier punten
    P1 = P(1,:); P2 = P(2,:); P3 = P(3,:); P4 = P(4,:);
    
    %% Functie aanroepen om ellipse te plotten
    % De vier coördinaten P1 tot P4 zijn hiervoor nodig, dit zijn rijvectoren
    % [X, Y].
    [r0, rx, ry, semi_major, semi_minor] = ellipseFourPoints(P1, P2, P3, P4);
    
    eccentricity = sqrt(1 - ((semi_minor^2) / (semi_major^2)));
    alpha = rad2deg(asin((2 * semi_minor) / (2 * semi_major)));
    inclination_theta = 90 - alpha;

    Oppervlakte_ellips = pi * semi_major * semi_minor;
    W_L_Verhouding = semi_minor/semi_major;
    D_eq = sqrt(4*semi_minor*semi_major);
    
    % Oppvervlakte en diameter van de ellips
    semi_major_mm = semi_major / lengte_kalibratie_pixels_mean;
    semi_minor_mm = semi_minor / lengte_kalibratie_pixels_mean;

    Oppervlakte_mm = pi * semi_minor_mm * semi_major_mm;
    D_eq_mm = sqrt(4*semi_minor_mm*semi_major_mm);
    
    %% Alles plotten
    figure(3)
    imshow(img)
    axis on; axis equal; hold on
    
    h1 = plot(rx, ry, 'g', 'LineWidth', 1);
    h2 = plot(P(:,1), P(:,2), 'ko', 'MarkerSize',5, 'MarkerFaceColor','r');
    h3 = plot(transpose(r0(1)), transpose(r0(2)), 'r+', 'MarkerSize',9);
    
    title('Fitting an ellipse')
    legend([h1, h2, h3], {'EllipseFit', 'Four selected points' , 'Ellipse-center'})
    hold off

    %% Na alle berekeningen: opslaan in de table
    results(rep,:) = table( ...
        semi_major, semi_minor, ...
        semi_major_mm, semi_minor_mm, ...
        Oppervlakte_mm, D_eq_mm, ...
        W_L_Verhouding, eccentricity, alpha, inclination_theta, ...
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

save('ellipse_results_Kamertemp.mat', 'results');
load('ellipse_results_Kamertemp.mat');