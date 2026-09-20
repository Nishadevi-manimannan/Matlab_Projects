% =========================================================================
%  LUNG NODULE DETECTION 
% =========================================================================
clear; clc; close all;

% Step 1: Image Acquisition & Pre-processing (With Safe Fallback)
fprintf('Opening file selector... (Select your CT image or press Cancel to load demo image)\n');
[file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.dcm', 'Image Files'}, 'Select a CT Scan Slice');

if isequal(file, 0)
    warning('No file selected. Loading built-in MATLAB sample chest CT scan for demonstration...');
    rawImage = imread('chest.png'); 
    grayImage = im2double(rawImage);
else
    fullPath = fullfile(path, file);
    if endsWith(file, '.dcm', 'IgnoreCase', true)
        rawImage = dicomread(fullPath);
        grayImage = double(rawImage);
        grayImage = (grayImage - min(grayImage(:))) / (max(grayImage(:)) - min(grayImage(:)));
    else
        rawImage = imread(fullPath);
        if size(rawImage, 3) == 3
            grayImage = rgb2gray(rawImage);
        else
            grayImage = rawImage;
        end
        grayImage = im2double(grayImage);
    end
end

% Noise reduction using an Edge-Preserving Filter
filteredImage = wiener2(grayImage, [5 5]);

% Step 2: Lung Parenchyma Segmentation (Masking out background/ribcage)
lungThresh = graythresh(filteredImage);
binaryMask = ~imbinarize(filteredImage, lungThresh); 
binaryMask = imclearborder(binaryMask); 
binaryMask = imfill(binaryMask, 'holes'); 
binaryMask = bwareaopen(binaryMask, 2000); 

% Step 3: Isolate High-Density Nodule Candidates Inside Lungs
lungTissue = filteredImage;
lungTissue(~binaryMask) = 0; 

% Local threshold to isolate high-density entities within the air spaces
noduleThresh = mean(lungTissue(binaryMask)) + 1.2 * std(lungTissue(binaryMask));
candidateMask = lungTissue > noduleThresh;
candidateMask = bwareaopen(candidateMask, 15); 

% Step 4: Component Morphological Cleaning
candidateMask = imclose(candidateMask, strel('disk', 2));
candidateMask = imfill(candidateMask, 'holes');

% Step 5: Shape Analysis & Geometric Filtering (Circularity Check)
[labeledImage, numRegions] = bwlabel(candidateMask);
regionProperties = regionprops(labeledImage, 'Area', 'Perimeter', 'Centroid', 'EquivDiameter');

% Architectural Threshold Constraints for Nodules
minArea = 30;         
maxArea = 1200;       
minCircularity = 0.65; 

% Step 6: Initialize Side-by-Side Visualization Display
fig = figure; 
fig.Name = 'Lung Nodule Detection Suite';
fig.NumberTitle = 'off';
fig.Units = 'pixels';
fig.Position; % Clean, fixed size window layout

% Left Panel: Original Image
subplot(1, 2, 1);
imshow(grayImage, []);
title('1. Original CT Scan Image', 'FontSize', 11, 'FontWeight', 'bold');
axis image; 

% Right Panel: Detection Overlay
subplot(1, 2, 2);
imshow(grayImage, []);
title('2. Detected Nodule Overlay (Clean View)', 'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold');
axis image; 
hold on;

detectedCount = 0;

% Initialize Command Window Print Header
fprintf('\n=========================================================\n');
fprintf('        LUNG NODULE DETECTION METRICS REPORT             \n');
fprintf('=========================================================\n');
fprintf('%-10s | %-16s | %-13s | %-10s\n', 'Nodule ID', 'Center (X, Y)', 'Diameter (px)', 'Area (px)');
fprintf('---------------------------------------------------------\n');

for k = 1:numRegions
    area = regionProperties(k).Area;
    perimeter = regionProperties(k).Perimeter;
    centroid = regionProperties(k).Centroid;
    equivDiameter = regionProperties(k).EquivDiameter;
    
    if perimeter > 0
        circularity = (4 * pi * area) / (perimeter ^ 2);
    else
        circularity = 0;
    end
    
    if (area >= minArea) && (area <= maxArea) && (circularity >= minCircularity)
        detectedCount = detectedCount + 1;
        
        % Plotting a perfect mathematical circle around the centroid
        radius = (equivDiameter / 2) + 4; 
        theta = linspace(0, 2*pi, 100);
        xCircle = centroid(1) + radius * cos(theta);
        yCircle = centroid(2) + radius * sin(theta);
        
        % Draw ONLY the green circle layer over the image (NO overlapping text blocks)
        plot(xCircle, yCircle, 'g-', 'LineWidth', 2.0);
        plot(centroid(1), centroid(2), 'g+', 'MarkerSize', 5, 'LineWidth', 1.0);
        
        % Print values cleanly down to the MATLAB command window instead
        fprintf('Nodule #%-2d | (%6.1f, %6.1f) | %11.1f px | %7d px\n', ...
            detectedCount, centroid(1), centroid(2), equivDiameter, area);
    end
end

hold off;
fprintf('=========================================================\n');

if detectedCount == 0
    fprintf('No pulmonary nodule candidates found matching structural criteria.\n');
    msgbox('No pulmonary nodule candidates matching circular geometry parameters were found.', 'Detection Notification');
else
    fprintf('Analysis Completed successfully: Found %d localized nodule candidate(s).\n\n', detectedCount);
end
