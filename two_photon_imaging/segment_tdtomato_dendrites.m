%% SEGMENT_TDTOMATO_DENDRITES Select active labeled pixels and export traces
% Input: prepare_roi_projections output (*_ROIprep.mat), requiring
% stdImg, Ffilt and meanImg_Red. Recordings without meanImg_Red are skipped.
% Choose red mean-image and green std-image thresholds with thresholdGUI.
% The final mask is the intersection of the two masks. Traces are taken
% from Ffilt without further baseline subtraction or normalization.
% Output: Segmentation/<base>_Segmentation.mat includes masks, thresholds,
% zip_RedGreen / zip_Green [row col], and RedGreen_pixel_resp /
% Green_pixel_resp [Frames x nPixels], with matching pixel-column order.
% Mask PNGs are saved in Segmentation/Figures/<base>/. Existing result
% files are skipped; remove them before intentionally reprocessing.
% Run as a script from the MATLAB base workspace; thresholdGUI writes
% the accepted value there. Closing its window cancels that recording.
% Next: analyze_spot_size_responses. See README.md for limits and settings.

clear; clc; close all;

%% ===== SELECT INPUT / OUTPUT FOLDERS =====
inputDir = uigetdir(pwd, 'Select folder containing prepare_roi_projections output (*_ROIprep.mat)');
if isequal(inputDir,0)
    error('No input folder selected. Script stopped.');
end

outputDir = fullfile(inputDir, 'Segmentation');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
figBaseDir = fullfile(outputDir, 'Figures');
if ~exist(figBaseDir,'dir')
    mkdir(figBaseDir);
end

%% ===== FIND FILES =====
files = dir(fullfile(inputDir, '*_ROIprep.mat'));
if isempty(files)
    error('No *_ROIprep.mat files found in inputDir.');
end

fprintf('Found %d ROIprep files.\n', numel(files));

%% ===== PROCESS EACH FILE =====
for f = 1:numel(files)

    fileName = files(f).name;
    baseName = erase(fileName, '_ROIprep.mat');
    SegName = [baseName '_Segmentation.mat'];

    if exist(fullfile(outputDir, SegName), 'file')
        fprintf('\n%s already segmented. Skipping.\n', baseName);
        continue;
    end

    S = load(fullfile(inputDir, fileName));

    if ~isfield(S,'meanImg_Red')
        fprintf('\n%s has no second/tdTomato channel (meanImg_Red missing) -- skipping.\n', baseName);
        continue;
    end
    if ~isfield(S,'stdImg') || ~isfield(S,'Ffilt')
        warning('%s missing stdImg/Ffilt. Skipping.', fileName);
        continue;
    end

    Mean_Red  = S.meanImg_Red;
    Std_Green = S.stdImg;
    Ffilt     = S.Ffilt;   % filtered (high-pass [if enabled] + median-filtered) GCaMP movie -- used for trace extraction
    if isfield(S,'xUse'); xUse = S.xUse; else; xUse = []; end
    if isfield(S,'yUse'); yUse = S.yUse; else; yUse = []; end

    fprintf('\nProcessing %d/%d: %s\n', f, numel(files), baseName);

    %% ----- THRESHOLD RED (anatomy) -----
    close all;
    thresholdGUI(Mean_Red);
    if exist('threshold', 'var')
        Red_Threshold = threshold;
        clear threshold;
        disp(['Using red threshold: ' num2str(Red_Threshold)]);
    else
        warning('Red threshold not set for %s (GUI closed without saving). Skipping.', baseName);
        continue;
    end

    Mask_mean_Red = Mean_Red;
    Mask_mean_Red(Mask_mean_Red<Red_Threshold) = 0;
    Mask_mean_Red(Mask_mean_Red ~= 0) = 1;

    [row_R, col_R] = find(Mask_mean_Red == 1);
    zip_Red = [row_R(:),col_R(:)];   % [row col] of all red-masked pixels

    close all

    %% ----- THRESHOLD GREEN (activity) -----
    figure; imshow(Mask_mean_Red);   % show red mask for visual reference while thresholding green
    thresholdGUI(Std_Green);
    if exist('threshold', 'var')
        Green_Threshold = threshold;
        clear threshold;
        disp(['Using green threshold: ' num2str(Green_Threshold)]);
    else
        warning('Green threshold not set for %s (GUI closed without saving). Skipping.', baseName);
        continue;
    end

    mask_std_Green = Std_Green;
    mask_std_Green(mask_std_Green<Green_Threshold) = 0;
    mask_std_Green(mask_std_Green~=0) = 1;

    [row, col] = find(mask_std_Green == 1);
    zip_Green = [row(:),col(:)];   % [row col] of all green-masked pixels

    close all;

    %% ----- EXTRACT FILTERED TRACES FOR ALL GREEN-MASKED PIXELS -----
    nGreenPix = size(zip_Green,1);
    nFrames = size(Ffilt,3);
    Green_pixel_resp = zeros(nFrames, nGreenPix, 'single');
    for p = 1:nGreenPix
        Green_pixel_resp(:,p) = squeeze(Ffilt(zip_Green(p,1), zip_Green(p,2), :));
    end

    %% ----- COMBINE: KEEP ONLY RED PIXELS THAT ARE ALSO GREEN-ACTIVE, EXTRACT FILTERED TRACES -----
    RedGreen_pixel_resp = [];
    x_red = []; y_red = [];
    for k = 1:size(zip_Red,1)
        if mask_std_Green(row_R(k), col_R(k))>0
            x_red = [x_red,row_R(k)]; y_red = [y_red,col_R(k)];

            % Filtered GCaMP trace for this pixel (no baseline subtraction).
            RedGreen_pixel_resp = [RedGreen_pixel_resp, squeeze(Ffilt(row_R(k), col_R(k),:))];
        end
    end
    x_red = x_red'; y_red = y_red';
    zip_RedGreen = [x_red(:),y_red(:)];   % [row col] of red AND green-active pixels; column order matches RedGreen_pixel_resp

    % Final mask: pixel must belong to the labeled cell's anatomy (red)
    % AND show Ca2+ activity (green).
    Mask = Mask_mean_Red & mask_std_Green;

    %% ----- SAVE QC FIGURES -----
    figDir = fullfile(figBaseDir, baseName);
    if ~exist(figDir,'dir'); mkdir(figDir); end

    figure;imagesc(Mask_mean_Red);colormap hot;
    saveas(gcf,fullfile(figDir,'maskRED'),'png');

    figure;imagesc(mask_std_Green);colormap hot;
    saveas(gcf,fullfile(figDir,'maskGreen'),'png');

    figure;imagesc(Mask);colormap hot;
    saveas(gcf,fullfile(figDir,'finalmask'),'png');

    close all;

    %% ----- SAVE RESULTS -----
    save(fullfile(outputDir, SegName), ...
        'zip_RedGreen', 'RedGreen_pixel_resp', 'zip_Green', 'Green_pixel_resp', ...
        'Red_Threshold','Green_Threshold', ...
        'Mask_mean_Red','mask_std_Green','Mask','xUse','yUse');

    fprintf('Saved: %s\n', SegName);
end

fprintf('\nAll files finished.\n');
