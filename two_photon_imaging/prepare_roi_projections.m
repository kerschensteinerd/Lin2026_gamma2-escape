%% PREPARE_ROI_PROJECTIONS Filter registered movies and create selection images
% Input: register_channels outputs, including registered_Gcamp and
% optional registered_<SECOND_CHANNEL_NAME> in matching files.
% Crop yUse/xUse. When APPLY_HIGHPASS_FILTER is true, subtract a 50-frame
% moving mean from each GCaMP pixel, then spatially median-filter every
% frame (3 x 3, symmetric padding). The second channel is median-filtered
% without the temporal operation. Edit all settings for the acquisition.
% Output: ROI_detection/<base>_ROIprep.mat contains Fcrop (unfiltered),
% Ffilt (filtered), stdImg = std(Ffilt,0,3), meanImg = mean(Fcrop,3),
% and meanImg_Red (mean of median-filtered second channel), if available.
% Requires Image Processing Toolbox. The displayed QC projections are
% not saved. Next: segment_tdtomato_dendrites. See README.md.

clear; clc; close all;

%% ===== DATASET-SPECIFIC SETTINGS =====
% Fixed "safe" region to keep after registration, excluding shifting/
% wrapped borders introduced by frame-by-frame registration (register_channels).
yUse = 20:110;
xUse = 20:110;

% High-pass filter (GCaMP only): moving-average window (frames)
% subtracted from each pixel's time trace to remove slow drift/bleaching.
movWin = 50;

% 2-D spatial median filter kernel size [rows cols], applied per frame.
medFiltSize = [3 3];

% Whether to subtract the moving mean before spatial filtering.
% This changes both the activity-selection image and exported traces.
APPLY_HIGHPASS_FILTER = true;   % >>> EDIT for the intended acquisition/analysis

% Optional second channel to also process (median filter + mean
% projection only). Must match register_channels's SECOND_CHANNEL_NAME.
% Set to '' if this dataset has no second channel.
SECOND_CHANNEL_NAME = 'Tdtomato';   % >>> EDIT: '' if not applicable

%% ===== SELECT INPUT / OUTPUT FOLDERS =====
inputDir = uigetdir(pwd, 'Select folder containing register_channels output (*_Gcamp_registered.mat)');
if isequal(inputDir,0)
    error('No input folder selected. Script stopped.');
end

outputDir = fullfile(inputDir, 'ROI_detection');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

hasSecondChannel = ~isempty(SECOND_CHANNEL_NAME);

%% ===== FIND FILES =====
files = dir(fullfile(inputDir,'*_Gcamp_registered.mat'));

if isempty(files)
    error('No *_Gcamp_registered.mat files found.');
end

fprintf('Found %d registered files\n',length(files));

%% ===== PROCESS EACH FILE =====
for f = 1:length(files)

    fileName = files(f).name;
    filePath = fullfile(inputDir,fileName);

    fprintf('\nProcessing %d/%d : %s\n',f,length(files),fileName);

    % ----- LOAD GCAMP -----
    S = load(filePath);

    if ~isfield(S,'registered_Gcamp')
        warning('%s missing registered_Gcamp variable',fileName);
        continue
    end

    F = single(S.registered_Gcamp);

    [H,W,T] = size(F);
    fprintf('Movie size: %d x %d x %d\n',H,W,T);

    % ----- REMOVE REGISTRATION BORDER -----
    Fcrop = F(yUse,xUse,:);
    [H2,W2,T2] = size(Fcrop);

    % ----- HIGH PASS FILTER (optional) -----
    if APPLY_HIGHPASS_FILTER
        Fhp = Fcrop - movmean(Fcrop,movWin,3);
    else
        Fhp = Fcrop;
    end

    % ----- MEDIAN FILTER -----
    % 'symmetric' border padding (mirrors nearby real pixels) instead of
    % medfilt2's default zero-padding, which otherwise biases pixels
    % near the frame edge toward zero every frame -- that bias survives
    % averaging (visible as an edge artifact in a mean projection) even
    % though it mostly cancels out in a std projection across time.
    Ffilt = zeros(size(Fhp),'single');
    for t = 1:T2
        Ffilt(:,:,t) = medfilt2(Fhp(:,:,t),medFiltSize,'symmetric');
    end

    % ----- STD / MEAN PROJECTIONS (GCaMP) -----
    stdImg  = std(Ffilt,0,3);
    meanImg = mean(Fcrop,3);

    baseName = erase(fileName,'_Gcamp_registered.mat');

    % ----- SECOND CHANNEL (optional): median filter + mean projection only -----
    hasSecondForThisFile = hasSecondChannel;
    if hasSecondForThisFile
        secondFile = fullfile(inputDir, [baseName '_' SECOND_CHANNEL_NAME '_registered.mat']);
        if exist(secondFile,'file')
            S2 = load(secondFile);
            varName = ['registered_' SECOND_CHANNEL_NAME];
            if isfield(S2, varName)
                F2 = single(S2.(varName));
                F2crop = F2(yUse,xUse,:);

                % 'symmetric' padding -- see the GCaMP median-filter
                % note above; this matters MORE here, since this
                % channel's meanImg_Red is a mean projection, which
                % (unlike a std projection) does not cancel out the
                % border bias that zero-padding would otherwise leave.
                F2filt = zeros(size(F2crop),'single');
                for t = 1:size(F2crop,3)
                    F2filt(:,:,t) = medfilt2(F2crop(:,:,t),medFiltSize,'symmetric');
                end
                meanImg_Red = mean(F2filt,3);   % named generically "_Red" -- see header
            else
                warning('%s missing variable %s. Proceeding without second channel for this recording.', secondFile, varName);
                hasSecondForThisFile = false;
            end
        else
            hasSecondForThisFile = false;
        end
    end

    % ----- DISPLAY QUICK CHECK -----
    figure(1); clf

    subplot(1,3,1)
    imagesc(meanImg)
    axis image off
    colormap gray
    title('Mean (GCaMP)')

    subplot(1,3,2)
    imagesc(stdImg)
    axis image off
    colormap gray
    title('STD (GCaMP)')

    subplot(1,3,3)
    % Keep QC display valid for flat or very sparse activity images.
    % This changes display limits only; saved fluorescence is unchanged.
    stdUpper = prctile(stdImg(:),99);
    if ~isfinite(stdUpper) || stdUpper <= 0
        stdUpper = max(stdImg(:));
    end
    if ~isfinite(stdUpper) || stdUpper <= 0
        stdUpper = 1;
    end
    imagesc(stdImg,[0 stdUpper])
    axis image off
    colormap gray
    title('STD clipped')

    drawnow

    % ----- SAVE -----
    % Ffilt is the filtered GCaMP movie (high-pass, if
    % APPLY_HIGHPASS_FILTER, then median-filtered) -- same processing
    % used to compute stdImg. Saved so downstream segment_tdtomato_dendrites can
    % extract PER-PIXEL TRACES from this filtered signal (no baseline
    % subtraction applied), rather than from the unfiltered Fcrop.
    if hasSecondForThisFile
        save(fullfile(outputDir,[baseName '_ROIprep.mat']), ...
            'stdImg','meanImg','Fcrop','Ffilt','xUse','yUse','movWin','medFiltSize','APPLY_HIGHPASS_FILTER','SECOND_CHANNEL_NAME','meanImg_Red','-v7.3');
    else
        save(fullfile(outputDir,[baseName '_ROIprep.mat']), ...
            'stdImg','meanImg','Fcrop','Ffilt','xUse','yUse','movWin','medFiltSize','APPLY_HIGHPASS_FILTER','SECOND_CHANNEL_NAME','-v7.3');
    end

    fprintf('Saved: %s\n',[baseName '_ROIprep.mat']);

end

fprintf('\nAll preprocessing finished\n');
