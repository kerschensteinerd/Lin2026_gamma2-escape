%% REGISTER_CHANNELS Apply anatomical-channel translations to matched channels
% Input: output of deinterleave_channels (<base>_<channel>.mat).
% The mean of the first min(nRef,nFrames) reference-channel frames is
% the template. Estimate translations on yReg/xReg within the cropped
% image, and apply the same transform to every configured channel.
% Output: registered_crop_centerReg/<base>_<channel>_registered.mat,
% with registered_Gcamp, registered_Ref or registered_<second channel>,
% plus shifts [dy dx] and the crop/registration windows.
% Requires Image Processing Toolbox. Confirm equal channel frame counts,
% optical alignment, crop bounds and adequate border margins for shifts.
% Next: prepare_roi_projections. See README.md for configuration.

clear; clc; close all;

%% ===== DATASET-SPECIFIC SETTINGS =====
% Anatomical reference channel used to drive registration.
REFERENCE_CHANNEL_NAME = 'TransmittedLight';   % >>> EDIT: e.g. 'TransmittedLight' or 'Structural'

% Optional second channel to co-register using the same transform.
% Set to '' if this dataset has no second channel.
SECOND_CHANNEL_NAME = 'Tdtomato';              % >>> EDIT: '' if not applicable

% Fixed "safe" region kept after cropping out fixed-position imaging-
% system edge noise (see header above).
xKeep = 8:122;   % x range (2nd array dim) kept after cropping edges
yKeep = 8:122;    % y range (1st array dim) kept after cropping edges;

% Center sub-region (within the CROPPED frame) used ONLY to compute
% registration shifts.
yReg = 20:105;
xReg = 20:105;

nRef = 20;   % number of initial frames averaged to build the reference template

%% ===== SELECT INPUT / OUTPUT FOLDERS =====
inputDir = uigetdir(pwd, 'Select folder containing deinterleave_channels output (per-channel .mat files)');
if isequal(inputDir,0)
    error('No input folder selected. Script stopped.');
end

outputDir = fullfile(inputDir, 'registered_crop_centerReg');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

hasSecondChannel = ~isempty(SECOND_CHANNEL_NAME);

%% ===== FIND FILES =====
ch1_files = dir(fullfile(inputDir, '*_Gcamp.mat'));
refFiles  = dir(fullfile(inputDir, ['*_' REFERENCE_CHANNEL_NAME '.mat']));

if isempty(ch1_files) || isempty(refFiles)
    error('Cannot find *_Gcamp.mat or *_%s.mat in input folder.', REFERENCE_CHANNEL_NAME);
end

fprintf('Found %d Gcamp files and %d %s files.\n', ...
    numel(ch1_files), numel(refFiles), REFERENCE_CHANNEL_NAME);

%% ===== PROCESS EACH RECORDING =====
for f = 1:numel(refFiles)

    ref_name = refFiles(f).name;
    baseName = erase(ref_name, ['_' REFERENCE_CHANNEL_NAME '.mat']);

    ch1_name = [baseName '_Gcamp.mat'];
    ch1_path = fullfile(inputDir, ch1_name);
    ref_path = fullfile(inputDir, ref_name);

    if ~exist(ch1_path, 'file')
        warning('Matching Gcamp file missing for %s. Skipping.', baseName);
        continue;
    end

    thisHasSecond = hasSecondChannel;
    if thisHasSecond
        ch2_name = [baseName '_' SECOND_CHANNEL_NAME '.mat'];
        ch2_path = fullfile(inputDir, ch2_name);
        if ~exist(ch2_path, 'file')
            warning('%s channel missing for %s -- registering without it for this recording.', ...
                SECOND_CHANNEL_NAME, baseName);
            thisHasSecond = false;
        end
    end

    fprintf('\nProcessing %d/%d: %s\n', f, numel(refFiles), baseName);

    % ----- LOAD -----
    S1 = load(ch1_path);
    Sref = load(ref_path);

    if ~isfield(S1, 'Gcamp')
        warning('%s does not contain variable "Gcamp". Skipping.', ch1_name);
        continue;
    end
    if ~isfield(Sref, REFERENCE_CHANNEL_NAME)
        warning('%s does not contain variable "%s". Skipping.', ref_name, REFERENCE_CHANNEL_NAME);
        continue;
    end

    Gcamp = S1.Gcamp;
    RefCh = Sref.(REFERENCE_CHANNEL_NAME);

    if thisHasSecond
        S2 = load(ch2_path);
        if ~isfield(S2, SECOND_CHANNEL_NAME)
            warning('%s does not contain variable "%s". Registering without it.', ch2_name, SECOND_CHANNEL_NAME);
            thisHasSecond = false;
        else
            SecondCh = S2.(SECOND_CHANNEL_NAME);
        end
    end

    % ----- CROP EDGES FIRST -----
    Gcamp = Gcamp(yKeep, xKeep, :);
    RefCh = RefCh(yKeep, xKeep, :);
    if thisHasSecond
        SecondCh = SecondCh(yKeep, xKeep, :);
    end

    [H, W, nFrames] = size(RefCh);
    fprintf('  Cropped size: %d x %d x %d\n', H, W, nFrames);

    if max(yReg) > H || max(xReg) > W
        error('yReg/xReg exceeds cropped image size for %s.', baseName);
    end

    % ----- MAKE REFERENCE TEMPLATE FROM REFERENCE CHANNEL CENTER REGION -----
    nRefUse = min(nRef, nFrames);
    refFull = mean(single(RefCh(:,:,1:nRefUse)), 3);
    refReg  = refFull(yReg, xReg);

    % ----- PREALLOCATE -----
    registered_Gcamp = zeros(H, W, nFrames, class(Gcamp));
    registered_Ref   = zeros(H, W, nFrames, class(RefCh));
    if thisHasSecond
        registered_Second = zeros(H, W, nFrames, class(SecondCh));
    end

    shifts = zeros(nFrames, 2);   % [dy dx] per frame
    Rfixed = imref2d([H W]);

    % ----- REGISTER FRAME BY FRAME -----
    for t = 1:nFrames
        movingFull = single(RefCh(:,:,t));
        movingReg  = movingFull(yReg, xReg);

        tform = imregcorr(movingReg, refReg, 'translation');

        dx = tform.Translation(1);
        dy = tform.Translation(2);
        shifts(t,:) = [dy, dx];

        registered_Ref(:,:,t) = imwarp(RefCh(:,:,t), tform, ...
            'OutputView', Rfixed, 'FillValues', 0);

        registered_Gcamp(:,:,t) = imwarp(Gcamp(:,:,t), tform, ...
            'OutputView', Rfixed, 'FillValues', 0);

        if thisHasSecond
            registered_Second(:,:,t) = imwarp(SecondCh(:,:,t), tform, ...
                'OutputView', Rfixed, 'FillValues', 0);
        end
    end

    % ----- SAVE -----
    save(fullfile(outputDir, [baseName '_Gcamp_registered.mat']), ...
        'registered_Gcamp', 'shifts', 'xKeep', 'yKeep', 'xReg', 'yReg', '-v7.3');

    save(fullfile(outputDir, [baseName '_' REFERENCE_CHANNEL_NAME '_registered.mat']), ...
        'registered_Ref', 'shifts', 'xKeep', 'yKeep', 'xReg', 'yReg', '-v7.3');

    if thisHasSecond
        % Save under a DYNAMICALLY named variable, "registered_<SECOND_CHANNEL_NAME>"
        % (e.g. "registered_Tdtomato"), matching what prepare_roi_projections.m
        % looks for -- NOT the literal name "registered_Second" used above
        % while processing, since prepare_roi_projections does not know that generic name.
        dynVarName = ['registered_' SECOND_CHANNEL_NAME];
        saveStruct = struct();
        saveStruct.(dynVarName) = registered_Second;
        saveStruct.shifts = shifts;
        saveStruct.xKeep  = xKeep;
        saveStruct.yKeep  = yKeep;
        saveStruct.xReg   = xReg;
        saveStruct.yReg   = yReg;

        save(fullfile(outputDir, [baseName '_' SECOND_CHANNEL_NAME '_registered.mat']), ...
            '-struct', 'saveStruct', '-v7.3');
    end

    fprintf('  Saved registered files for %s\n', baseName);

    clear S1 Sref Gcamp RefCh registered_Gcamp registered_Ref shifts
    if exist('S2','var'); clear S2 SecondCh registered_Second; end
end

fprintf('\nAll done.\n');
