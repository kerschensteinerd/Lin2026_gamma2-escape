%% EXTRACT_STIMULUS_TIMING Detect complete stimulus-marker square waves
% Input: deinterleave_channels output <base>_Digit.mat (Digit movie).
% Collapse every frame to its spatial standard deviation. Baseline is
% the single baselineFrame value. Values > baseline+spikeThreshold are
% replaced by their following frame (preceding frame at recording end).
% Detect rising/falling edges after thresholding at baseline+onThreshold;
% discard incomplete boundary events. All indices are 1-based and align
% with the unchanged time dimension of the other channels.
% Output: StimTiming/<base>_StimTiming.mat and a QC PNG. Stim_ON is the
% first ON frame; Stim_OFF the last ON frame. Durations are OFF-ON,
% not inclusive sample counts. Existing outputs are skipped.
% Run after deinterleave_channels, independently of motion correction.
% Next: analyze_spot_size_responses. Inspect every QC plot. See README.md.

clear; clc; close all;

%% ===== DATASET-SPECIFIC SETTINGS =====
baselineFrame     = 200;   % frame at the low baseline (no square wave, no spike)
spikeThreshold    = 60;    % > baseline + this  -> overhead spike
onThreshold       = 5;     % >= baseline + this -> ON
expectedNumEvents = 5;     % square waves per recording (one per spot size); [] to skip the check

%% ===== SELECT INPUT / OUTPUT FOLDERS =====
inputDir = uigetdir(pwd, 'Select folder containing deinterleave_channels output (*_Digit.mat)');
if isequal(inputDir,0)
    error('No input folder selected. Script stopped.');
end

outputDir = fullfile(inputDir, 'StimTiming');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% ===== FIND DIGIT FILES =====
files = dir(fullfile(inputDir,'*_Digit.mat'));
if isempty(files)
    error('No *_Digit.mat files found in input folder.');
end

fprintf('Found %d Digit files.\n', numel(files));

%% ===== PROCESS EACH FILE =====
for f = 1:numel(files)

    fileName = files(f).name;
    baseName = erase(fileName, '_Digit.mat');
    outName  = [baseName '_StimTiming.mat'];

    fprintf('\nProcessing %d/%d: %s\n', f, numel(files), fileName);

    if exist(fullfile(outputDir, outName), 'file')
        fprintf('  %s already done. Skipping.\n', outName);
        continue;
    end

    % ----- LOAD -----
    S = load(fullfile(inputDir, fileName));
    if ~isfield(S,'Digit')
        warning('%s does not contain variable "Digit". Skipping.', fileName);
        continue;
    end

    Digit = double(S.Digit);

    % ----- 1-D STIMULUS TRACE: spatial std of every frame -----
    Stim_raw = squeeze(std(Digit,0,[1 2]));
    Stim_raw = Stim_raw(:);
    N = numel(Stim_raw);

    if N < baselineFrame
        warning('%s has only %d frames (< baselineFrame = %d). Skipping.', fileName, N, baselineFrame);
        continue;
    end

    % ----- BASELINE LEVEL -----
    baseline = Stim_raw(baselineFrame);

    % ----- REMOVE OVERHEAD SPIKES -----
    cleaned = Stim_raw;
    for i = 1:N
        if cleaned(i) > baseline + spikeThreshold
            if i < N
                cleaned(i) = cleaned(i+1);
            else
                cleaned(i) = cleaned(i-1);
            end
        end
    end

    % ----- BINARIZE: square wave ON / OFF -----
    Stim_binary = cleaned >= baseline + onThreshold;

    % ----- FIND ON / OFF EDGES -----
    d = diff(double(Stim_binary));
    Stim_ON  = find(d ==  1) + 1;   % first frame of each square wave
    Stim_OFF = find(d == -1);       % last frame of each square wave
    Stim_ON  = Stim_ON(:);
    Stim_OFF = Stim_OFF(:);

    % ----- KEEP ONLY COMPLETE ON->OFF PAIRS -----
    nBefore = [numel(Stim_ON) numel(Stim_OFF)];
    while ~isempty(Stim_OFF) && (isempty(Stim_ON) || Stim_OFF(1) < Stim_ON(1))
        Stim_OFF(1) = [];           % signal was already ON at the start
    end
    while ~isempty(Stim_ON) && (isempty(Stim_OFF) || Stim_ON(end) > Stim_OFF(end))
        Stim_ON(end) = [];          % signal still ON at the end
    end
    if any(nBefore ~= [numel(Stim_ON) numel(Stim_OFF)])
        warning('%s: dropped an unmatched ON/OFF edge (signal already ON at the start, or still ON at the end). Check the QC figure.', baseName);
    end

    % ----- DURATIONS -----
    stim_durations = Stim_OFF - Stim_ON;
    if isempty(stim_durations)
        stim_duration = NaN;
        warning('%s: no square waves detected. Check baselineFrame / thresholds and the QC figure.', baseName);
    else
        stim_duration = max(stim_durations);
    end

    % ----- EXPECTED COUNT CHECK -----
    if ~isempty(expectedNumEvents) && numel(Stim_ON) ~= expectedNumEvents
        warning('%s: detected %d square waves but expectedNumEvents = %d. Check the QC figure.', ...
            baseName, numel(Stim_ON), expectedNumEvents);
    end

    fprintf('  Detected %d square waves. Duration (OFF-ON) in frames: ', numel(Stim_ON));
    fprintf('%d ', stim_durations);
    fprintf('\n');

    % ----- SAVE -----
    save(fullfile(outputDir, outName), ...
        'Stim_ON','Stim_OFF','stim_duration','stim_durations', ...
        'Stim_raw','Stim_binary','baseline', ...
        'baselineFrame','spikeThreshold','onThreshold');

    % ----- QC FIGURE -----
    fig = figure('Visible','off','Position',[100 100 1000 600]);

    subplot(2,1,1)
    plot(Stim_raw,'k'); hold on
    yline(baseline,'b--','baseline');
    yline(baseline + onThreshold,'r--','ON threshold');
    xlabel('Frame'); ylabel('Digit spatial std (raw)');
    title(strrep(baseName,'_','\_'));
    box off

    subplot(2,1,2)
    plot(cleaned,'k'); hold on
    for i = 1:numel(Stim_ON)
        xline(Stim_ON(i),'g-','LineWidth',1.2);
        xline(Stim_OFF(i),'b-','LineWidth',1.2);
    end
    xlabel('Frame'); ylabel('Cleaned signal (spikes removed)');
    title(sprintf('%d square waves  (green = ON, blue = OFF)', numel(Stim_ON)));
    box off

    saveas(fig, fullfile(outputDir, [baseName '_StimTiming_check.png']));
    close(fig);

    fprintf('  Saved: %s\n', outName);
end

fprintf('\nAll Digit files finished.\n');
