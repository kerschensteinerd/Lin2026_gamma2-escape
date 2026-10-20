%% ANALYZE_SPOT_SIZE_RESPONSES Qualify pixels by repeat reliability
% Select the Segmentation folder, StimTiming folder, and a folder with
% *spot_order*.mat files. Use <base>_spot_order.mat when present; a single
% fallback order file applies to all recordings. Each file must contain
% spot_order, listing each size once in presentation order.
% For each square wave, subtract the pixel-wise median over the inclusive
% Stim_ON + baselineWindow(1):baselineWindow(2) interval. The default
% [-34 0] has 35 samples and includes the ON frame. Check that this
% interval is appropriate for the recorded stimulus marker.
% Cut three 67-frame repeat windows starting 33 frames after ON and
% spaced by 67 frames. Smooth each with a fifth-order, 13-frame
% Savitzky-Golay filter. RR is the mean pairwise Pearson correlation
% among repeats. A pixel qualifies if RR > 0.4 for any spot size.
% Average qualified-pixel responses over pixels and repeats for each
% size; these are baseline-subtracted fluorescence, without F0 division.
% Output: Segmentation/SpotAnalysis/<base>_SpotAnalysis.mat and mean-response
% PNGs. Results are overwritten on rerun. Requires Signal Processing
% Toolbox. See README.md for arrays, configuration, and analysis limits.

clear; clc; close all;

%% ===== SETTINGS =====
baselineWindow = [-34 0];   % frames relative to EACH square wave's ON frame used for that wave's baseline (median)
nRepeats       = 3;          % repeats per square wave
repeatOffset   = 33;         % frames after Stim_ON at which the first repeat window starts
repeatPeriod   = 67;         % frames between the starts of consecutive repeat windows
repeatLength   = 67;         % frames in each repeat window
sgolayOrder    = 5;          % Savitzky-Golay polynomial order
sgolayFrameLen = 13;         % Savitzky-Golay window length (odd, > order)
rrThreshold    = 0.4;        % qualify if RR > this for at least one spot size
fps            = 16.67;      % imaging frame rate (Hz); the plot's time axis is in seconds. [] = plot in frames
offWindow      = [2 4];      % seconds from repeat-window start shaded as the stimulus-OFF period (plot only; needs fps)

if mod(sgolayFrameLen,2) == 0 || sgolayFrameLen <= sgolayOrder || sgolayFrameLen > repeatLength
    error('sgolayFrameLen must be odd, larger than sgolayOrder, and not longer than repeatLength.');
end

%% ===== SELECT FOLDERS =====
segDir = uigetdir(pwd, 'Select folder containing segment_tdtomato_dendrites output (*_Segmentation.mat)');
if isequal(segDir,0)
    error('No segmentation folder selected. Script stopped.');
end

stimDir = uigetdir(segDir, 'Select folder containing extract_stimulus_timing output (*_StimTiming.mat)');
if isequal(stimDir,0)
    error('No stimulus-timing folder selected. Script stopped.');
end

orderDir = uigetdir(segDir, 'Select folder containing the spot_order .mat file');
if isequal(orderDir,0)
    error('No spot-order folder selected. Script stopped.');
end

outDir = fullfile(segDir, 'SpotAnalysis');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ===== FIND FILES =====
segFiles = dir(fullfile(segDir, '*_Segmentation.mat'));
if isempty(segFiles)
    error('No *_Segmentation.mat files found in %s', segDir);
end

orderFiles = dir(fullfile(orderDir, '*spot_order*.mat'));
if isempty(orderFiles)
    error('No *spot_order*.mat file found in %s', orderDir);
end

fprintf('Found %d segmentation files.\n', numel(segFiles));

%% ===== PROCESS EACH RECORDING =====
for f = 1:numel(segFiles)

    segName  = segFiles(f).name;
    baseName = erase(segName, '_Segmentation.mat');

    fprintf('\nProcessing %d/%d: %s\n', f, numel(segFiles), baseName);

    % ----- LOAD PIXEL TRACES (segmentation) -----
    S = load(fullfile(segDir, segName), 'RedGreen_pixel_resp', 'zip_RedGreen', 'xUse', 'yUse');
    if ~isfield(S,'RedGreen_pixel_resp') || ~isfield(S,'zip_RedGreen')
        warning('%s is missing RedGreen_pixel_resp / zip_RedGreen. Skipping.', segName);
        continue;
    end
    if isempty(S.RedGreen_pixel_resp)
        warning('%s has no red+green pixels. Skipping.', baseName);
        continue;
    end

    pixelResp = double(S.RedGreen_pixel_resp);    % [Frames x nPixels]
    pixelPos  = S.zip_RedGreen;                   % [nPixels x 2] (row, col)
    if isfield(S,'xUse'); xUse = S.xUse; else; xUse = []; end
    if isfield(S,'yUse'); yUse = S.yUse; else; yUse = []; end

    [nFrames, nPixels] = size(pixelResp);

    % ----- LOAD STIMULUS TIMING (stimulus timing) -----
    stimFile = fullfile(stimDir, [baseName '_StimTiming.mat']);
    if ~exist(stimFile,'file')
        warning('Missing %s_StimTiming.mat. Skipping.', baseName);
        continue;
    end
    T = load(stimFile, 'Stim_ON', 'Stim_OFF');
    if ~isfield(T,'Stim_ON') || ~isfield(T,'Stim_OFF')
        warning('%s_StimTiming.mat is missing Stim_ON / Stim_OFF. Skipping.', baseName);
        continue;
    end
    Stim_ON  = T.Stim_ON(:);
    Stim_OFF = T.Stim_OFF(:);

    % ----- LOAD SPOT ORDER -----
    orderFile = fullfile(orderDir, [baseName '_spot_order.mat']);
    if ~exist(orderFile,'file')
        if numel(orderFiles) == 1
            orderFile = fullfile(orderDir, orderFiles(1).name);
        else
            warning('No %s_spot_order.mat, and the folder has %d *spot_order*.mat files. Skipping %s.', ...
                baseName, numel(orderFiles), baseName);
            continue;
        end
    end
    O = load(orderFile);
    if ~isfield(O,'spot_order')
        warning('%s has no variable "spot_order". Skipping %s.', orderFile, baseName);
        continue;
    end
    spotOrder = O.spot_order(:);

    % ----- CHECKS -----
    if numel(spotOrder) ~= numel(Stim_ON)
        warning('%s: spot_order has %d entries but %d square waves were detected. Skipping.', ...
            baseName, numel(spotOrder), numel(Stim_ON));
        continue;
    end
    if numel(unique(spotOrder)) ~= numel(spotOrder)
        warning('%s: spot_order contains a repeated spot size; this script expects each size once. Skipping.', baseName);
        continue;
    end

    % baseline window of every square wave: must be inside the recording and
    % must start after the previous square wave has ended
    baselineOK = true;
    for i = 1:numel(Stim_ON)
        winStart = Stim_ON(i) + baselineWindow(1);
        winEnd   = Stim_ON(i) + baselineWindow(2);
        if winStart < 1 || winEnd > nFrames
            warning('%s: baseline window of square wave %d falls outside the recording. Skipping.', baseName, i);
            baselineOK = false;
            break;
        end
        if i > 1 && winStart <= Stim_OFF(i-1)
            gapFrames = Stim_ON(i) - Stim_OFF(i-1);
            warning(['%s: baseline window of square wave %d starts at frame %d, but the previous square wave ends at frame %d ' ...
                     '(only %d frames before this ON). baselineWindow(1) must be greater than %d. Skipping.'], ...
                baseName, i, winStart, Stim_OFF(i-1), gapFrames, -gapFrames);
            baselineOK = false;
            break;
        end
    end
    if ~baselineOK
        continue;
    end

    lastFrameNeeded = max(Stim_ON) + repeatOffset + (nRepeats-1)*repeatPeriod + repeatLength - 1;
    if lastFrameNeeded > nFrames
        warning('%s: repeat windows need frame %d but the recording has %d frames. Skipping.', ...
            baseName, lastFrameNeeded, nFrames);
        continue;
    end

    spotSizes = sort(spotOrder);          % ascending order, [nSizes x 1]
    nSizes    = numel(spotSizes);

    % ----- 1) BASELINE OF EVERY PIXEL BEFORE EVERY SQUARE WAVE -----
    % baselineLevel(i,:) = median of each pixel over the window before square wave i
    nStim = numel(Stim_ON);
    baselineLevel = nan(nStim, nPixels);
    for i = 1:nStim
        baselineIdx = Stim_ON(i) + (baselineWindow(1):baselineWindow(2));
        baselineLevel(i,:) = median(pixelResp(baselineIdx,:), 1);
    end

    % ----- 2) REPEAT RESPONSES OF EVERY PIXEL, PER SPOT SIZE -----
    % repeatResp(time, repeat, spot size, pixel)
    repeatResp = nan(repeatLength, nRepeats, nSizes, nPixels);

    for k = 1:nSizes
        stimIdx = find(spotOrder == spotSizes(k));    % which square wave showed this spot size
        onFrame = Stim_ON(stimIdx);

        for r = 1:nRepeats
            firstFrame = onFrame + repeatOffset + (r-1)*repeatPeriod;
            frameIdx   = firstFrame : (firstFrame + repeatLength - 1);

            % subtract the baseline measured before THIS square wave
            respBS = pixelResp(frameIdx,:) - baselineLevel(stimIdx,:);

            % smooth every pixel's window (sgolayfilt works down each column)
            smoothed = sgolayfilt(respBS, sgolayOrder, sgolayFrameLen);
            repeatResp(:, r, k, :) = reshape(smoothed, repeatLength, 1, 1, nPixels);
        end
    end

    % ----- 3) REPEAT RELIABILITY (RR) -----
    % RR = mean pairwise correlation between the repeats
    RR = nan(nPixels, nSizes);
    pairMask = tril(true(nRepeats), -1);      % each pair of repeats counted once

    for p = 1:nPixels
        for k = 1:nSizes
            C = corrcoef(repeatResp(:,:,k,p));   % [nRepeats x nRepeats]
            RR(p,k) = mean(C(pairMask));
        end
    end

    % ----- 4) QUALIFIED PIXELS -----
    isQualified = any(RR > rrThreshold, 2);      % at least one spot size above threshold
    nQualified  = sum(isQualified);

    fprintf('  %d / %d pixels qualified (RR > %.2f for at least one spot size)\n', ...
        nQualified, nPixels, rrThreshold);
    fprintf('  Pixels with RR > %.2f at each spot size (%s): %s\n', rrThreshold, ...
        num2str(spotSizes'), num2str(sum(RR > rrThreshold, 1)));

    % ----- 5) MEAN RESPONSE OF QUALIFIED PIXELS -----
    meanResp = nan(repeatLength, nSizes);
    if nQualified > 0
        for k = 1:nSizes
            meanResp(:,k) = mean(repeatResp(:,:,k,isQualified), [2 4]);
        end
    else
        warning('%s: no pixel passed RR > %.2f. Nothing to plot.', baseName, rrThreshold);
    end

    % ----- SAVE -----
    save(fullfile(outDir, [baseName '_SpotAnalysis.mat']), ...
        'repeatResp','RR','isQualified','pixelPos','meanResp','spotSizes','spotOrder','Stim_ON','Stim_OFF', ...
        'baselineLevel','baselineWindow','nRepeats','repeatOffset','repeatPeriod','repeatLength', ...
        'sgolayOrder','sgolayFrameLen','rrThreshold','fps','offWindow','xUse','yUse');

    % ----- PLOT: mean response of qualified pixels, one panel per spot size (ascending) -----
    if nQualified > 0
        figure(1); clf
        set(gcf, 'Color', 'w', 'Position', [100 100 1600 380]);

        tAxis = 0:(repeatLength-1);
        xLabelText = 'Frame';
        if ~isempty(fps)
            tAxis = tAxis / fps;
            xLabelText = 'Time (s)';
        end
        shadeOff = ~isempty(fps) && ~isempty(offWindow);

        % same y-axis limits in every panel, so responses can be compared across spot sizes
        yLimits = [min(meanResp(:)), max(meanResp(:))];
        yPad = 0.08 * (yLimits(2) - yLimits(1));
        yLimits = yLimits + [-yPad yPad];
        if yLimits(1) >= yLimits(2)
            yLimits = yLimits + [-1 1];
        end

        for k = 1:nSizes
            subplot(1, nSizes, k)
            hold on

            % shaded stimulus-OFF window (drawn first, so the trace is on top)
            if shadeOff
                offStart = offWindow(1);
                offEnd   = min(offWindow(2), tAxis(end));
                patch([offStart offEnd offEnd offStart], ...
                      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
                      [0.85 0.85 0.85], 'EdgeColor', 'none');
            end

            plot(tAxis, meanResp(:,k), 'k', 'LineWidth', 2);

            xlim([tAxis(1) tAxis(end)])
            ylim(yLimits)
            xlabel(xLabelText)
            if k == 1
                ylabel('Baseline-subtracted response (a.u.)')
            end
            title(sprintf('Spot size %g', spotSizes(k)))

            box off
            set(gca, 'FontSize', 11)
        end

        titleText = sprintf('%s  (%d of %d pixels qualified)', strrep(baseName,'_','\_'), nQualified, nPixels);
        if shadeOff
            titleText = [titleText '   gray = stimulus OFF'];
        end
        sgtitle(titleText)

        saveas(gcf, fullfile(outDir, [baseName '_SpotAnalysis_meanResp.png']));
    end

    fprintf('  Saved: %s_SpotAnalysis.mat\n', baseName);
end

fprintf('\nAll recordings finished.\n');
