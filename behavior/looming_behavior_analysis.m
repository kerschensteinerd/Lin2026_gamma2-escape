%% LOOMING_BEHAVIOR_ANALYSIS  Looming defensive-behavior analysis: heatmaps, reaction time, max speed, escape probability
%
%   PURPOSE
%   Behavior-analysis figures for the looming-stimulus defensive-response
%   assay, comparing WT vs. KO mice across 3 stimulus contrast levels
%   (22%, 64%, 91%). Contains 4 analysis/plotting sections:
%     1. Per-mouse speed heatmaps (one panel per genotype x contrast;
%        WT panels use a white->black colormap, KO panels white->orange)
%     2. Reaction time (first frame speed exceeds a threshold, relative
%        to the looming trigger) vs. contrast, WT vs. KO
%     3. Max speed (within a post-trigger analysis window) vs. contrast, WT vs. KO
%     4. Escape probability vs. contrast, WT vs. KO -- COMPUTED from the
%        data (fraction of mice whose speed ever exceeds the same
%        threshold used in Section 2)
%
%
%   All four sections use a supplied merged MAT file containing HeatmapData.
%   The input schema and analysis settings are described in README.md.


clear; clc; close all;

%% ===================================================================
%  LOAD MERGED DATASETS (once, shared by all sections below)
%  ===================================================================
[dataFile, dataDir] = uigetfile('*.mat', 'Select merged behavior data (HeatmapData)');
if isequal(dataFile,0)
    error('No behavior data file selected. Script stopped.');
end
fig6File = fullfile(dataDir, dataFile);
S = load(fig6File, 'HeatmapData');
HeatmapData_fig6 = S.HeatmapData;

clear S

% Shared helper: pull one genotype x contrast heatmap out of a
% HeatmapData struct array (errors if no match / more than one match).
getHeatmap = @(HD, genotype, contrast) HD( ...
    strcmp({HD.Genotype}, genotype) & [HD.Contrast] == contrast ...
    ).heatmap;

fps = 30;
contrastList = [22 64 91];
genotypeList = {'WT','KO'};

KO_Color = [235 117 46] / 255;
WT_Color = [0 0 0];


%% ===================================================================
%  SECTION 1: PER-MOUSE SPEED HEATMAPS (WT vs KO, 3 contrasts each)
%  ===================================================================

% ----- USER INPUTS -----
upsampleFactor = 5;     % 1 = no smoothing, 5~10 = smoother

maxDisplayValue = 90;   % values >= this become max color

backgroundShift = 10;
% 0 will use the same color as old value = 10

% looming event times (seconds relative to trace start)
loomTimes = [2 2.8 3.6 4.4 5.2];

% display window
xWindow = [1 8];

% ----- COLORMAPS: one per genotype -----
% KO: white -> orange (KO_Color). WT: white -> black (WT_Color).
cmapN = 256;
cmapKO = buildBackgroundShiftedCmap(KO_Color, backgroundShift, maxDisplayValue, cmapN);
cmapWT = buildBackgroundShiftedCmap(WT_Color, backgroundShift, maxDisplayValue, cmapN);

% ----- FIGURE -----
figure('Color','w', ...
       'Position',[100 100 1800 900]);

% 2 rows (genotype) x 3 columns (contrast)
panelIdx = 0;
for g = 1:numel(genotypeList)
    genotype = genotypeList{g};

    if strcmp(genotype,'WT')
        cmapThis = cmapWT;
    else
        cmapThis = cmapKO;
    end

    for c = 1:numel(contrastList)
        contrast = contrastList(c);
        panelIdx = panelIdx + 1;

        match = strcmp({HeatmapData_fig6.Genotype}, genotype) & [HeatmapData_fig6.Contrast] == contrast;
        if ~any(match)
            warning('No data found for %s contrast %d. Skipping panel.', genotype, contrast);
            continue;
        end
        data = HeatmapData_fig6(find(match,1)).heatmap;

        nFrames = size(data,2);
        nMice   = size(data,1);

        % original time axis
        t = (0:nFrames-1) / fps;

        % ----- upsample for smoother display -----
        t_up = linspace(t(1), t(end), nFrames * upsampleFactor);
        data_up = interp1(t, data', t_up, 'pchip')';

        % ----- plot heatmap -----
        ax = subplot(2,3,panelIdx);

        imagesc(t_up, 1:nMice, data_up);

        set(ax, 'YDir', 'normal');

        colormap(ax, cmapThis);
        colorbar;

        caxis([0 maxDisplayValue]);

        % ----- axis settings -----
        xlim(xWindow)

        % remove all text
        xticks([])
        yticks([])

        xlabel('')
        ylabel('')
        title(sprintf('%s, %d%%', genotype, contrast))

        box off

        % ----- looming event lines -----
        hold on
        for j = 1:length(loomTimes)
            xline(loomTimes(j), '--k', 'LineWidth',1.5);
        end
    end
end


%% ===================================================================
%  SECTION 2: REACTION TIME vs CONTRAST (WT vs KO)
%  ===================================================================

% ----- USER INPUTS -----
speedThreshold = 35;
triggerTime = 2;         % zero-reference point for reaction time
analysisWindow = [2 10]; % window searched for a threshold crossing

KO_heatmapList = { getHeatmap(HeatmapData_fig6, 'KO', contrastList(1)), ...
                   getHeatmap(HeatmapData_fig6, 'KO', contrastList(2)), ...
                   getHeatmap(HeatmapData_fig6, 'KO', contrastList(3)) };

WT_heatmapList = { getHeatmap(HeatmapData_fig6, 'WT', contrastList(1)), ...
                   getHeatmap(HeatmapData_fig6, 'WT', contrastList(2)), ...
                   getHeatmap(HeatmapData_fig6, 'WT', contrastList(3)) };

% ----- CALCULATE KO REACTION TIME -----
KO_RT_all = cell(1, length(KO_heatmapList));

for i = 1:length(KO_heatmapList)

    data = KO_heatmapList{i};

    nMice   = size(data,1);
    nFrames = size(data,2);

    t = (0:nFrames-1) / fps;

    % Frames after looming trigger
    validIdx = find(t >= analysisWindow(1) & t <= analysisWindow(2));

    reactionTime = nan(nMice,1);

    for m = 1:nMice

        trace = data(m,:);

        % Find first time speed exceeds threshold
        crossIdx = find(trace(validIdx) > speedThreshold, 1, 'first');

        if ~isempty(crossIdx)

            actualFrame = validIdx(crossIdx);

            % Reaction time relative to looming trigger
            reactionTime(m) = t(actualFrame) - triggerTime;

        end

    end

    KO_RT_all{i} = reactionTime;

end

% ----- KO MEAN + SEM -----
KO_meanRT = nan(1, length(KO_heatmapList));
KO_semRT  = nan(1, length(KO_heatmapList));

for i = 1:length(KO_heatmapList)

    rt = KO_RT_all{i};

    KO_meanRT(i) = mean(rt, 'omitnan');

    n = sum(~isnan(rt));

    KO_semRT(i) = std(rt, 'omitnan') / sqrt(n);

end

% ----- CALCULATE WT REACTION TIME -----
WT_RT_all = cell(1, length(WT_heatmapList));

for i = 1:length(WT_heatmapList)

    data = WT_heatmapList{i};

    nMice   = size(data,1);
    nFrames = size(data,2);

    t = (0:nFrames-1) / fps;

    % Frames after looming trigger
    validIdx = find(t >= analysisWindow(1) & t <= analysisWindow(2));

    reactionTime = nan(nMice,1);

    for m = 1:nMice

        trace = data(m,:);

        % Find first time speed exceeds threshold
        crossIdx = find(trace(validIdx) > speedThreshold, 1, 'first');

        if ~isempty(crossIdx)

            actualFrame = validIdx(crossIdx);

            % Reaction time relative to looming trigger
            reactionTime(m) = t(actualFrame) - triggerTime;

        end

    end

    WT_RT_all{i} = reactionTime;

end

% ----- WT MEAN + SEM -----
WT_meanRT = nan(1, length(WT_heatmapList));
WT_semRT  = nan(1, length(WT_heatmapList));

for i = 1:length(WT_heatmapList)

    rt = WT_RT_all{i};

    WT_meanRT(i) = mean(rt, 'omitnan');

    n = sum(~isnan(rt));

    WT_semRT(i) = std(rt, 'omitnan') / sqrt(n);

end

% ----- PLOT -----
figure('Color','w');

hold on

errorbar(contrastList, WT_meanRT, WT_semRT, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'CapSize',0, ...
    'Color',WT_Color, ...
    'MarkerEdgeColor',WT_Color, ...
    'MarkerFaceColor',WT_Color);

errorbar(contrastList, KO_meanRT, KO_semRT, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'CapSize',0, ...
    'Color',KO_Color, ...
    'MarkerEdgeColor',KO_Color, ...
    'MarkerFaceColor',KO_Color);

xlabel('Contrast (%)');
ylabel('Reaction Time (s)');

ylim([0 5.4])
xlim([20 100])
xticks(contrastList)

legend({'WT','KO'}, ...
    'Location','best', ...
    'Box','off');

box off
set(gca,'FontSize',14)


%% ===================================================================
%  SECTION 3: MAX SPEED vs CONTRAST (WT vs KO)
%  ===================================================================

% ----- USER INPUTS -----
analysisWindow = [2 10];

KO_heatmapList = { getHeatmap(HeatmapData_fig6, 'KO', contrastList(1)), ...
                   getHeatmap(HeatmapData_fig6, 'KO', contrastList(2)), ...
                   getHeatmap(HeatmapData_fig6, 'KO', contrastList(3)) };

WT_heatmapList = { getHeatmap(HeatmapData_fig6, 'WT', contrastList(1)), ...
                   getHeatmap(HeatmapData_fig6, 'WT', contrastList(2)), ...
                   getHeatmap(HeatmapData_fig6, 'WT', contrastList(3)) };

% ----- CALCULATE MAX SPEED -----
KO_all = cell(1,3);
WT_all = cell(1,3);

KO_meanMax = nan(1,3);
KO_semMax  = nan(1,3);

WT_meanMax = nan(1,3);
WT_semMax  = nan(1,3);

for i = 1:3

    % ----- KO -----
    data = KO_heatmapList{i};

    nFrames = size(data,2);

    t = (0:nFrames-1) / fps;

    validIdx = find(t >= analysisWindow(1) & ...
                    t <= analysisWindow(2));

    maxVals = nan(size(data,1),1);

    for m = 1:size(data,1)

        trace = data(m,:);

        maxVals(m) = max(trace(validIdx));

    end

    KO_all{i} = maxVals;

    KO_meanMax(i) = mean(maxVals, 'omitnan');

    KO_semMax(i) = std(maxVals, 'omitnan') / ...
        sqrt(sum(~isnan(maxVals)));

    % ----- WT -----
    data = WT_heatmapList{i};

    nFrames = size(data,2);

    t = (0:nFrames-1) / fps;

    validIdx = find(t >= analysisWindow(1) & ...
                    t <= analysisWindow(2));

    maxVals = nan(size(data,1),1);

    for m = 1:size(data,1)

        trace = data(m,:);

        maxVals(m) = max(trace(validIdx));

    end

    WT_all{i} = maxVals;

    WT_meanMax(i) = mean(maxVals, 'omitnan');

    WT_semMax(i) = std(maxVals, 'omitnan') / ...
        sqrt(sum(~isnan(maxVals)));

end

% ----- PLOT -----
figure('Color','w');

hold on

% WT
errorbar(contrastList, WT_meanMax, WT_semMax, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'CapSize',0, ...
    'Color',WT_Color, ...
    'MarkerEdgeColor',WT_Color, ...
    'MarkerFaceColor',WT_Color);

% KO
errorbar(contrastList, KO_meanMax, KO_semMax, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'CapSize',0, ...
    'Color',KO_Color, ...
    'MarkerEdgeColor',KO_Color, ...
    'MarkerFaceColor',KO_Color);

xlabel('Contrast (%)');
ylabel('Max Speed');

xlim([20 100])
ylim([0 110])

xticks(contrastList)

box off
set(gca,'FontSize',14)

legend({'WT','KO'}, ...
    'Box','off', ...
    'Location','best');


%% ===================================================================
%  SECTION 4: ESCAPE PROBABILITY vs CONTRAST (WT vs KO)
%  Computed from data: escape probability = (# mice whose speed ever
%  exceeds speedThreshold after the trigger) / (total mice), using the
%  SAME threshold/trigger definition as Section 2. A mouse "escapes"
%  here exactly when Section 2 would have given it a non-NaN reaction
%  time. Uses the same HeatmapData_fig6 dataset as every other section.
%  ===================================================================

% ----- USER INPUTS (matches Section 2's definition) -----
speedThreshold = 35;
triggerTime = 2;
analysisWindow = [2 10];

KO_heatmapList = { getHeatmap(HeatmapData_fig6, 'KO', contrastList(1)), ...
                   getHeatmap(HeatmapData_fig6, 'KO', contrastList(2)), ...
                   getHeatmap(HeatmapData_fig6, 'KO', contrastList(3)) };

WT_heatmapList = { getHeatmap(HeatmapData_fig6, 'WT', contrastList(1)), ...
                   getHeatmap(HeatmapData_fig6, 'WT', contrastList(2)), ...
                   getHeatmap(HeatmapData_fig6, 'WT', contrastList(3)) };

KO_escapeProb = nan(1,3);
WT_escapeProb = nan(1,3);
KO_escapeCount = nan(1,3); KO_nMice = nan(1,3);
WT_escapeCount = nan(1,3); WT_nMice = nan(1,3);

for i = 1:3

    % ----- KO -----
    data = KO_heatmapList{i};
    nMice   = size(data,1);
    nFrames = size(data,2);
    t = (0:nFrames-1) / fps;
    validIdx = find(t >= analysisWindow(1) & t <= analysisWindow(2));

    escaped = false(nMice,1);
    for m = 1:nMice
        trace = data(m,:);
        escaped(m) = any(trace(validIdx) > speedThreshold);
    end

    KO_escapeCount(i) = sum(escaped);
    KO_nMice(i) = nMice;
    KO_escapeProb(i) = KO_escapeCount(i) / KO_nMice(i);

    % ----- WT -----
    data = WT_heatmapList{i};
    nMice   = size(data,1);
    nFrames = size(data,2);
    t = (0:nFrames-1) / fps;
    validIdx = find(t >= analysisWindow(1) & t <= analysisWindow(2));

    escaped = false(nMice,1);
    for m = 1:nMice
        trace = data(m,:);
        escaped(m) = any(trace(validIdx) > speedThreshold);
    end

    WT_escapeCount(i) = sum(escaped);
    WT_nMice(i) = nMice;
    WT_escapeProb(i) = WT_escapeCount(i) / WT_nMice(i);

end

% Print a quick summary, e.g. "WT 22%: 8/23 = 0.348"
fprintf('\nEscape probability summary:\n');
for i = 1:3
    fprintf('  WT %d%%: %d/%d = %.3f\n', contrastList(i), WT_escapeCount(i), WT_nMice(i), WT_escapeProb(i));
    fprintf('  KO %d%%: %d/%d = %.3f\n', contrastList(i), KO_escapeCount(i), KO_nMice(i), KO_escapeProb(i));
end

% ----- PLOT -----
figure('Color','w');

hold on

% WT
plot(contrastList, WT_escapeProb, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'Color',WT_Color, ...
    'MarkerEdgeColor',WT_Color, ...
    'MarkerFaceColor',WT_Color);

% KO
plot(contrastList, KO_escapeProb, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'Color',KO_Color, ...
    'MarkerEdgeColor',KO_Color, ...
    'MarkerFaceColor',KO_Color);

% ----- AXIS -----
xlabel('Contrast (%)');
ylabel('Escape Probability');

xlim([20 100]);
ylim([0 1.05]);

xticks(contrastList);

box off
set(gca,'FontSize',14)

legend({'WT','KO'}, ...
    'Location','best', ...
    'Box','off');


%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================
function cmap = buildBackgroundShiftedCmap(maxColor, backgroundShift, maxDisplayValue, cmapN)
% Builds an Nx3 colormap that interpolates from a background-shifted
% near-white color up to maxColor, used by Section 1's heatmaps.
% backgroundShift/maxDisplayValue sets how far from pure white the
% color scale starts (0 = pure white start).
    startFrac = backgroundShift / maxDisplayValue;
    baseColor = (1-startFrac)*[1 1 1] + startFrac*maxColor;

    cmap = [linspace(baseColor(1), maxColor(1), cmapN)', ...
            linspace(baseColor(2), maxColor(2), cmapN)', ...
            linspace(baseColor(3), maxColor(3), cmapN)'];
end
