function thresholdGUI(data)
%THRESHOLDGUI  Interactive slider tool to pick a threshold and binarize an image
%
%   thresholdGUI(data)
%
%   PURPOSE
%   Shared helper function used during manual dendrite/ROI segmentation
%   (see segment_tdtomato_dendrites.m). Displays a 2-D image next to a
%   live thresholded (binary) preview, lets the user drag a slider to
%   adjust the threshold while watching the preview update, and saves
%   the chosen value to the BASE workspace when "Save Threshold" is
%   clicked.
%
%   INPUT
%   data - any 2-D numeric image (e.g., a mean or std projection image)
%          to threshold.
%
%   OUTPUT
%   None returned directly. When the user clicks "Save Threshold", the
%   chosen scalar value is written to a variable named "threshold" in
%   the BASE workspace via assignin('base', ...) -- i.e., the variable
%   becomes available in the script that called thresholdGUI once this
%   function returns (this only works when the caller is itself a
%   script running at the base workspace level, not another function).
%   Execution of the calling script is paused (uiwait) until this GUI's
%   figure is closed, either via "Save Threshold" or by closing the
%   window directly (in which case "threshold" is never created --
%   callers should check with `exist('threshold','var')` afterward).
%
%   USAGE
%       thresholdGUI(stdImg);
%       if exist('threshold','var')
%           myThreshold = threshold;
%       else
%           disp('Threshold was not set.');
%       end
%
%   DEPENDENCIES
%   Base MATLAB only (uicontrol, imagesc).
%
%   Part of: Vg3-AC-Ca2imaging-pipeline (helper for segmentation)

    % ----- CREATE FIGURE -----
    fig = figure('Position', [100, 100, 1600, 300], 'Name', 'Thresholding GUI');

    % Left panel: original (unthresholded) data
    subplot(1, 2, 1);
    imagesc(data);
    axis off;
    title('Original Data');

    % Slider for interactively adjusting the threshold, initialized to
    % the image's mean value.
    uicontrol('Style', 'slider', 'Min', min(data(:)), 'Max', max(data(:)), ...
        'Value', mean(data(:)), 'Position', [600, 10, 600, 20], ...
        'Callback', @thresholdSliderCallback);

    % Right panel: live preview of the binary mask at the current threshold
    maskedData = applyThreshold(data, mean(data(:)));
    subplot(1, 2, 2);
    imagesc(maskedData);
    axis off;
    title('Masked Data');

    % Button to commit the current threshold value and close the GUI
    uicontrol('Style', 'pushbutton', 'String', 'Save Threshold', ...
        'Position', [1210, 5, 200, 30], 'Callback', @saveThreshold);

    thresholdValue = mean(data(:));  % current threshold, updated live by the slider

    % ----- CALLBACKS -----
    function thresholdSliderCallback(hObject, ~)
        % Update the preview mask whenever the slider moves.
        thresholdValue = get(hObject, 'Value');
        maskedData = applyThreshold(data, thresholdValue);
        subplot(1, 2, 2);
        imagesc(maskedData);
        axis off;
        title(['Masked Data (Threshold = ' num2str(thresholdValue) ')']);
    end

    function maskedData = applyThreshold(data, threshold)
        % Binarize: pixels >= threshold become 1 (true), else 0 (false).
        maskedData = data >= threshold;
    end

    function saveThreshold(~, ~)
        % Commit the current slider value to the BASE workspace as
        % "threshold" so the calling script can read it, then close.
        assignin('base', 'threshold', thresholdValue);
        disp(['Threshold saved: threshold = ' num2str(thresholdValue)]);
        close(fig);
    end

    % Block the caller until this figure is closed (Save Threshold, or
    % the window's close button).
    uiwait(fig);

end
