%% DEINTERLEAVE_CHANNELS Split interleaved TIFFs into per-channel MATLAB stacks
% Select a folder of multi-page .tif/.tiff files. Edit nChannels,
% channelsToExtract and channelNames to match the acquisition order.
% Each output <base>_<channel>.mat contains a Height x Width x Frames
% array named for the channel, in the source pixel type (-v7.3 format).
% The default channel order is Gcamp, Tdtomato, TransmittedLight, Digit.
% Incomplete trailing channel cycles are dropped. All TIFFs must share
% their channel order and geometry. Movies are loaded fully into RAM.
% Next: register_channels and extract_stimulus_timing.
% See README.md for configuration, file flow and execution order.

clear; clc;

%% ===== USER INPUT =====
nChannels = 4;                       % total number of interleaved channels in each TIFF
channelsToExtract = [1 2 3 4];       % 1-based indices of channels to keep (default: extract all channels)

% Name for each possible channel index (index i = name of channel i).
% Used both as the saved variable name and as the "_<name>.mat" suffix.
channelNames = {'Gcamp','Tdtomato','TransmittedLight','Digit'};

% Name of the subfolder (created inside inputDir) where extracted
% per-channel .mat files will be saved.
OUTPUT_SUBFOLDER_NAME = sprintf('extract_%dchannel', numel(channelsToExtract));

%% ===== SELECT INPUT FOLDER =====
% A dialog box opens for you to choose the folder containing your raw,
% multi-channel TIFF files. To skip the dialog (e.g., for batch/automated
% runs), comment out the two lines below and instead hardcode a path:
%   inputDir = 'C:\path\to\your\data';
inputDir = uigetdir(pwd, 'Select folder containing raw multi-channel TIFFs');
if isequal(inputDir,0)
    error('No input folder selected. Script stopped.');
end

%% ===== MAKE OUTPUT FOLDER =====
outputDir = fullfile(inputDir, OUTPUT_SUBFOLDER_NAME);
if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

%% ===== FIND TIFF FILES =====
files = [dir(fullfile(inputDir,'*.tif')); dir(fullfile(inputDir,'*.tiff'))];

fprintf('Found %d files\n', length(files));

%% ===== PROCESS FILES =====
for f = 1:length(files)

    fileName = files(f).name;
    filePath = fullfile(inputDir,fileName);
    [~,baseName,~] = fileparts(fileName);

    fprintf('\nProcessing %s\n', fileName);

    % Read TIFF page metadata (one entry per page/frame in the file)
    info = imfinfo(filePath);
    nImages = numel(info);
    nFrames = floor(nImages / nChannels);   % frames per channel; drops any incomplete trailing cycle

    H = info(1).Height;
    W = info(1).Width;

    % Detect pixel datatype (e.g., uint16) from the first page so the
    % preallocated arrays below match the source data's class exactly.
    tmp = imread(filePath,1,'Info',info);
    imgClass = class(tmp);

    % Preallocate one H x W x nFrames array per requested channel
    % (only for the channels in channelsToExtract, to save memory).
    data = cell(1,length(channelsToExtract));

    for c = 1:length(channelsToExtract)
        data{c} = zeros(H,W,nFrames,imgClass);
    end

    % Walk through the file frame-by-frame; for each frame, pull out
    % only the pages belonging to the requested channels.
    % Page index formula: interleaved order means the page for channel
    % "ch" at timepoint "t" sits at position (t-1)*nChannels + ch.
    for t = 1:nFrames

        for c = 1:length(channelsToExtract)

            ch = channelsToExtract(c);

            pageIdx = (t-1)*nChannels + ch;

            data{c}(:,:,t) = imread(filePath,pageIdx,'Info',info);

        end

    end

    % Save each extracted channel to its own .mat file, named after the
    % channel and tagged to the source filename.
    for c = 1:length(channelsToExtract)

        ch = channelsToExtract(c);
        varName = channelNames{ch};

        eval([varName ' = data{c};']);   % dynamically name the variable, e.g. Gcamp = data{c};

        save(fullfile(outputDir,[baseName '_' varName '.mat']), ...
             varName,'-v7.3');            % '-v7.3' supports arrays >2GB

        clear(varName)

    end

    clear data

end

fprintf('\nFinished.\n');
