function [plottingWindow, shadedRegions] = createPlottingWindow(sampleRate, timeWindow_sec, stimuliParams)
% createPlottingWindow creates a time window and shaded region info for multiple stimuli
%
% Inputs:
% - sampleRate: sample rate in Hz (samples/sec)
% - timeWindow_sec: [start, end] seconds around stimulus (e.g. [-0.5, 1.5])
% - stimuliParams: struct array with fields:
%       .duration (sec)
%       .start (sec, optional, default 0)
%       .color (char or RGB, optional, default 'b')
%       .alpha (scalar 0-1, optional, default 0.2)
%       .label (string, optional)
%
% Outputs:
% - plottingWindow: [startSamples, endSamples] (relative to stimulus)
% - shadedRegions: struct array with fields for each stimulus:
%       .startSamples
%       .endSamples
%       .color
%       .alpha
%       .label

% Default start times and colors if missing
nStim = length(stimuliParams);
for i = 1:nStim
    if ~isfield(stimuliParams(i), 'start') || isempty(stimuliParams(i).start)
        stimuliParams(i).start = 0;
    end
    if ~isfield(stimuliParams(i), 'color') || isempty(stimuliParams(i).color)
        stimuliParams(i).color = 'b';
    end
    if ~isfield(stimuliParams(i), 'alpha') || isempty(stimuliParams(i).alpha)
        stimuliParams(i).alpha = 0.2;
    end
    if ~isfield(stimuliParams(i), 'label')
        stimuliParams(i).label = sprintf('Stimulus %d', i);
    end
end

% Convert plotting window to samples
plottingWindow = round(timeWindow_sec * sampleRate);

% Initialize output struct for shaded regions
shadedRegions = struct('startSamples', [], 'endSamples', [], 'color', [], 'alpha', [], 'label', []);

for i = 1:nStim
    shadedRegions(i).startSamples = round(stimuliParams(i).start * sampleRate);
    shadedRegions(i).endSamples = round((stimuliParams(i).start + stimuliParams(i).duration) * sampleRate);
    shadedRegions(i).color = stimuliParams(i).color;
    shadedRegions(i).alpha = stimuliParams(i).alpha;
    shadedRegions(i).label = stimuliParams(i).label;
end

end
