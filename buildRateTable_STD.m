function rate_table = buildRateTable(spikeMatrix, sampleRate, alignmentWindow, binWidth)
    nTrials = size(spikeMatrix,1);
    nSamples = size(spikeMatrix,2);

    % Convert binWidth to samples
    binWidth_samples = round(binWidth * sampleRate);

    % Define bin edges in samples
    binEdges = 1:binWidth_samples:nSamples+1;  % +1 so last bin includes all remaining samples

    % Bin spikes
    spikeCounts = zeros(nTrials, length(binEdges)-1);
    for i = 1:length(binEdges)-1
        spikeCounts(:,i) = sum(spikeMatrix(:, binEdges(i):binEdges(i+1)-1), 2);
    end

    % Compute firing rate (Hz) for each trial and time bin
    firingRate_trials = spikeCounts ./ binWidth;

    % Mean and SD across trials
    firingRate_mean = mean(firingRate_trials, 1);
    firingRate_std = std(firingRate_trials, 0, 1);

    % Time vector (center of bins in seconds)
    timeVect = ((binEdges(1:end-1) + binEdges(2:end)-1)/2) / sampleRate;

    % Baseline z-scoring
    baselineMask = timeVect < 0;
    baselineFR = firingRate_mean(baselineMask);
    if isempty(baselineFR)
        baselineMean = mean(firingRate_mean);
        baselineStd = std(firingRate_mean);
    else
        baselineMean = mean(baselineFR);
        baselineStd = std(baselineFR);
    end
    if baselineStd == 0
        baselineStd = 1e-6;
    end
    firingRate_z = (firingRate_mean - baselineMean) ./ baselineStd;

    % Build table including mean, std, baseline mean/std
    rate_table = table( ...
        timeVect(:), ...
        firingRate_mean(:), ...
        firingRate_std(:), ...
        firingRate_z(:), ...
        repmat(baselineMean, numel(timeVect), 1), ...
        repmat(baselineStd, numel(timeVect), 1), ...
        'VariableNames', {'Time_sec', 'FR_mean_Hz', 'FR_std_Hz', 'FR_z', 'Baseline_mean', 'Baseline_std'});
end