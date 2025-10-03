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

    % Compute firing rate (Hz)
    firingRate_raw = sum(spikeCounts,1) ./ nTrials ./ binWidth;

    % Time vector (center of bins in seconds)
    timeVect = ((binEdges(1:end-1) + binEdges(2:end)-1)/2) / sampleRate;

    % Baseline z-scoring
    baselineMask = timeVect < 0;
    baselineFR = firingRate_raw(baselineMask);
    if isempty(baselineFR)
        baselineMean = mean(firingRate_raw);
        baselineStd = std(firingRate_raw);
    else
        baselineMean = mean(baselineFR);
        baselineStd = std(baselineFR);
    end
    if baselineStd == 0
        baselineStd = 1e-6;
    end
    firingRate_z = (firingRate_raw - baselineMean) ./ baselineStd;

    % Build table
    rate_table = table(timeVect(:), firingRate_raw(:), firingRate_z(:), ...
        'VariableNames', {'Time_sec', 'FR_raw_Hz', 'FR_z'});
end
