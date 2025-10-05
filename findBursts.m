function [onsetTime, peakTime, maxFR] = findBursts(timeSec, frHz, window, thresh)
    onsetTime = NaN; peakTime = NaN; maxFR = NaN;

    % Mask within window
    mask = timeSec >= window(1) & timeSec <= window(2);

    % Identify burst bins
    isBurst = frHz > thresh & mask;

    if any(isBurst)
        % Find contiguous burst segments
        cc = bwconncomp(isBurst);
        % Pick the burst with the highest peak FR
        bestIdx = NaN; bestPeak = -Inf;
        for i = 1:cc.NumObjects
            idx = cc.PixelIdxList{i};
            [peak, pIdx] = max(frHz(idx));
            if peak > bestPeak
                bestPeak = peak;
                bestIdx = idx(pIdx);
                onsetIdx = idx(1);
                offsetIdx = idx(end);
            end
        end
        % Return burst metrics
        onsetTime = timeSec(onsetIdx);
        peakTime  = timeSec(bestIdx);
        maxFR     = bestPeak;
    end
end