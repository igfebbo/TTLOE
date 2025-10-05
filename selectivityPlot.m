function selectivityPlot(rateTables, blockParams, spikedata, n_clusters)

% Next try
% Cluster selectivity plot with cluster labels
% Assumes: rateTables, blockParams, sampleRate, n_clusters, spikedata available

responseWindow = [0 4]; % seconds post stimulus

% Pre-allocate vectors
rough_vector = nan(n_clusters,1);
smooth_vector = nan(n_clusters,1);

% Loop over blocks and alignments to extract mean firing rate in response window
for b = 1:numel(blockParams)
    ttlNames = blockParams(b).alignmentTTLnames; % e.g. {'smoothTTL_off', 'roughTTL_off'}
    
    for t = 1:numel(ttlNames)
        ttlName = ttlNames{t};
        
        for c = 1:n_clusters
            rateTable = rateTables(c).RateTable{b}{t};
            if isempty(rateTable)
                continue;
            end
            
            % Find indices within the response window
            timeIdx = rateTable.Time_sec >= responseWindow(1) & rateTable.Time_sec <= responseWindow(2);
            
            meanFR = mean(rateTable.FR_raw_Hz(timeIdx));
            
            % Assign mean firing rates to rough or smooth vectors
            if contains(lower(ttlName), 'rough')
                rough_vector(c) = meanFR;
            elseif contains(lower(ttlName), 'smooth')
                smooth_vector(c) = meanFR;
            end
        end
    end
end

% Identify valid clusters (non-NaN for both rough and smooth)
valid_idx = ~isnan(rough_vector) & ~isnan(smooth_vector);

% Get cluster IDs for valid clusters
cluster_ids = arrayfun(@(x) x.ClusterID, spikedata);
valid_cluster_ids = cluster_ids(valid_idx);

% Extract valid firing rates for plotting and labeling
rough_valid = rough_vector(valid_idx);
smooth_valid = smooth_vector(valid_idx);

% Create scatter plot
figure;
scatter(rough_valid, smooth_valid, 40, 'filled','MarkerFaceAlpha', 0.3);
xlabel('Rough firing rate (Hz)');
ylabel('Smooth firing rate (Hz)');
title('Cluster Selectivity');
hold on;

% Set axis limits with some padding
x_min = min(rough_valid) - 0.2;
x_max = max(rough_valid) + 0.2;
y_min = min(smooth_valid) - 0.2;
y_max = max(smooth_valid) + 0.2;
xlim([x_min, x_max]);
ylim([y_min, y_max]);

% Plot diagonal reference line
plot([x_min x_max], [y_min y_max], 'k--');

% Label each point with the cluster ID
for i = 1:length(valid_cluster_ids)
    text(rough_valid(i) + 0.02, smooth_valid(i), num2str(valid_cluster_ids(i)), 'FontSize', 6);
end

drawnow;

end 