% This is the function version of Cluster Localization on the probe
function clusterLocalizationOnProbefunction(mouseID, day, ksortPath, analysis_folder)
% clusterLocalizationOnProbe
% Visualizes the location of best channels for each cluster on the probe
%
% Inputs:
%   - mouseID: string (e.g., 'kl16')
%   - day: integer (e.g., 1)
%   - ksortPath: string, path to Kilosort output
%   - analysis_folder: string, folder to save the plots

    % Load data
    templates = readNPY(fullfile(ksortPath, 'templates.npy'));  
    jsonText = fileread('/mnt/multiverse/homes/kathi/Channel Map/H9-ASSY-77.json');
    probeMap = jsondecode(jsonText);

    channelPositionKS = readNPY(fullfile(ksortPath, 'channel_positions.npy'));
    clusterInfo = readtable(fullfile(ksortPath, 'cluster_info.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    
    % Get geometry info
    chanMapZeroBased = probeMap.chanMap;
    xc = probeMap.xc;
    yc = probeMap.yc;
    chanMap = chanMapZeroBased + 1; % convert to 1-based indexing

    % Select good clusters
    goodClustersTable = clusterInfo(strcmp(clusterInfo.group, 'good'), :);
    BestChannels = goodClustersTable.ch;
    depths = goodClustersTable.depth;

    GoodClusterInfo = table(goodClustersTable.cluster_id, ...
                            goodClustersTable.ch, ...
                            goodClustersTable.depth, ...
                            'VariableNames', {'ClusterID', 'BestChannel', 'Depth'});

    % Sort by depth
    GoodClusterInfo = sortrows(GoodClusterInfo, "Depth");

    % Convert to 1-based
    BestChannels_1 = BestChannels + 1;

    %% ---- HEATMAP PLOT ----
    counts = zeros(length(chanMap), 1);
    for i = 1:length(chanMap)
        chNum = chanMap(i);
        counts(i) = sum(BestChannels_1 == chNum);
    end

    figure; hold on;
    set(gcf, 'Position', [100, 100, 400, 800]);
    scatter(xc, yc, 50, counts, 'filled');

    for i = 1:length(chanMap)
        text(xc(i), yc(i), num2str(chanMap(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'Color', 'k');
    end

    xlabel('X (µm)');
    ylabel('Y (µm)');
    title(['Best Channels ', mouseID, ', Day ', num2str(day)]);
    xlim([-200 200]);
    ylim([0 1500]);
    axis equal;
    grid on;
    colorbar;
    caxis([0 max(counts)]);

    % Adjust layout to avoid cutting off elements
    set(gca, 'Position', [0.25, 0.1, 0.55, 0.85]);

    saveas(gcf, fullfile(analysis_folder, 'heatmap_best_channels.png'));
    close(gcf);

    %% ---- PLOT MAIN CHANNEL FOR EACH CLUSTER ----
    for cIdx = 1:height(GoodClusterInfo)
        figure; hold on;
        set(gcf, 'Position', [100, 100, 400, 800]);

        % Grey background dots
        scatter(xc, yc, 70, [0.6 0.6 0.6], 'filled');

        % Highlight main channel in red
        mainCh = GoodClusterInfo.BestChannel(cIdx) + 1;
        chIndex = find(chanMap == mainCh);
        scatter(xc(chIndex), yc(chIndex), 100, 'r', 'filled');

        % Add channel labels
        for i = 1:length(chanMap)
            text(xc(i), yc(i), num2str(chanMap(i)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'Color', 'k');
        end

        % Title with cluster ID and depth
        depth = GoodClusterInfo.Depth(cIdx);
        title(sprintf('Cluster %d - Main Channel %d - Depth %.1f µm', ...
            GoodClusterInfo.ClusterID(cIdx), mainCh, depth));

        xlabel('X (µm)');
        ylabel('Y (µm)');
        xlim([-200 200]);
        ylim([0 1500]);
        axis equal; grid on;

        set(gca, 'Position', [0.25, 0.1, 0.55, 0.85]);

        % Save figure
        set(gcf, 'PaperPositionMode', 'auto');
        filename = sprintf('cluster_%d_main_channel.png', GoodClusterInfo.ClusterID(cIdx));
        saveas(gcf, fullfile(analysis_folder, filename));
        close(gcf);
    end
end