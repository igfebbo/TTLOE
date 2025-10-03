function[alignedSpikeMatrices, rateTables] = clustersByBlock(spikedata, blockParams,sampleRate, binWidth)
% Analyze neural clusters across blocks based on a third TTL signal
%
% Inputs:
% - spikedata: struct with SpikeTimes (in samples) and ClusterID
% - stimTTL_on: stimulus TTL times (in samples)
% - blockTTL_on: block TTL times (in samples)
% - sampleRate: Hz
% - alignmentWindow_sec: [start, end] around stim in seconds
% - binWidth: PSTH bin width in seconds
% - output_folder: folder to save figures


    %% Initialization
    nClusters = length(spikedata);
    nBlocks = numel(blockParams);
    alignedSpikeMatrices(nClusters) = struct('ClusterID', [],'SpikeMatrix', [], 'BlockLabels', [], 'StimulusLabels', [], ...
                                         'AlignmentWindows', [] );
    rateTables(nClusters) = struct('ClusterID', [], 'BlockLabels', [], 'StimulusLabels', [], ...
                               'AlignmentWindows', [], 'RateTable', []);
    %peakStats = table();


%% Loop over clusters

for i = 1:nClusters
    clusterID = spikedata(i).ClusterID;
    spikeTimes_samples = double(spikedata(i).SpikeTimes);
    
    % Initialize struct for this cluster
    alignedSpikeMatrices(i).ClusterID = clusterID;
    alignedSpikeMatrices(i).SpikeMatrix = cell(nBlocks, 1);
    alignedSpikeMatrices(i).BlockLabels = cell(nBlocks, 1);
    alignedSpikeMatrices(i).StimulusLabels = cell(nBlocks, 1);
    alignedSpikeMatrices(i).AlignmentWindows = cell(nBlocks, 1);

    %% loop through blocks

    for b = 1:nBlocks
        blockLabel = blockParams(b).label;
        ttlList    = blockParams(b).alignmentTTLs;
        nAligns    = numel(ttlList);

        % Pre-allocate block level
        alignedSpikeMatrices(i).SpikeMatrix{b} = cell(1, nAligns);
        alignedSpikeMatrices(i).StimulusLabels{b} = cell(1, nAligns);
        alignedSpikeMatrices(i).AlignmentWindows{b} = cell(1, nAligns);

        
    %% loop over stimuli per block to generate spike matrix
        for t = 1:nAligns
            alignmentWindow = blockParams(b).alignmentWindows{t};
            alignmentTTLs   = blockParams(b).alignmentTTLs{t};
            stimLabel       = blockParams(b).alignmentTTLnames{t};
 
            % Align spikes
            [~,~,~,~,spikeMatrix] = alignEvents(alignmentTTLs, spikeTimes_samples, alignmentWindow);

            
            % Store results
            alignedSpikeMatrices(i).SpikeMatrix{b}{t}      = spikeMatrix;
            alignedSpikeMatrices(i).BlockLabels{b}{t}      = blockLabel;
            alignedSpikeMatrices(i).StimulusLabels{b}{t}   = stimLabel;
            alignedSpikeMatrices(i).AlignmentWindows{b}{t} = alignmentWindow;   
                
            if isempty(spikeMatrix), continue; end


            rateTable = buildRateTable(spikeMatrix, sampleRate, alignmentWindow, binWidth);

            rateTables(i).ClusterID = clusterID;
            rateTables(i).BlockLabels{b}{t} = blockLabel;
            rateTables(i).StimulusLabels{b}{t} = stimLabel;
            rateTables(i).AlignmentWindows{b}{t} = alignmentWindow;
            rateTables(i).RateTable{b}{t} = rateTable;

        end
    end
end
%% commented fully out because doing figures in another function, keeping for reference though
% for i = 27:29 %1:nClusters
%     clusterID = spikedata(i).ClusterID;
%     for b = 1:nBlocks
%     [rate_table] = buildRateTable(alignedSpikeMatrix{i,b}, sampleRate, alignmentWindow_sec, binWidth);
%     allrateTables{i,b} = rate_table;
%     end
% %end
% 
% for i = 27:29 %1:nClusters
%     clusterID = spikedata(i).ClusterID;
%     maxFR = 30;
% 
%     % Create a new figure for this cluster
%     fig = figure('Visible', 'on', 'Name', sprintf('Cluster %d', clusterID));
% 
%     legendHandlesAll = []; % collect unique handles across blocks
%     legendLabelsAll = {};
% 
%     %% --- Loop over blocks ---
%     for b = 1:nBlocks
%         % --- Raster ---
%         subplot(nBlocks, 2, (b-1)*2 + 1);
%         spikeMat_raster(alignedSpikeMatrix{i,b}, sampleRate, [], ...
%                         'Offset', alignmentWindow_sec(1));
%         title(sprintf('%s Raster', blockLabels{b}));
%         xlabel('Time (s)'); 
%         ylabel('Trial'); 
%         hold on;
%         shadeStimuli(stimuliParams, TTL_struct, TTLnames, block(b), n_trials);
% 
%         % --- PSTH ---
%         subplot(nBlocks, 2, (b-1)*2 + 2);
%         timeVect = allrateTables{i,b}.Time_sec + alignmentWindow_sec(1);
%         plot(timeVect, allrateTables{i,b}.FR_raw_Hz, 'k', 'LineWidth', 2);
%         title(sprintf('%s PSTH', blockLabels{b}));
%         xlabel('Time (s)'); ylabel('FR');
%         ylim([0 maxFR]); hold on;
% 
%         % --- Shade stimuli for this block ---
%         [yHandles, yLabels] = shadeStimuli(stimuliParams, TTL_struct, TTLnames, block(b), maxFR);
% 
%         % collect unique handles/labels for the *global* legend
%         for k = 1:numel(yLabels)
%             if ~ismember(yLabels{k}, legendLabelsAll)
%                 legendHandlesAll(end+1) = yHandles(k); %#ok<AGROW>
%                 legendLabelsAll{end+1} = yLabels{k}; %#ok<AGROW>
%             end
%         end
%     end
% 
%     % --- Add one legend for the whole figure ---
%     if ~isempty(legendHandlesAll)
%     lgd = legend(legendHandlesAll, legendLabelsAll, ...
%                  'Orientation', 'horizontal');
%     % Anchor legend to the figure, not a subplot
%     lgd.Units = 'normalized';
%     lgd.Position = [0.3, 0.01, 0.4, 0.05];  % [x y w h] in normalized figure units
%     %set(lgd, 'Box', 'off');
% end
% 
%     % --- Add a title for the whole figure ---
%     sgtitle(sprintf('Cluster %d', clusterID), 'FontSize', 14, 'FontWeight', 'bold');
% end
% 
% 
% %         %% --- Compute statistics: peak response vs baseline ---
% %         baselineWindow = [-0.5 0];  
% %         responseWindow = [0 0.5];  
% %         baselineIdx = timeVect >= baselineWindow(1) & timeVect <= baselineWindow(2);
% %         responseIdx = timeVect >= responseWindow(1) & timeVect <= responseWindow(2);
% % 
% %         % FR per trial in Hz
% %         FR_per_bin = spikeMatrix/binWidth;
% %         baselineFR = mean(FR_per_bin(:,baselineIdx),2);
% %         responseFR = mean(FR_per_bin(:,responseIdx),2);
% % 
% %         % Wilcoxon signed-rank test
% %         [p,h] = signrank(responseFR, baselineFR);
% % 
% %         % Peak FR in response window
% %         meanResponseFR = mean(FR_per_bin(:,responseIdx),1);
% %         [peakFR, peakIdx] = max(meanResponseFR);
% %         peakTime = timeVect(responseIdx);
% %         peakTime = peakTime(peakIdx);
% % 
% %         % Add significance marker on PSTH
% %         if h==1
% %             hold on
% %             plot(peakTime, peakFR+5, 'r*', 'MarkerSize',8)
% %         end
% % 
% %         % Append stats to table
% %         peakStats = [peakStats; table(clusterID, b, peakFR, peakTime, mean(baselineFR), p, h, ...
% %             'VariableNames', {'ClusterID','Block','PeakFR','PeakTime','BaselineFR','pValue','Significant'})];
% % 
% %       Optionally save figure
% %saveas(fig, fullfile(output_folder, sprintf('cluster%d_block%d.fig', clusterID,b)));
% %close(gcf)
%     end
% %end
% % 
% % % Save peak stats table
% % %writetable(peakStats, fullfile(output_folder,'PeakFR_Stats.csv'));
% 
% 
% %end
