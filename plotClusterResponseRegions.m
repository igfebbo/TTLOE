function plotClusterResponseRegions(clusterIdx, alignedSpikeMatrices, rateTables,responseTable, stimuliParams, blockParams,TTL_struct, sampleRate, analysis_folder)
    % Extract metadata for this cluster
    ClusterID = alignedSpikeMatrices(clusterIdx).ClusterID;  
    nBlocks = numel(blockParams);
    maxAligns = max(cellfun(@numel, {blockParams.alignmentTTLnames}));
    % nTiles = 0;
    % for b = 1:nBlocks
    % nTiles = nTiles + numel(blockParams(b).alignmentTTLnames);
    % end
    % nRows = ceil(sqrt(nTiles));
    % nCols = ceil(nTiles / nRows);
 

    %% ----- RASTER FIGURE -----
    figRaster = figure('Name', sprintf('Cluster %d – Rasters', ClusterID), 'Color', 'w');
    %tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(sprintf('Cluster %d – Rasters', ClusterID));
    

    for b = 1:nBlocks
        blockLabel = blockParams(b).label;
        ttlNames   = blockParams(b).alignmentTTLnames;
        nAligns    = numel(ttlNames);
        

        for t = 1:nAligns
            stim = rateTables(clusterIdx).StimulusLabels{1,1}{1,t};
            offset     = blockParams(b).alignmentWindows{t}(1)/sampleRate;
            nexttile((b-1)*maxAligns + t);  % tile index
            hold on;
            %subplot(nBlocks, nAligns, (b-1)*nAligns + t);
            spikeMatrix = alignedSpikeMatrices(clusterIdx).SpikeMatrix{b}{t};
            if ~isempty(spikeMatrix)
                spikeMat_raster(spikeMatrix, sampleRate, 'Offset', offset);
            end
            hold on;
            title(sprintf('%s – %s', blockLabel, ttlNames{t}));
            xlabel('Time (s)');
            ylabel('Trial');
            n_trials = size(spikeMatrix, 1);
            % shade burst times

            [approachOnset, retractOnset]= shadeBurstRegion(responseTable, ClusterID, stim, n_trials, 0.050, offset);

            % Debug: vertical line
            if ~isnan(approachOnset)
                xline(approachOnset, 'r', 'LineWidth', 2);
            end
            if ~isnan(retractOnset)
                xline(retractOnset, 'b', 'LineWidth', 2);
            end

  
            % lgd = legend(yHandles, yLabels, ...
            %      'Orientation', 'horizontal');
            % % Anchor legend to the figure, not a subplot
            %lgd.Units = 'normalized';
            %lgd.Position = [0.3, 0.01, 0.4, 0.05];  % [x y w h] in normalized figure units
            % %set(lgd, 'Box', 'off');
% Save as editable MATLAB figure
figFolder = fullfile(analysis_folder, 'figures');
% Check if it exists; if not, create it
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end
%savefig(figRaster, fullfile(figFolder, sprintf('Cluster%d_Rasters.fig', ClusterID)));
        end
    end

%     %% ----- PSTH FIGURE -----
%     figPSTH = figure('Name', sprintf('Cluster %d – PSTHs', ClusterID), 'Color', 'w');
%     %tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
%     sgtitle(sprintf('Cluster %d – PSTHs', ClusterID));
%     maxFR = 50;
% 
%     for b = 1:nBlocks
%         blockLabel = blockParams(b).label;
%         ttlNames   = blockParams(b).alignmentTTLnames;
%         nAligns    = numel(ttlNames);
% 
%         for t = 1:nAligns
%             offset     = blockParams(b).alignmentWindows{t}(1)/sampleRate;
%             nexttile((b-1)*maxAligns + t);  % tile index
%             %subplot(nBlocks, nAligns, (b-1)*nAligns + t);
%             rateTable = rateTables(clusterIdx).RateTable{b}{t};
%             if ~isempty(rateTable)
%                 plot(rateTable.Time_sec+offset, rateTable.FR_mean_Hz, 'LineWidth', 1.5); %hold on;
%                 %plot(rateTable.Time_sec, rateTable.FR_z, '--', 'LineWidth', 1.2);
%                 %hold off;
%             end
%             title(sprintf('%s – %s', blockLabel, ttlNames{t}));
%             xlabel('Time (s)');
%             ylabel('Firing Rate (Hz)');
%             shadeStimuli(stimuliParams, TTL_struct,blockParams(b).alignmentTTLs{t},blockParams(b).alignmentWindows{t}, maxFR);
% %savefig(figPSTH, fullfile(figFolder, sprintf('Cluster%d_PSTHs.fig', ClusterID)));
%         end
%     end
end
