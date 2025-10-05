function heatmap_new(rateTables, blockParams, sampleRate, blockIdx)

if nargin < 4
        blockIdx = 1; % default to block 1
end

nClusters = numel(rateTables); % total clusters
blockIdx = 1; % assuming your data is in block 1
alignNames = blockParams(blockIdx).alignmentTTLnames; % e.g. {'smoothTTL_off', 'roughTTL_off'}
nTimePoints = length(rateTables(1).RateTable{blockIdx}{1}.Time_sec);

% Pre-allocate heatmaps with NaNs
smoothZ = NaN(nClusters, nTimePoints);
roughZ  = NaN(nClusters, nTimePoints);

for c = 1:nClusters
    for t = 1:numel(alignNames)
        rateTable = rateTables(c).RateTable{blockIdx}{t};
        if ~isempty(rateTable)
            % Get TTL name for this alignment
            ttlName = alignNames{t};
            
            % Assign z-scored firing rate based on TTL name
            if contains(lower(ttlName), 'smooth')
                smoothZ(c, :) = rateTable.FR_z(:)';
            elseif contains(lower(ttlName), 'rough')
                roughZ(c, :) = rateTable.FR_z(:)';
            end
        end
    end
end

% Time vector for plotting (add offset if needed, here assuming no offset)
%timeVec = rateTables(1).RateTable{blockIdx}{1}.Time_sec;
offset = blockParams(blockIdx).alignmentWindows{1}(1) / sampleRate;
timeVec = rateTables(1).RateTable{blockIdx}{1}.Time_sec + offset;
% Apply alignment offset

% Plot heatmaps
figure;
subplot(2,1,1);
imagesc(timeVec, 1:nClusters, smoothZ);
axis xy;
xlim([-1, 6]);
colorbar;
title('Smooth Texture Response (Z-scored FR)');
xlabel('Time (s)');
ylabel('Cluster');
%add vertical lines for presentation and withdrawal of texture stimuli
%xline(0, 'w-', 'LineWidth', 1.5);   % TTL onset

[~, zeroIdx] = min(abs(timeVec - 0));
zeroTime = timeVec(zeroIdx);
xline(0, 'w-', 'LineWidth', 1.5);
xline(3, 'w--', 'LineWidth', 1.2);   % texture withdrawal (adjust as needed)

subplot(2,1,2);
imagesc(timeVec, 1:nClusters, roughZ);
axis xy;
xlim([-1, 6]);
colorbar;
title('Rough Texture Response (Z-scored FR)');
xlabel('Time (s)');
ylabel('Cluster');
% lines for texture
%xline(0, 'w-', 'LineWidth', 1.5);    % TTL onset
[~, zeroIdx] = min(abs(timeVec - 0));
zeroTime = timeVec(zeroIdx);
xline(0, 'w-', 'LineWidth', 1.5);
xline(3, 'w--', 'LineWidth', 1.2);   % texture withdrawal
end 