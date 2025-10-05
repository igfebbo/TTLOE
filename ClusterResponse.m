function responseTable = ClusterResponse(rateTables, analysis_folder, mouseID)
% ClusterOnsetResponse detects bursts of activity (≥ threshold SD above baseline)
% in two time windows: approach and retract.
%
% Inputs:
%   rateTables      - cell array of rate tables (from buildRateTable)
%   analysis_folder - path to save results (optional)
%   mouseID         - identifier string for the animal
%   sampleRate      - sampling rate (Hz)
%
% Output:
%   responseTable   - summary table with burst info per cluster

    % === Parameters ===
    approachWindow = [1 4];     % time window for approach burst (s)
    retractWindow  = [4 7];     % time window for retract burst (s)
    threshold      = 2;         % SD threshold above baseline for burst

    nClusters = numel(rateTables);

    allClusterIDs = {};
    allStimNames  = {};
    allApproachOnset = [];
    allApproachPeak  = [];
    allApproachMaxFR = [];
    allRetractOnset  = [];
    allRetractPeak   = [];
    allRetractMaxFR  = [];

    % % Preallocate results
    % clusterID      = (1:nClusters)';
    % approach_onset = nan(nClusters,1);
    % approach_peak  = nan(nClusters,1);
    % retract_onset  = nan(nClusters,1);
    % retract_peak   = nan(nClusters,1);
    % approach_maxFR = nan(nClusters,1);
    % retract_maxFR  = nan(nClusters,1);

    % === Loop over clusters ===
    for c = 1:nClusters
        clusterTable = rateTables(c);

        % Get cluster ID (assuming it's numeric or string)
        clusterID = clusterTable.ClusterID; 
     
        %======Loop over stimuli===
        for s = 1: numel(clusterTable.StimulusLabels{1,1})
            stim = rateTables(c).StimulusLabels{1, 1}{1, s}; 
            rate_table = rateTables(c).RateTable{1,1}{1,s};

        % Get baseline mean & std (same for all rows)
        baseMean = rate_table.Baseline_mean(1);
        baseStd  = rate_table.Baseline_std(1);
        thresh   = baseMean + (baseStd*threshold);

        % --- Helper anonymous function to detect bursts in a window ---
        detectBurst = @(window) findBursts(rate_table.Time_sec, rate_table.FR_mean_Hz, window, thresh);

        % Detect bursts in each window
        [approach_onset, approach_peak, approach_maxFR] = detectBurst(approachWindow);
        [retract_onset,  retract_peak, retract_maxFR]   = detectBurst(retractWindow);

        % Append results to arrays
        allClusterIDs   = [allClusterIDs; clusterID];
        allStimNames    = [allStimNames; stim];
        allApproachOnset = [allApproachOnset; approach_onset];
        allApproachPeak  = [allApproachPeak; approach_peak];
        allApproachMaxFR = [allApproachMaxFR; approach_maxFR];
        allRetractOnset  = [allRetractOnset; retract_onset];
        allRetractPeak   = [allRetractPeak; retract_peak];
        allRetractMaxFR  = [allRetractMaxFR; retract_maxFR];
        end
    end

% === Build output table ===
responseTable = table(allClusterIDs, allStimNames, ...
    allApproachOnset, allApproachPeak, allApproachMaxFR, ...
    allRetractOnset,  allRetractPeak, allRetractMaxFR, ...
    'VariableNames', {'ClusterID','Stimulus', ...
    'ApproachOnset_s','ApproachPeak_s','ApproachMaxFR_Hz', ...
    'RetractOnset_s','RetractPeak_s','RetractMaxFR_Hz'});

    % % === Optionally save ===
    % if nargin >= 2 && ~isempty(analysis_folder)
    %     outFile = fullfile(analysis_folder, sprintf('%s_OnsetResponse.mat', mouseID));
    %     save(outFile, 'responseTable');
    %     fprintf('Saved response table to %s\n', outFile);
    %end
%end
