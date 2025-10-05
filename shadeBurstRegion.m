function [approachOnset, retractOnset] = shadeBurstRegion(responseTable, clusterID, stim, nTrials, shadeDur, shift)
    approachOnset = NaN;
    retractOnset  = NaN;

    if nargin < 5
        shadeDur = 0.010; 
    end

    % Find matching row
    row = find(strcmp(string(responseTable.ClusterID), string(clusterID)) & ...
               strcmp(string(responseTable.Stimulus), string(stim)), 1);
    if isempty(row)
        return; 
    end

    % Approach burst
    onset = responseTable.ApproachOnset_s(row) + shift;
    if ~isnan(onset)
        offset = onset + shadeDur;
        patch([onset offset offset onset], [0 nTrials nTrials 0], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        approachOnset = onset;  % return for debugging
    end

    % Retract burst
    rOnset = responseTable.RetractOnset_s(row) + shift;
    if ~isnan(rOnset)
        rOffset = rOnset + shadeDur;
        patch([rOnset rOffset rOffset rOnset], [0 nTrials nTrials 0], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        retractOnset = rOnset;
    end
end
