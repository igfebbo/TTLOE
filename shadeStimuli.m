function [legendHandles, legendLabels] = shadeStimuli(stimuliParams, TTL_struct, alignTTLs, alignWindow, yMax)
    hold on;
    legendHandles = [];
    legendLabels = {};

    for s = 1:length(stimuliParams)
        ttlField = [stimuliParams(s).name stimuliParams(s).onORoff];
        stimTimesThisBlock = TTL_struct.(ttlField);
        startTime = alignTTLs(1) + alignWindow(1);
        endTime   = alignTTLs(1) + alignWindow(2);
        
        % Keep only stimulus events within this trial's alignment window
        stimTimesPerStim = stimTimesThisBlock(stimTimesThisBlock >= startTime & stimTimesThisBlock <= endTime);

        % stimTimesThisBlock = stimTimesThisBlock( ...
        %     stimTimesThisBlock >= blockParams.edges(1) & ...
        %     stimTimesThisBlock <= blockParams.edges(end));
        % 
        if ~isempty(stimTimesPerStim)
            shadeStart = stimuliParams(s).delay;
            shadeEnd   = shadeStart + stimuliParams(s).duration;

            % Call original function
            handleStruct = makeShadedRect(shadeStart, shadeEnd, 0, yMax, ...
                                          stimuliParams(s).color, stimuliParams(s).alpha);

                        % --- Create a small transparent legend "dot" (proxy) ---
            % Try scatter first (preferred, supports MarkerFaceAlpha)
            try
                hLegend = scatter(nan, nan, 80, 'MarkerFaceColor', stimuliParams(s).color, ...
                                  'MarkerEdgeColor', 'none');
                % set marker alpha if supported
                if isprop(hLegend, 'MarkerFaceAlpha')
                    hLegend.MarkerFaceAlpha = stimuliParams(s).alpha;
                end
            catch
                % Fallback: small circular patch as legend proxy
                theta = linspace(0, 2*pi, 30);
                r = 0.08; % radius of legend symbol (small)
                x = r*cos(theta);
                y = r*sin(theta);
                hLegend = patch(x, y, stimuliParams(s).color, ...
                                'FaceAlpha', stimuliParams(s).alpha, ...
                                'EdgeColor', 'none', ...
                                'Visible', 'off'); % invisible in axes
            end

            % store the legend handle & label
            legendHandles(end+1) = hLegend; %#ok<AGROW>
            legendLabels{end+1} = stimuliParams(s).label; %#ok<AGROW>


            % theta = linspace(0, 2*pi, 30); % circle
            % x = cos(theta) * 0.1;  % scale radius
            % y = sin(theta) * 0.1;
            % 
            % hDot = patch(x, y, stimuliParams(s).color, ...
            %  'FaceAlpha', stimuliParams(s).alpha, ...
            %  'EdgeColor', 'none', ...
            %  'Visible', 'off');  % don't show in axes
            % 
            % for k = 1:length(handleStruct)
            %     legendHandles(end+1) = hDot; %handleStruct(k).P;
            %     legendLabels{end+1} = stimuliParams(s).label;
            end
end


