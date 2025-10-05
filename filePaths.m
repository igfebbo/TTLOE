function [ksortPath, continuousPath, acquisitionPath, ttlFolderPath, dataPath, savePath] = filePaths(mouseID, day, experiment, recording, excelFilePath)
    % This function generates the file path for continuous data and TTL signal based on inputs.
    
    % Read the Excel sheet into MATLAB
    data = readtable(excelFilePath);

    % Filter the table to find the row matching the user inputs
    filteredData = data(strcmp(data.mouse_id, mouseID) & data.day == day, :);

    % Check if the filtered data is empty
    if isempty(filteredData)
        error('No matching data found for mouse ID: %s, day: %d', mouseID, day);
    end

    % Extract the full_file value
    fullFile = filteredData.full_file{1};  % Assuming there's one match; adjust if multiple matches are possible

    % Format the inputs into appropriate folder names
    dayFolder = sprintf('day_%d', day);
    experimentFolder = sprintf('experiment%d', experiment);
    recordingFolder = sprintf('recording%d', recording);

    % Define the root directory of your data files
    rootDir = '/mnt/multiverse/homes/kathi/data';  % Adjust to your data directory
    
    % Construct the kilosort output file path
    ksortPath = fullfile(rootDir, mouseID, dayFolder, fullFile, ...
        'Record Node 101', experimentFolder, recordingFolder, ...
        'continuous', 'Acquisition_Board-100.Rhythm Data', 'kilosort4');

    % Construct the continuous file path
    continuousPath = fullfile(rootDir, mouseID, dayFolder, fullFile, ...
        'Record Node 101', experimentFolder, recordingFolder, ...
        'continuous', 'Acquisition_Board-100.Rhythm Data', 'continuous.dat');

    % Construct the continuous file path
    acquisitionPath = fullfile(rootDir, mouseID, dayFolder, fullFile, ...
        'Record Node 101', experimentFolder, recordingFolder, ...
        'continuous', 'Acquisition_Board-100.Rhythm Data');

    % Construct the TTL signal file path
    ttlFolderPath = fullfile(rootDir, mouseID, dayFolder, fullFile, ...
        'Record Node 101', experimentFolder, recordingFolder, ...
        'events', 'Acquisition_Board-100.Rhythm Data', 'TTL');
    
    
    % Construct the data folder path
    dataPath = fullfile(rootDir, mouseID, dayFolder, fullFile,'/');

    % Construct the data folder path
    savePath = fullfile(rootDir, mouseID, dayFolder);

%     % Output the continuous file path
% fprintf('The file path is:\n%s\n', continuousPath);
% 
% % Output the TTL file path
% fprintf('The file path is:\n%s\n', ttlFolderPath);

% Check if the kilosort folder exists
if isfolder(ksortPath)
    fprintf('Kilosort folder exists at the specified path.\n');
else
    fprintf('Kilosort folder does not exist at the specified path.\n');
end

% Check if the continuous file exists
if isfile(continuousPath)
    fprintf('Continuous file exists at the specified path.\n');
else
    fprintf('Continuous file does not exist at the specified path.\n');
end

% Check if the TTL file exists
if isfolder(ttlFolderPath)
    fprintf('TTL folder exists at the specified path.\n');
else
    fprintf('TTL folder does not exist at the specified path.\n');
end
% Check if the data folder exists
if isfolder(dataPath)
    fprintf('data folder exists at the specified path.\n');
else
    fprintf('data folder does not exist at the specified path.\n');
end