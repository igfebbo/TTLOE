%% Input recording parameters

%This information needed to build file paths from the excel sheet
mouseID='if46';
day=2;
experiment=1;
rec=1;
excelFilePath='/mnt/multiverse/homes/izzy/data/Mouse_Bible.xlsx';
recNum=num2str(rec); %converts recording number to a string for saving files later
sampleRate=30000;

%% Input TTL and stimulus information. 
% Note: TTL name and number, and stimuli parameters may seem redundant,
% but they are kept separate because in some instances ,TTLs are referenced as hardware inputs, and stimuli are referenced for plotting etc.
% Define all TTLs that were used during recording.
% nums= the port number that your TTL was plugged into in OpenEphys
% The name 'blockTTL' is used to splice a recording into different analysis segments

TTL_definitions = struct( ...
    'nums',  {1,3,5,6}, ...
    'names', {'orangeTTL','blockTTL','smoothTTL','roughTTL'} ...
);

TTLs = numel(TTL_definitions);
% Define analysis/plotting properties for TTLs that you wish to appear as stimuli on your plots. 
% 
% Duration, delay and alignment window are in seconds.
% 
% If you are going to align your spikes to a stimulus, then the delay for
% that stimulus is zero. Define other stimulus delays relative to your
% alignment stimulus
% 
% alpha is the transparency of the shaded region that will overlay on your plot during the stimulus
%
%for onORoff, syntax needs to start with an underscore, no spaces. example: _on, _off
stimuliParams = struct( ...
    'name', {'orangeTTL','smoothTTL','roughTTL'}, ...
    'label', {'orange laser','Smooth Texture','Rough Texture'}, ...
    'duration', {0.05,4,4}, ...
    'delay', {0,0,0}, ...
    'color', {[1,0.5,0],'b','r'}, ...
    'alpha', {0.4,0.1,0.1}, ...
    'alignmentWindow', {[-1.0,1.0],[-1.0,6.0],[-1.0,6.0]}, ...
    'onORoff', {'_on', '_off','_off'} ...
    );
    
%% Define file paths and create analysis folder

[ksortPath, continuousPath, acquisitionPath, ttlFolderPath, dataPath, savePath] = filePaths(mouseID, day, experiment, rec, excelFilePath);


% Create 'analysis' folder if it doesn't exist
analysis_folder= fullfile(savePath, ['analysis_',recNum]);
if ~exist(analysis_folder, 'dir')
    mkdir(analysis_folder);
end

%% Create 'session' object that holds information about your recordings, does not load data
% run only if using openEphys MATLAB analysis tools. Otherwise unnecessary.
% Useful if you want to analyze multiple recordings at once.

%directory = dataPath; 
%session = Session(directory);

%% Loading and aligning TTLs and recordings

spikedata=loadGoodSpikesFromPhy(ksortPath);
states=readNPY(fullfile(ttlFolderPath,'states.npy'));
sample_numbers_TTL=readNPY(fullfile(ttlFolderPath,'sample_numbers.npy'));
sample_numbers_cont=readNPY(fullfile(acquisitionPath,'sample_numbers.npy'));
full_words=readNPY(fullfile(ttlFolderPath,'full_words.npy'));
sampleNumbersTTL_zeroed=sample_numbers_TTL-sample_numbers_cont(1);
n_clusters = length(spikedata); % Number of clusters in our data set

for i = 1:TTLs
    name = TTL_definitions(i).names;  % e.g., 'orangeTTL'
    
    on_field  = [name '_on'];
    off_field = [name '_off'];

    TTL_struct.(on_field)  = sampleNumbersTTL_zeroed(states == TTL_definitions(i).nums);
    TTL_struct.(off_field) = sampleNumbersTTL_zeroed(states == -TTL_definitions(i).nums);
end

%% Create and save spikeMatrices and firing rate tables aligned to different stimuli (no block splitting)
% raw firing rate is basleine subtracted

%[SpikeMatrices,rateTables] = plotClusterRasters(spikedata, TTL_struct.orangeTTL_on, sampleRate, plottingWindow, stimuliParams, '' , analysis_folder);

%% Create and save raster plots 

%% Create and save spikeMatrices and firing rate tables aligned to different TTLs per block (blocks are sectioned by a TTL signal, default: TTL3)
%set desired bin width
binWidth=0.05;

% Define blocks: user specifies labels and which TTL fields to align to. Can add multiple. Example: blocks 1–2 align to single TTL, blocks 3–4 align to two stimuli
% One alignment window per block.

userBlocks = struct( ...
    'label', { ...
        'Rough and Smooth Texture Response', ...
        'Baseline-No laser', ...
        'orange 1mW', 'orange 2mw', 'orange 4mw', ...
        'orange 6mW', 'orange 7mw', 'orange 10mw', 'orange 12mw' ...
    }, ...
    'alignTTLs', { ...
        {'smoothTTL_off','roughTTL_off'}, ... % Block 1 aligns to OFFs
        {'orangeTTL_on'}, ...                 % Block 2 aligns to ON
        {'orangeTTL_on'}, {'orangeTTL_on'}, {'orangeTTL_on'}, ...
        {'orangeTTL_on'}, {'orangeTTL_on'}, {'orangeTTL_on'}, {'orangeTTL_on'} ...
    } ...
);

nblocks = numel(userBlocks);

if nblocks == 1
    blockEdges = [0;sample_numbers_cont(end)];
else
blockEdges = [0; TTL_struct.blockTTL_on];
end

%nBlocks = numel(blockEdges) - 1;

blockParams = struct([]);

for b = 1:nBlocks
    blockParams(b).label = userBlocks(b).label;
    blockParams(b).edges = [blockEdges(b), blockEdges(b+1)];
    
    % alignment TTLs (user defined explicitly)
    ttlNames = userBlocks(b).alignTTLs;
    blockParams(b).alignmentTTLnames = userBlocks(b).alignTTLs;
    blockParams(b).alignmentTTLs = cell(size(ttlNames));
    blockParams(b).alignmentWindows = cell(size(ttlNames));
    
    
    for t = 1:numel(ttlNames)
        ttlName = ttlNames{t};
        if ~isfield(TTL_struct, ttlName)
            error('TTL_struct does not contain field "%s"', ttlName);
        end
        blockParams(b).alignmentTTLs{t} = TTL_struct.(ttlName)(TTL_struct.(ttlName)>=blockEdges(b) & TTL_struct.(ttlName)<blockEdges(b+1));

        % look up corresponding alignment window from stimuliParams
        stimIdx = find(strcmp({stimuliParams.name}, erase(ttlName, {'_on','_off'})));
        if isempty(stimIdx)
            error('No stimuliParams entry found for "%s"', ttlName);
        end
        blockParams(b).alignmentWindows{t} = stimuliParams(stimIdx).alignmentWindow*sampleRate;
    end
end

[SpikeMatrices, rateTables] = clustersByBlock(spikedata, blockParams,sampleRate, binWidth);
 
%% Save rateTables and spikeMAtrices

save(fullfile(analysis_folder, [mouseID '_SpikeMatrices.mat']), 'SpikeMatrices', '-v7.3');
save(fullfile(analysis_folder, [mouseID '_rateTables.mat']), 'rateTables', '-v7.3');

%% Create rasters and PSTHs from spikeMatrices and firing rate tables aligned to different TTLs per block

for i = 1:length(spikedata) 
plotClusterRastersAndPSTHs(i, SpikeMatrices, rateTables, stimuliParams,blockParams, TTL_struct, sampleRate, analysis_folder);
end

%% Create and save heatmap from rateTables (must have rate tables loaded)

%heatMap()