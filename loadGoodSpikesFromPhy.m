function S = loadGoodSpikesFromPhy(dataPath)
% S = loadGoodSpikesFromPhy
% Load in and output data about the clusters labelled "good" in phy. 
% Input the path to the Phy output directory.
% Outputs a 1xN struct with fields ClusterID and SpikeTimes. Spike times
% are in terms of samples. 

% IN PROGRESS
% - May add additional outputs later
% GHP March 2022
% 
% GHP June 2022: Adding amplitude as an output. This is the amplitude read
% directly from Phy, and NOT in any particular units. 
if nargin == 0
    dataPath = uigetdir(ksortPath);
end
allSpikeTimes = readNPY(fullfile(dataPath,'spike_times.npy'));
allSpikeClusters = readNPY(fullfile(dataPath,'spike_clusters.npy'));
[clusterIDs,clusterGroups] = readClusterGroupsCSV(fullfile(dataPath,'cluster_group.tsv'));
goodClusters = clusterIDs(clusterGroups == 2);

try
    allAmplitudes = readNPY(fullfile(dataPath,'amplitudes.npy'));
catch
    warning('amplitudes.npy not found on path')
    allAmplitudes = nan(size(allSpikeTimes));
end
nClusts = length(goodClusters);
S = struct;
for i = 1:nClusts
    S(i).ClusterID = goodClusters(i);
    S(i).SpikeTimes = allSpikeTimes(allSpikeClusters == goodClusters(i));
    S(i).Amplitudes = allAmplitudes(allSpikeClusters == goodClusters(i));
end