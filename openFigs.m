figFiles = dir(fullfile('/mnt/multiverse/homes/izzy/data/if46/day_2/analysis_1/figures', '*.fig'));

for k = 1:length(figFiles)
    open(fullfile(figFiles(k).folder, figFiles(k).name));
end
