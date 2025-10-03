figFiles = dir(fullfile('/mnt/multiverse/homes/izzy/data/if39/day_1/analysis_recording__07_08_251/', '*.fig'));

for k = 1:length(figFiles)
    open(fullfile(figFiles(k).folder, figFiles(k).name));
end
