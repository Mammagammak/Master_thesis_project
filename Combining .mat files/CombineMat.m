function combined_results = CombineMat(JSD)
%% this function combines the .mat files from different experiments for further data analysis
oldFolder = cd(JSD);

fList = dir('*.mat'); 

fList = {fList.name}'; %'// Gets only the name info from the mat files and makes the list from them

%% // Processing the list of files:
if isempty(fList), combined_results = []; return, end %// Check that some files exist

%// Initialize the result struct by loading the last file (also preallocates the struct):
nFiles = size(fList,1);
% combined_results(nFiles) = load(fullfile(fList{nFiles}));

%// See if there is only 1 file, and return if so:
if nFiles == 1, return, end

%// Process any additional files 
for ind1 = 1:nFiles
    combined_results(ind1) = load(fullfile(fList{ind1}));
end

%% Cleanup - changing the current directory back:
cd(oldFolder);