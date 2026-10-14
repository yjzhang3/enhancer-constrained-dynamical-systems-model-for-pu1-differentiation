function [bestValues,bestXX,minLastValue] = read_best_job_arr
% Current directory
currentDir = pwd;

% Get only the current folder name
[~, currentFolderName] = fileparts(currentDir);

% Construct the target folder name
folderPath = fullfile(currentDir, [currentFolderName '_NOTSoutputs']);

% Find matching files
files = dir(fullfile(folderPath, 'ODENOTSfit_PS_*.txt'));

minLastValue = inf;
bestXX = '';
bestValues = [];

for k = 1:numel(files)

    filename = files(k).name;
    filepath = fullfile(files(k).folder, filename);

    % Read file contents
    txt = strtrim(fileread(filepath));

    % Convert comma-separated values to numeric array
    values = str2double(strsplit(txt, ','));

    % Last value
    lastValue = values(end);

    if lastValue < minLastValue
        minLastValue = lastValue;

        % Store the numerical array
        bestValues = values;

        % Extract xx from filename
        token = regexp(filename, '^ODENOTSfit_PS_(\d+)\.txt$', 'tokens');
        bestXX = token{1}{1};
    end
end

bestValues = bestValues(1:end-1);
