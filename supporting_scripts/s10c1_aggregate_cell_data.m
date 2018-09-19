%% import/aggregate 3d cell data for a particular genotype
% this should come from the celldata, e.g.
% <genotype><replicate_number>celldata.mat

% aggregate 3d cell data

[file_array, ~, ~] = uigetfile('*.*', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
end

fprintf('\nSave aggregate data as? (Leave blank to skip saving)\n');
save_name = input('','s');

% file_array = input('Files array {string1 ,string2...}: ');
lm97 = [];
for i = 1:numofFiles
    clear lm98
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    if size(lm98,2) ~= 11
        warning('Missing data. Usually this means it was not normalized. Skipping...');
        continue
    end
    lm97 = [lm97;lm98];
end
clear lm98;

% remove out of bounds cells
lm97(lm97(:,11)==0,:)=[];

% normalize orientations
lm97(lm97(:,7)>90,7)=180-lm97(lm97(:,7)>90,7);

if ~isempty(save_name)
save(save_name,'lm97');    
end

fprintf('Done.\n');

% save... lm97 as <genotype>celldata.mat