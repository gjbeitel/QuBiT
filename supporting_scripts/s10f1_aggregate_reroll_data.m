%% import/aggregate reroll data for a particular genotype
% this should come from the celldata, e.g.
% <genotype><replicate_number>reroll.mat

% aggregate data



[file_array, ~, ~] = uigetfile('*.mat', 'Select .mat file(s).', 'MultiSelect', 'on');

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

reroll_agg = [];

for i = 1:numofFiles
    clear reroll_data
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)));
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    
    
    
    aggregate = cell2mat(reroll_data(2:end,[6 4 3 5]));
    aggregate(isnan(aggregate(:,1)),:)=[];
    
    reroll_agg = [reroll_agg;aggregate];
    
    
end

save(save_name,'reroll_agg');

fprintf('Done.\n');
