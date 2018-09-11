%% import/aggregate data for tube stat plots for a particular genotype
% this should come from the allFiles data, e.g.
% <genotype><replicate_number>T.mat

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

length_data = 1:9;
cross_data = 1:9;
perim_data = 1:9;
cs_aggregate = [];
for i = 1:numofFiles
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    % load necessary data into workspace
    
    try branchid = allFiles.(varname(1:end-4)).branchid;
    catch
        warning('File does not contain branch normalization data. Skipping...');
        continue
    end
    
    try seglen = allFiles.(varname(1:end-4)).seglen;
        segcs = allFiles.(varname(1:end-4)).segcs;
        segper = allFiles.(varname(1:end-4)).segper;

    missing_index = 0;
    for j = 1:size(branchid,2)
        if j+1 > size(branchid,2)
            break
        end
        
        fprintf('  Iteration %i of %i...\n',j,size(branchid,2));
        % check segment length, if <1, normalize, if >1, skip.
        if branchid(j+1) - branchid(j) > 1
            fprintf('     Segment length exceeds 1. Interpolating...\n');
            missing_segments = ceil(branchid(j+1) - branchid(j) - 1);
            seglen(:,j+missing_index+1+missing_segments:size(seglen,2)+missing_index+missing_segments) = seglen(:,j+1:end);
            
                         segcs(:,j+missing_index+1+missing_segments:size(segcs,2)+missing_index+missing_segments) = segcs(:,j+1:end);
                         segcs(:,j+missing_index+1:j+missing_index+missing_segments) = zeros(size(segcs,1),missing_segments);
                         segper(:,j+missing_index+1+missing_segments:size(segper,2)+missing_index+missing_segments) = segper(:,j+1:end);
                         segper(:,j+missing_index+1:j+missing_index+missing_segments) = zeros(size(segper,1),missing_segments);
            
                         val = ceil(seglen(1,j+missing_index))-seglen(1,j+missing_index);
                         if val == 0
                             val = 1;
                         end
                         base_el = size(segcs,1)*val/(branchid(j+1)-branchid(j));
            
            k = 1;
            while k <= missing_segments
                seglen(1,j+missing_index+k) = floor(seglen(1,j+missing_index))+k;
                seglen(2,j+missing_index+k) = seglen(2,j+missing_index)/(branchid(j+1)-branchid(j));
                
                                 move_el = size(segcs,1)/(branchid(j+1)-branchid(j));
                                 if round(base_el+move_el*k) > size(segcs,1)
                                     range_el = round(base_el+move_el*(k-1)):size(segcs,1);
                                 else
                                 range_el = round(base_el+move_el*(k-1)):round(base_el+move_el*k);
                                 end
                                 segcs(1:numel(range_el),j+missing_index+k) = segcs(range_el,j+missing_index);
                                 segper(1:numel(range_el),j+missing_index+k) = segper(range_el,j+missing_index);
                
                %segcs(:,j+missing_index+k) = segcs(:,j+missing_index);
                %segper(:,j+missing_index+k) = segper(:,j+missing_index);
                k = k + 1;
            end
            
            seglen(2,j+missing_index) = seglen(2,j+missing_index)*val/(branchid(j+1)-branchid(j));
            
                         segcs(round(base_el):end,j+missing_index) = zeros(size(segcs,1)-round(base_el)+1,1);
                         segcs(all(segcs==0,2),:)=[];
                         segper(round(base_el):end,j+missing_index) = zeros(size(segper,1)-round(base_el)+1,1);
                         segper(all(segper==0,2),:)=[];
                         
            missing_index = missing_segments+missing_index;
            % j = j+missing_segments;
            
        elseif branchid(j+1) - branchid(j) < 1
            fprintf('     Extrapolating segment length...\n');
            seglen(2,j) = seglen(2,j)/(branchid(j+1) - branchid(j));
        else
        end
        
        % j = j+1;
    end
    while seglen(1,end) < seglen(1,end-1)
        seglen(:,end) = [];
    end
    while size(segcs,2) > size(seglen,2)
        segcs(:,end) = [];
    end
    while size(segper,2) > size(seglen,2)
        segper(:,end) = [];
    end
    % compile into data arrays
    for j = 1:size(seglen,2)
        if isempty(find(length_data(:,floor(seglen(1,j)))==0, 1)) == 1
            length_data(end+1,floor(seglen(1,j))) = seglen(2,j);
        else
            length_data(find(length_data(:,floor(seglen(1,j)))==0, 1 ),floor(seglen(1,j))) = seglen(2,j);
        end
    end
    for j = 1:size(seglen,2)
        try cross_data(find(cross_data(:,floor(seglen(1,j)))==0, 1 ):...
                find(cross_data(:,floor(seglen(1,j)))==0, 1 )+size(segcs,1)-1,floor(seglen(1,j))) = segcs(1:end,j);
            perim_data(find(perim_data(:,floor(seglen(1,j)))==0, 1 ):...
                find(perim_data(:,floor(seglen(1,j)))==0, 1 )+size(segper,1)-1,floor(seglen(1,j))) = segper(1:end,j);
            % 180531
        catch
            cross_data(end+1:end+size(segcs,1),floor(seglen(1,j))) = segcs(1:end,j);
            perim_data(end+1:end+size(segper,1),floor(seglen(1,j))) = segper(1:end,j);
        end
    end
    
    catch
        warning('File does not contain all necessary tube segment components. Skipping...');
        continue
    end
    
    
    try csdata = allFiles.(userDefinedallFilesName).csdata;
        cs_aggregate = [cs_aggregate;cell2mat(csdata([false;cell2mat(csdata(2:end,2))~=0],[2 4 5]))];
    catch
        warning('CS v2 data not found. Skipping...');
    end
    
end

length_data(length_data==0)=NaN;
cross_data(cross_data==0)=NaN;
perim_data(perim_data==0)=NaN;

% find and remove outliers? am i justified in doing this?
length_data([zeros(1,size(length_data,2),'logical');isoutlier(length_data(2:end,:),1)]==1)=NaN;
cross_data([zeros(1,size(cross_data,2),'logical');isoutlier(cross_data(2:end,:),1)]==1)=NaN;
perim_data([zeros(1,size(perim_data,2),'logical');isoutlier(perim_data(2:end,:),1)]==1)=NaN;


if ~isempty(save_name)
save(save_name,'length_data','cross_data','perim_data','cs_aggregate');    
end




%% 180531

[file_array, ~, ~] = uigetfile('*.mat', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
end

s_segcs = nan(10,10,numofFiles);
s_norpo = nan(10,10,numofFiles);
for i = 1:numofFiles
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    % load necessary data into workspace
    
    try branchid = allFiles.(varname(1:end-4)).branchid;
    catch
        warning('File does not contain branch normalization data. Skipping...');
        continue
    end
    try segcs = allFiles.(varname(1:end-4)).segcs;
    catch
        warning('File does not contain cross-section data. Skipping...');
        continue
    end
    
    normalized_position = zeros(size(segcs));
    for j = 1:size(branchid,2)
        if j+1 > size(branchid,2)
            break
        end
        
        normalized_position(:,j) = (branchid(j)+(branchid(j+1)-branchid(j))/11:(branchid(j+1)-branchid(j))/11:branchid(j+1)-(branchid(j+1)-branchid(j))/11)';
    end
    
    if size(normalized_position) ~= size(segcs)
        warning('Debug: Normalized position array is not the same size as cross-section data array in file %s.',i);
    end
    s_segcs(1:size(segcs,1),1:size(segcs,2),i) = segcs;
    s_norpo(1:size(normalized_position,1),1:size(normalized_position,2),i) = normalized_position;
    
end
fprintf('Done.\n');

% save... length_data cross_data perim_data as <genotype>tubedata.mat