%% Bulk reapply filters

[file_array, ~, ~] = uigetfile('*.mat', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
end

if w1118 == true
if stage == 14
    % Filter by size, range: 1um to 4um, calculated from hist analysis
    size_min = (1/.151)^2; size_max = (4/.151)^2; 
    % Filter by radial dist. (max allowed values, conical A->P), range calculated from 2um^2 to 12um^2 surface area
    filter_start = 4.1585*1.2+1; filter_end = 11.5924*1.2+1;
elseif stage == 16
    % Filter by size, range: 1.5um to 6um, calculated from hist analysis
    size_min = (2/.151)^2; size_max = (6.5/.151)^2;
    % Filter by radial dist. (max allowed values, conical A->P). range calculated form 10um^2 to 50um^2 surface area
    filter_start = sqrt(10/.151^2/pi)+5; filter_end = sqrt(30/.151^2/pi)+5;
end
end

if src42 == true
if stage == 14
    % Filter by size, range: 0.6um to 3um, calculated from hist analysis
    size_min = (0.6/.151)^2; size_max = (3/.151)^2; 
    % Filter by radial dist. (max allowed values, conical A->P), range calculated from 2um^2 to 12um^2 surface area
    filter_start = sqrt(2/.151^2/pi)+5; filter_end = sqrt(12/.151^2/pi)+5;
elseif stage == 16
    % Filter by size, range: 1.3um to 5.5um, calculated from hist analysis
    size_min = (1.3/.151)^2; size_max = (5.5/.151)^2;
    % Filter by radial dist. (max allowed values, conical A->P). range calculated form 10um^2 to 50um^2 surface area
    filter_start = sqrt(10/.151^2/pi)+5; filter_end = sqrt(50/.151^2/pi)+5;
end
end

fprintf('--Filter parameters:--\nw1118: %i\nsrc42: %i\nStage: %i\n',w1118,src42,stage);
fprintf('--Details:--\nSize: [%2.1f, %2.1f]\nRadial: [%2.1f, %2.1f]\n',size_min,size_max,filter_start,filter_end);
for j = 1:numofFiles
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(j)),j,numofFiles);
        load(char(file_array(j)))
        varname = char(file_array(j));
        load([varname(1:end-12),'.mat'])
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,j,numofFiles);
        load(file_array);
        varname = file_array;
        load([varname(1:end-12),'.mat'])
    end
    % load necessary data into workspace
    
    try
        lm98 = cell2mat(lm99(2:end,:));
    catch
        warning('Data formatted incorrectly.');
        continue
    end


    % lm98 = cell2mat(lm99(2:end,:));
    filter_max = max(cell2mat(lm99(2:end,1)));
    
    % Filter by rel. pos.

    lm98(lm98(:,10)==1,:) = [];
    lm98(lm98(:,10)==max(lm98(:,10)),:) = [];

    % Filter by size
    
    lm98(lm98(:,11)<=size_min,:) = [];
    lm98(lm98(:,11)>=size_max,:) = [];

    
    i = 1;
    while i <= size(lm98,1)
        
        filter_value = filter_start - (filter_start-filter_end)*lm98(i,1)/filter_max;
        if lm98(i,12) >= filter_value
            lm98(i,:) = [];
        else
            i = i+1;
        end
        
    end

    
    if size(lm98,1) == 0
        warning('No data found after filtering! Skipping...');
        continue
    end
    
    fprintf('Creating filtered cell image...\n');
    lm14 = bwconncomp(allFiles.(varname(1:end-12)).bwcellimage);
bwcellimagef = false(size(allFiles.(varname(1:end-12)).bwcellimage));
for i = 1:size(lm98,1)
    bwcellimagef(lm14.PixelIdxList{lm98(i,1)}) = true;
end
allFiles.(varname(1:end-12)).bwcellimagef = bwcellimagef;
    
    fprintf('Saving %s...\n',char(file_array(j)));
    save(char(file_array(j)),'lm98','-append');
    save([varname(1:end-12),'.mat'],'allFiles','-append');
    
%     fprintf('Plot filtered cells 3D...\n');
% lm14 = bwlabeln(allFiles.(userDefinedallFilesName).bwcellimagef);
% figure();axis equal;axis([-inf inf -inf inf -inf inf]);hold on;
% fprintf('Plotting cells...\n');
% colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
%     [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
%     [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1];[0.2 0.5 0.8];[0.2 0.8 0.5];...
%     [0.8 0.2 0.5];[0.8 0.5 0.2];[0.5 0.8 0.2];[0.5 0.2 0.8]];
% for kk = 1:max(lm14(:))
% % for kk = 1:5
%     fprintf('Plotting object with index %3.0f (%3.0f of %3.0f)...\n',lm98(kk,1),kk,max(lm14(:)));
%     % lm15 = lm14.PixelIdxList{lm98(kk,1)};
%     % lm16 = false(size(allFiles.(userDefinedallFilesName).bwcellimage));
%     % lm16(lm15) = true;
%     blah = isosurface(lm14==kk,0.5);
%     patch(blah, 'facecolor', colorspec(mod(kk,size(colorspec,1))+1,:), 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
% end
    
end

fprintf('Complete.\n');