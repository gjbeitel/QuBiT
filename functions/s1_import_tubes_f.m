function allFiles = s1_import_tubes_f

%% Import datasets

fprintf('Select the HDF5 files.\n');
[FileName, PathName, ~] = uigetfile('*.h5', 'Select the HDF5 files.', 'MultiSelect', 'on');

if FileName == 0
    error('No files selected.');
end

if ischar(FileName)
    FileName = cellstr(FileName);
end

for i = 1:length(FileName)
    [~, name, ~] = fileparts(FileName{i});
    name = matlab.lang.makeValidName(name);
    allFiles.(name).h5import = h5read(fullfile(PathName, FileName{i}), '/exported_data');
end

%% Batch Skeleton

FileList = fieldnames(allFiles);
for i = 1:length(FileName)
    
    fprintf('Managing File %i of %i...\n', i, length(FileName));
    
    
    allFiles.(FileList{i}).h5p = squeeze(allFiles.(FileList{i}).h5import(1,:,:,:));
    
    % Have not been able to test if this works yet...
    if (range(allFiles.(FileList{i}).h5p(:)) ~= 0)
        
        allFiles.(FileList{i}).h5bw = uint8(permute(allFiles.(FileList{i}).h5p, [3 2 1]));
        allFiles.(FileList{i}).h5bwfilled = zeros(size(allFiles.(FileList{i}).h5bw));
        
        for ii = 1:size(allFiles.(FileList{i}).h5bw,3)
            allFiles.(FileList{i}).h5bwfilled(:,:,ii) = imfill(allFiles.(FileList{i}).h5bw(:,:,ii));
        end
        
        allFiles.(FileList{i}).h5bw = permute(allFiles.(FileList{i}).h5bwfilled, [1 3 2]);  % x z y
        allFiles.(FileList{i}).h5bwfilled = zeros(size(allFiles.(FileList{i}).h5bw));
        for ii = 1:size(allFiles.(FileList{i}).h5bw,3)
            allFiles.(FileList{i}).h5bwfilled(:,:,ii) = imfill(allFiles.(FileList{i}).h5bw(:,:,ii));
        end
        
        allFiles.(FileList{i}).h5bw = permute(allFiles.(FileList{i}).h5bwfilled, [3 2 1]);  % y z x
        allFiles.(FileList{i}).h5bwfilled = zeros(size(allFiles.(FileList{i}).h5bw));
        for ii = 1:size(allFiles.(FileList{i}).h5bw,3)
            allFiles.(FileList{i}).h5bwfilled(:,:,ii) = imfill(allFiles.(FileList{i}).h5bw(:,:,ii));
        end
        
        allFiles.(FileList{i}).h5bwfilled = permute(allFiles.(FileList{i}).h5bwfilled, [1 3 2]); % y x z
        % allFiles.(varname{fileNo}).h5bwfilled = permute(allFiles.(varname{fileNo}).h5bwfilled, [3 1 2]); % x y z
        % At some point, use this to avoid further confusion with permute
        
        allFiles.(FileList{i}).h5bwfilled = logical(allFiles.(FileList{i}).h5bwfilled);
        
        fprintf('\nConstructing the skeleton for file %i...\n', i);
        allFiles.(FileList{i}).S = skeleton(allFiles.(FileList{i}).h5bwfilled);
        
        % Flip x-y dimensions in S{i}
        for ii = 1:numel(allFiles.(FileList{i}).S)
            a = allFiles.(FileList{i}).S{ii};
            a(:,[1 2])=a(:,[2 1]);
            allFiles.(FileList{i}).S{ii} = a;
        end
    else
        warning('Blank stack in %s. Skipping...\n', FileList{i});
    end
    
end

%% Save

splitallfiles(allFiles);

fprintf('Data saved.\n');

end

%% Supporting function(s)

function splitallfiles(allFiles)

a = fieldnames(allFiles);

if length(a) > 1
    
    d = allFiles;
    
    for i = 1:length(a)
        fprintf('Saving %s (file %i of %i)...\n',a{i},i,length(a));
        clear allFiles
        allFiles.(a{i}) = d.(a{i});
        save(a{i},'allFiles','-v7.3');
        clear allFiles
    end
    
else
    
    save(a{1}, 'allFiles', '-v7.3');
    
end

end