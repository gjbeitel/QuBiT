function allFiles = s4_calc_tube_length_f(allFiles)

if ~isstruct(allFiles)
    error('The input argument must be a struct of the H5 data.')
end

FileList = fieldnames(allFiles);

%% Gather segment piece dist values
S = allFiles.(FileList{1}).S;
skelnums = allFiles.(FileList{1}).skelnums;

clear segmentpieces
segmentpieces(1,:) = skelnums;
startindex = 0;
for i = 1:size(skelnums,2)
    pieces = size(S{skelnums(i)},1)-1;
    fprintf('Gathering segment %i from pieces [%i to %i]...\n',segmentpieces(1,i),startindex+1,startindex+pieces);
    % segmentpieces(3,i) = sum(lengthpieces(startindex+1:startindex+pieces));
    startindex = startindex+pieces;
    segmentpieces(2,i) = startindex;
end
allFiles.(FileList{1}).segmentpieces = segmentpieces;

return

%% Debug only

%% Save

if strcmp(overwrite,'Y')
    fprintf('Base file (%s.mat) updated.\n',FileList{1});
    save(FileList{1},'allFiles','-v7.3','-append');
end

%% Debug only

% allFilesnames = fieldnames(allFiles);
% if exist('FileList{1}','var') == 0
%     error('Working file name not found. Unable to auto-save.');
% else
%     allFilescheck = strcmp(FileList{1},allFilesnames);
%     if sum(allFilescheck(:)) == 1
%         fprintf('Saving cell image to allFiles in %s...\n',FileList{1});
%         save(FileList{1},'allFiles','-v7.3','-append');
%         fprintf('allFiles has been saved to %s successfully.\n', FileList{1});
%     else
%         warning('The working file name and save file name do not match. Skipping auto-save...');
%     end
% end


% end
%
%
% % If the user had decided to save the variable, the variable will
% % be saved here.  Note that the -v7.3 is required in the save
% % function due to the probable size of the data structure.
%
%
%
% % clearvars -except allFiles cmd2004 figstab FileList{1} varname savecmd

end
