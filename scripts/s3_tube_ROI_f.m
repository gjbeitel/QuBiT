function allFiles = s3_tube_ROI_f(allFiles)
% Code Requirements and Initial Actions
%
% Code requires the use of myupdatefcn.m
% unalteredcombined.m contains test trial version of combined.


if ~isstruct(allFiles)
    error('The input argument must be a struct of the H5 data.')
end


%% Combineskeleton.m and labeling on plot

FileList = fieldnames(allFiles);





[skel, skelnums] = combineskeleton(allFiles.(FileList{1}).S); % Function call to combine skeleton
% figure('name', FileList{1}); hold on;
% plot3(skel(:,1),skel(:,2),skel(:,3),'black', 'LineWidth', 1.5);
% text(skel(1,1),skel(1,2),skel(1,3),'   START OF SKELETON');
% text(skel(end,1),skel(end,2),skel(end,3), '   END OF SKELETON');

%     % PROMPT:  Select alternate start point?
%     cmd4 = '?';
%     while (~strcmp(cmd4,'Y') && ~strcmp(cmd4,'N'))
%         cmd4 = input('\nSelect alternate start point?  Y/N\n','s');
%         if isempty(cmd4)
%             cmd4 = '?';
%         end
%     end
%
%     switch cmd4
%         case 'Y'
%             fprintf('Instructions for selecting alternate point:\n');
%             fprintf('1) Click on a point so a coordinate text box shows up.\n');
%             fprintf('2) Right-click on the text box.\n');
%             fprintf('3) Select "Export Cursor Data to Workspace".\n');
%             fprintf('4) Save data as variable called "cursor_info" (default).\n');
%             fprintf('5) Press any key to continue.\n');
%
%             fprintf('\nNote:  Rotating the image will de-select datacursormode.\n');
%             fprintf('If this should happen, simply click on the icon on the toolbar directly to right of rotation icon.\n');
%             fprintf('(Should look like a line with a plus sign and a yellow box.)\n\n');
%             fprintf('Note:  If wrong point is selected, simply click on the correct point and continue.\n\n');
%
%             datacursormode on
%             pause
%             dcm_obj = datacursormode(fig);
%             start = getCursorInfo(dcm_obj);
%             start = start.Position;
%
%             % Note:  Coordinates at this point are already switched on plot.
%             % Do NOT switch them again.
%             plot3(start(1,1),start(1,2),start(1,3), '*', 'MarkerSize', 12);
%             text(start(1,1),start(1,2),start(1,3), '   ALTERNATE START');
%
%         otherwise
%             start = skel(1,:);
%     end
%
%     % PROMPT:  Select alternate end point?
%     cmd5 = '?';
%     while (~strcmp(cmd5,'Y') && ~strcmp(cmd5,'N'))
%         cmd5 = input('\nSelect alternate end point? Y/N\n','s');
%         if isempty(cmd5)
%             cmd5 = '?';
%         end
%     end
%
%     switch cmd5
%         case 'Y'
%             fprintf('Instructions for selecting alternate point:\n');
%             fprintf('1) Click on a point so a coordinate text box shows up.\n');
%             fprintf('2) Right-click on the text box.\n');
%             fprintf('3) Select "Export Cursor Data to Workspace".\n');
%             fprintf('4) Save data as variable called "cursor_info" (default).\n');
%             fprintf('5) Press any key to continue.\n');
%
%             fprintf('\nNote:  Rotating the image will de-select datacursormode.\n');
%             fprintf('If this should happen, simply click on the icon on the toolbar directly to right of rotation icon.\n');
%             fprintf('(Should look like a line with a plus sign and a yellow box.)\n\n');
%             fprintf('Note:  If wrong point is selected, simply click on the correct point and continue.\n\n');
%
%             datacursormode on
%             pause
%             dcm_obj = datacursormode(fig);
%             ender = getCursorInfo(dcm_obj);
%             ender = ender.Position;
%
%             % Note:  Coordinates at this point are already switched on plot.
%             % Do NOT switch them again.
%             plot3(ender(1,1),ender(1,2),ender(1,3), '*', 'MarkerSize', 12);
%             text(ender(1,1),ender(1,2),ender(1,3), '   ALTERNATE END');
%         otherwise
%             ender = skel(end,:);
%     end

%% Adjusting skeleton matrix to alternate ends (if necessary)
% This section makes Alternateends.m unnecessary
%     if strcmp(cmd4,'Y') || strcmp(cmd5,'Y')
%         temp = ismember(skel, start, 'rows');
%         ind = find(temp == 1);
%         skel = skel(ind:end,:);
%
%         temp = ismember(skel, ender, 'rows');
%         ind = find(temp == 1);
%         skel = skel(1:ind,:);
%
%         figure
%         hold on
%         plot3(skel(:,1), skel(:,2), skel(:,3), 'black');
%
%         % Note:  Coordinates at this point are already switched on plot.
%         % Do NOT switch them again.
%         plot3(start(1,1),start(1,2),start(1,3), '*', 'MarkerSize', 12);
%         text(start(1,1),start(1,2),start(1,3),['START', char(10), 'X: ', num2str(start(1,1)), char(10), 'Y: ', num2str(start(1,2)), char(10), 'Z: ', num2str(start(1,3))]);
%         plot3(ender(1,1),ender(1,2),ender(1,3), '*', 'MarkerSize', 12);
%         text(ender(1,1),ender(1,2),ender(1,3),['END', char(10), 'X: ', num2str(ender(1,1)), char(10), 'Y: ', num2str(ender(1,2)), char(10), 'Z: ', num2str(ender(1,3))]);
%     end


%% Append

allFiles.(FileList{1}).skel = skel;
allFiles.(FileList{1}).skelnums = skelnums;


fprintf('Creating the spline...\n');
splinedata = skel(1:3:end,:);
allFiles.(FileList{1}).scurve = cscvn(splinedata.');
syms t;
allFiles.(FileList{1}).scurvep = fnder(allFiles.(FileList{1}).scurve);


%% Gather segment piece dist values
% redundant
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
    fprintf('Workspace and base file (%s.mat) updated.\n',FileList{1});
    save(FileList{1},'allFiles','-v7.3','-append');
else
    fprintf('Workspace updated, but base file was not overwritten.\n');
end


%% Debug only

% allFilesnames = fieldnames(allFiles);
% if exist('userDefinedallFilesName','var') == 0
%     error('Working file name not found. Unable to auto-save.');
% else
%     allFilescheck = strcmp(userDefinedallFilesName,allFilesnames);
%     if sum(allFilescheck(:)) == 1
%         fprintf('Saving to allFiles in %s...\n',userDefinedallFilesName);
%         save(userDefinedallFilesName,'allFiles','-v7.3','-append');
%         fprintf('allFiles has been saved to %s successfully.\n', userDefinedallFilesName);
%     else
%         warning('The working file name and save file name do not match. Skipping auto-save...');
%     end
% end
%





% Editor Notes:
%        cursor_info.Position = skel(end,:); % Set Default (possibly incorrect)

%         datacursormode on
%         dcm_obj = datacursormode(fig);
%         set(dcm_obj,'UpdateFcn',@myupdatefcn);
%         pause
%         ender = cursor_info.Position;

%        ender = skel(end,:);

end

%% Supporting functions

function [skel, skelnums] = combineskeleton(S)

fprintf('\nSegment Numbers Must Be Input In Order\n');
cmd = 0;
skel = [];
skelnums = [];
segments = 0;
done = 100+rand;
restart = 101+rand;

while (cmd ~= done)
    fprintf('-------------------------\n');
    cmd = input('Input Segment # -OR- "restart" -OR- "done"... ');
    if isempty(cmd)
        cmd = 0;
    end
    if (cmd >= 1 && cmd <= length(S))
        if isempty(skel)
            skel = S{cmd};
            segments = segments + 1;
            skelnums(segments) = cmd;
            
            fprintf('\nCurrent Skeleton: \n');
            disp(skelnums);
            
        elseif (segments == 1)
            
            if (any(ismember(skelnums, cmd)))
                cmd2 = '?';
                while (~strcmp(cmd2,'Y') && ~strcmp(cmd2,'N'))
                    cmd2 = input('Input segment has already been used.  Are you sure you want to use this segment again? Y or N?\n','s');
                    if isempty(cmd2)
                        cmd2 = '?';
                    end
                end
                
                if (cmd2 == 'N')
                    fprintf('\nCurrent Skeleton: \n');
                    disp(skelnums);
                    continue
                end
            end
            
            head = skel(1,:);
            tail = skel(end,:);
            if (S{cmd}(1,:) == head) % head of new to head of old
                skel = flipud(skel);
                skel = skel(1:end-1,:);
                skel = vertcat(skel,S{cmd});
                segments = segments + 1;
                skelnums(segments) = cmd;
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            elseif (S{cmd}(1,:) == tail) % head of new to tail of old
                skel = skel(1:end-1,:);
                skel = vertcat(skel,S{cmd});
                segments = segments + 1;
                skelnums(segments) = cmd;
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            elseif (S{cmd}(end,:) == head) % tail of new to head of old
                skel = flipud(skel);
                skel = skel(1:end-1,:);
                skel = vertcat(skel,flipud(S{cmd}));
                segments = segments + 1;
                skelnums(segments) = cmd;
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            elseif (S{cmd}(end,:) == tail) % tail of new to tail of old
                skel = skel(1:end-1,:);
                skel = vertcat(skel,flipud(S{cmd}));
                segments = segments + 1;
                skelnums(segments) = cmd;
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            else
                fprintf('Input segment not connected to skeleton\n');
                
                cmd10 = '?';
                while (~strcmp(cmd10,'Y') && ~strcmp(cmd10,'N'))
                    cmd10 = input('Connect anyway? Y or N?\n','s');
                    if isempty(cmd10)
                        cmd10 = '?';
                    end
                end
                
                if (cmd10 == 'N')
                    fprintf('\nskelnums = \n');
                    disp(skelnums);
                    continue
                end
                
                % If Yes
                a = pdist([tail;(S{cmd}(1,:))]);
                b = pdist([tail;(S{cmd}(end,:))]);
                c = pdist([head;(S{cmd}(1,:))]);
                d = pdist([head;(S{cmd}(end,:))]);
                
                if (a == min([a b c d])) % head of new to tail of old
                    skel = skel(1:end-1,:);
                    skel = vertcat(skel,S{cmd});
                    segments = segments + 1;
                    skelnums(segments) = cmd;
                elseif (b == min([a b c d])) % tail of new to tail of old
                    skel = skel(1:end-1,:);
                    skel = vertcat(skel,flipud(S{cmd}));
                    segments = segments + 1;
                    skelnums(segments) = cmd;
                elseif (c == min([a b c d])) % head of new to head of old
                    skel = flipud(skel);
                    skel = skel(1:end-1,:);
                    skel = vertcat(skel,S{cmd});
                    segments = segments + 1;
                    skelnums(segments) = cmd;
                elseif (d == min([a b c d])) % tail of new to head of old
                    skel = flipud(skel);
                    skel = skel(1:end-1,:);
                    skel = vertcat(skel,flipud(S{cmd}));
                    segments = segments + 1;
                    skelnums(segments) = cmd;
                end
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            end
            
            clear head
            clear tail
            
        else
            
            if (any(ismember(skelnums, cmd)))
                cmd2 = '?';
                while (~strcmp(cmd2,'Y') && ~strcmp(cmd2,'N'))
                    cmd2 = input('Input segment has already been used.  Are you sure you want to use this segment again? Y or N?\n','s');
                    if isempty(cmd2)
                        cmd2 = '?';
                    end
                end
                
                if (cmd2 == 'N')
                    fprintf('\nskelnums = \n');
                    disp(skelnums);
                    continue
                end
            end
            
            tail = skel(end,:);
            if (S{cmd}(1,:) == tail) % head of new to tail of old
                skel = skel(1:end-1,:);
                skel = vertcat(skel,S{cmd});
                segments = segments + 1;
                skelnums(segments) = cmd;
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            elseif (S{cmd}(end,:) == tail) % tail of new to tail of old
                skel = skel(1:end-1,:);
                skel = vertcat(skel,flipud(S{cmd}));
                segments = segments + 1;
                skelnums(segments) = cmd;
                
                fprintf('\nCurrent Skeleton: \n\n');
                disp(skelnums);
            else
                fprintf('Input segment not connected to skeleton\n');
                
                cmd10 = '?';
                while (~strcmp(cmd10,'Y') && ~strcmp(cmd10,'N'))
                    cmd10 = input('Connect anyway? Y or N?\n','s');
                    if isempty(cmd10)
                        cmd10 = '?';
                    end
                end
                
                if (cmd10 == 'N')
                    fprintf('\nskelnums = \n');
                    disp(skelnums);
                    continue
                end
                
                % If Yes
                a = pdist([tail;(S{cmd}(1,:))]);
                b = pdist([tail;(S{cmd}(end,:))]);
                
                if (a < b) % head of new to tail of old
                    skel = skel(1:end-1,:);
                    skel = vertcat(skel,S{cmd});
                    segments = segments + 1;
                    skelnums(segments) = cmd;
                else % tail of new to tail of old
                    skel = skel(1:end-1,:);
                    skel = vertcat(skel,flipud(S{cmd}));
                    segments = segments + 1;
                    skelnums(segments) = cmd;
                end
                
                fprintf('\nCurrent Skeleton: \n');
                disp(skelnums);
            end
            
            clear tail
            
        end
    elseif (cmd == restart)
        skel = [];
        skelnums = [];
        segments = 0;
        
        fprintf('\nCurrent Skeleton: \n\n\t[empty]\n');
    elseif (cmd == done)
        continue
    else
        fprintf('Invalid Input\n\n');
        fprintf('\nCurrent Skeleton: \n');
        disp(skelnums);
    end
    
end

end
