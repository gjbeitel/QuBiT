

% This is just case 4.

% Required inputs:
% Tube data (allFiles.XXXXX.h5bwfilled)
% Cell data as image stack?

%% Load or import tube data

if exist('allFiles','var') == 1
    fprintf('\nallFiles found in workspace.\n');
else
    fprintf('\nSelect the allFiles .mat file.\n');
    [FileName,PathName] = uigetfile('*.mat', 'Select the allFiles .mat file.', 'MultiSelect', 'off');
    load(fullfile(PathName,FileName));
    if exist('allFiles','var') == 1
        fprintf('allFiles imported successfully.\n');
    else
        error('Invalid file or unsupported data type.');
    end
end

varname = fieldnames(allFiles);
fileNo = 1;

fprintf('Reading data from allFiles...\n');
try
    skel = allFiles.(varname{fileNo}).skel;
    scurve = allFiles.(varname{fileNo}).scurve;
    scurvep = allFiles.(varname{fileNo}).scurvep;
    fprintf('Done.\n');
catch
    error('allFiles does not include all necessary data.');
end



%% import cells (objects) .h5 file

% % Use this to import from h5 (code works fine, but method might be
% % outdated.) Otherwise, it will import directly from variables.
% fprintf('Select the cell (objects) .h5 file.\n');
% [FileName,PathName] = uigetfile('*.h5', 'Select the cell (objects) .h5 file.', 'MultiSelect', 'off');
% try
%     lm12 = h5read(fullfile(PathName,FileName),'/exported_data');
% catch
%     error('Invalid file or unsupported data type.');
% end
% if numel(size(squeeze(lm12))) > 3
%     error('Objects .h5 file has more than one data channel. Be sure that "Objects" was exported correctly from the Pixel+Object Classification module in Ilastik.');
% end


try bwcellimage = allFiles.(varname{fileNo}).bwcellimage;
    fprintf('Cell image found in allFiles.\n');
catch
    if exist('bwcellimage','var') == 1
        fprintf('Cell image found in workspace.\n');
    else
        error('Cell image not found in workspace or allFiles.');
    end
    
    
end

% lm13 = squeeze(bwcellimage);
%
%
% lm19 = bwlabeln(lm13);
%
% for i = 1:max(lm19(:))
%     filter = find(lm19==i);
%     if size(filter,1) < 20
%         fprintf('Removing dust... (%4.0f of %4.0f with size %2.0f)\n',i,max(lm19(:)),size(filter,1));
%         lm13(filter) = 0;
%     end
% end
%
% bwcellimage = lm13;
% allFiles.(varname{FileNo}).bwcellimage = bwcellimage;
% save(userDefinedallFilesName,'allFiles','-append');

%% Calculate cell orientations

fprintf('Finding and plotting cell orientation vectors...\n');

% filter out objects of size 1, or you will fail to calculate 3D
% regionprops for single points


lm14 = bwlabeln(bwcellimage);
% lm17 = regionprops3(bwlabeln(bwcellimage));
lm21 = cell(max(lm14(:)),1);

% % Cell array of all cell data:

celldata = cell(max(lm14(:))+1,10);
celldata{1,1} = 'Cell ID';
celldata{1,2} = 'Centroid';
% lm99{1,3} = 'Eccentricity';
celldata{1,3} = 'Aspect Ratio';
celldata{1,4} = 'Rel. size (pix)';
celldata{1,5} = 'Orientation (degrees)';
celldata{1,6} = 'Rel. position (centerline)';
celldata{1,7} = 'Dist. from centerline (pix)';
% lm99{1,9} = 'Orientation (vector)';
celldata{1,8} = 'Rel. position (circumference)';
% The next 3 lines are used to approximate rel. pos. (circ.) We'll
% calculate the actual ref vector below.
short_axis = find(size(bwcellimage)==min(size(bwcellimage)));
ref_vector = zeros(1,3);
ref_vector(short_axis) = 1;

for i = 1:max(lm14(:))

    if nnz(lm14==i) < 20
        continue
    end
    
    fprintf('Calculating for object %4.0f of %4.0f...\n',i,max(lm14(:)));
    %     lm16 = zeros(length(lm15), 3);
    %     for mm = 1:length(lm15)
    %         [I,J,K] = ind2sub(size(bwcellimage), lm15(mm));
    %         lm16(mm,:) = [I,J,K];
    %     end
%     lm16 = false(size(bwcellimage));
%     lm16(lm15) = true;
    
    [I,J,K] = ind2sub(size(lm14),find(lm14==i));
    lm17 = regionprops3([J,I,K],'IsPixList');
    lm18 = lm17.FirstAxis; % Extract primary axis
    lm20 = lm17.Centroid; % Extract centroid
    
    % Find local tube skeleton tangent vector
    % Distance list from centroid to skel
    for j = 1:size(skel,1)
        lm21{i}(j) = pdist([skel(j,:);lm20]);
    end
    [value, index] = min(lm21{i});
    point = skel(index,:);
    if index > size(skel,1)/2
        vector = skel(index,:) - skel(index-5,:);
    else
        vector = skel(index+5,:) - skel(index,:);
    end
    % Calculate angle between cell vector and skel tangent vector and store it
    % lm22 = acosd(dot(lm18, vector) / (norm(lm18)*norm(vector)));
    % ^old code that just randomly assigned the cell vector direction
    % and therefore didn't discriminate between <90 and >90. New code
    % below assigns vector direction correctly?
    if dot(lm18,vector) < dot(-lm18,vector)
        lm18 = -lm18;
    end
    lm22 = acosd(dot(lm18, vector) / (norm(lm18)*norm(vector)));
    positive = point - lm20;
    if dot(cross(lm18,vector),positive) < 0
        lm22 = -lm22;
    end
    % Calculate reference vector for each object centroid. This is the
    % reference vector projected onto the cross-sectional plane for the
    % particular location on the tube skeleton.
    proj = ref_vector - (dot(ref_vector,vector)/norm(vector)^2)*vector;

    celldata{i+1,1} = i; % Cell ID
    celldata{i+1,2} = lm20; % Centroid
%     lm99{i+1,3} = [lm17.EigenValues(1),lm17.EigenValues(2),... Eigenvalues
%         (1-lm17.EigenValues(2)/lm17.EigenValues(1))^(1/2),... Eccentricity of major to secondary axis (1-x1/x2)^(1/2)
%         (1-(lm17.EigenValues(2)/lm17.EigenValues(1))^2)^(1/2)]; % Eccentricity (1-(x1/x2)^2)^(1/2)
    celldata{i+1,3} = lm17.SecondAxisLength/lm17.FirstAxisLength; % Aspect ratio of major axis to secondary axis
    celldata{i+1,4} = nnz(lm14==i); % Rel. size (pix)
    celldata{i+1,5} = lm22; % Cell orientation (angle) rel to local skel vector
    celldata{i+1,6} = index; % Rel. position (centerline)
    celldata{i+1,7} = value; % Cell dist. from centerline (use this to find and eliminate cells on branches)
    % lm99{i+1,9} = lm18; % Cell orientation (vector form), for visuals
    
    % Circumferential angle using the {projected reference vector on the
    % local cross-sectional plane} and the {vector formed by the local tube
    % skeleton point to the object centroid}
    circ_angle = atan2(norm(cross(proj,lm20-point)),dot(proj,lm20-point));
    % Need to determine whether this is positive or negative - can take the
    % cross product and check it against the local tangent vector on the
    % tube skeleton.
    if dot(cross(proj,lm20-point),vector) < 0
        circ_angle = -circ_angle;
    end
    celldata{i+1,8} = circ_angle;
    
end


%% Collect and filter data
% 180718: This step has moved to s8, after segments are normalized.

%% save


% fprintf('Save?\n');
% preview = input('','s');
% 
% if preview == 'Y'


allFiles.(userDefinedallFilesName).celldata = celldata;
lm99 = celldata;
    
    
    if exist('userDefinedallFilesName','var') == 0
        error('Auto-save cell data failed because a working filename was not found.');
    else
        fprintf('Saving allFiles data as %s...\n',userDefinedallFilesName);
        save(userDefinedallFilesName,'allFiles','-append');
        save([userDefinedallFilesName,'celldata'],'lm99');
        fprintf('Complete.\n');
    end
% end
