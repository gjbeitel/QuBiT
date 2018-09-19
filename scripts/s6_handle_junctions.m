%% Import cell junction data using directory of image stack

% user defined parameters:
if exist('dilation_radius','var') == 0
    dilation_radius = 1;
end
if exist('dilation_value','var') == 0
    dilation_value = 1;
end
if 2*dilation_radius < dilation_value
    error('Surface parameters invalid. Check dilation_radius and dilation_value');
end


% case 3
fprintf('\nImport junctions stack directory...\n');
folder_name = uigetdir;
D = dir(folder_name);
D2 = D(not([D.isdir]));
DirSize = length(D2);

% 
% try
%     test_file = sprintf('%s/imj0000.tif',folder_name);
%     test = imread(test_file);
%     sequence_pattern = 'imj0';
%     sequence_begin = 0;
% catch
% 
% 
% sequence_pattern = input('Input image sequence pattern.\nFor example, "T00001C01Z" for raw images or "imj0" (no quotes).\n','s');
% sequence_begin = input('Input stack start value. This is usually 0 or 1.\n');
% 
% 
% test_file = sprintf('%s/%s%03i.tif',folder_name,sequence_pattern,sequence_begin);
% try
%     test = imread(test_file);
% catch
%     error('Unable to import image sequences.');
% end
% 
% end

fprintf('\nWorking...\n');
I = imread([D2(1).folder,'/',D2(1).name]);
lm12 = uint8(zeros(size(I,1),size(I,2),DirSize));

for i = 1:DirSize
    lm12(:,:,i) = imread([D2(i).folder,'/',D2(i).name]);
end

% for ii = sequence_begin:DirSize+sequence_begin-1
%     file_name = sprintf('%s/%s%03i.tif',folder_name,sequence_pattern,ii);
%     try
%         lm12(:,:,ii-sequence_begin+1) = imread(file_name);
%     catch
%         warning('Directory size exceeds image stack size. Usually this is due to hidden files in the directory.');
%     end
% end

% while range(lm12(:,:,end)) == 0
%     lm12(:,:,end) = [];
% end

fprintf('\nImage sequence imported successfully.');

base_int = input('\nQuick intensity threshold segmentation:\n');

if base_int == -1
    base_int = mean(lm12(lm12~=0));
elseif (0 < base_int) && (base_int < 1)
    percent_threshold = base_int*numel(lm12);
    [~,idx] = sort(lm12(:), 'descend');
    values = lm12(idx(1:ceil(percent_threshold)));
    base_int = min(values);
end


lm13 = lm12;
lm13(lm13 <= base_int) = 0;
lm13(lm13 ~= 0) = 1;

%% Create cells through inversion of cell boundaries on apical surface and export

fprintf('Mapping cells...\n');
varname = fieldnames(allFiles);
fileNo = 1;

% try xscale = allFiles.(varname{fileNo}).xscale;
% catch
%     warning('Image X scaling info not found in allFiles. Defaulted to 1.');
%     xscale = 1;
% end
% try yscale = allFiles.(varname{fileNo}).yscale;
% catch
%     warning('Image Y scaling info not found in allFiles. Defaulted to 1.');
%     yscale = 1;
% end
% try zscale = allFiles.(varname{fileNo}).zscale;
% catch
%     warning('Image Z scaling info not found in allFiles. Defaulted to 1.');
%     zscale = 1;
% end

juns = lm13;
tube = permute(allFiles.(varname{fileNo}).h5bwfilled, [1 3 2]); % switch XY
% expdims = zeros(ceil(size(tube,1)*zscale/xscale),ceil(size(tube,2)*zscale/yscale),size(tube,3));
% cubetube = expdims;
% cubejuns = expdims;
% 
% if (xscale ~= 1) && (yscale ~= 1) && (zscale ~= 1)
%     for i=1:size(tube,3)
%         cubetube(:,:,i) = imresize(tube(:,:,i),zscale/xscale);
%         cubejuns(:,:,i) = imresize(jp(:,:,i),zscale/xscale);
%     end
% else
%   cubetube = tube;
%   cubejuns = jp;
% end


ctouter = imdilate(tube,ones(dilation_radius*2+dilation_value,dilation_radius*2+dilation_value,dilation_radius*2+dilation_value));
ctinner = imdilate(tube,ones(dilation_radius*2-dilation_value,dilation_radius*2-dilation_value,dilation_radius*2-dilation_value));
ctshell = ctouter - ctinner;

% cubejuns = imdilate(cubejuns,ones(3,5,5));
% cubejuns = smooth3(cubejuns); % 170217: attempting this with initial smoothing to reduce bg noise and comp time

% juns_inv = imcomplement(juns);
% [sx,sy,sz] = ind2sub(size(ctshell),find(ctshell==1));
% apicaljunctions = juns_inv(sx,sy,sz);

apicaljunctions = and(ctshell,juns);

[shellx, shelly, shellz] = ind2sub(size(ctshell), find(ctshell == 1));
for i = 1:size(shellx,1)
    if apicaljunctions(shellx(i),shelly(i),shellz(i)) == 1
        apicaljunctions(shellx(i),shelly(i),shellz(i)) = 0;
    else
        apicaljunctions(shellx(i),shelly(i),shellz(i)) = 1;
    end
end

bwcellimage = permute(apicaljunctions,[1 3 2]);

if exist('FileNo','var') == 0
    FileNo = 1;
end
allFiles.(varname{FileNo}).bwcellimage = bwcellimage;


fprintf('Preview cell image and tube?\n');
preview = input('','s');

if preview == 'Y'
    figure;axis equal;axis([-inf inf -inf inf -inf inf]);hold on;
    blah = isosurface(allFiles.(userDefinedallFilesName).h5bwfilled, 0.5);
    patch(blah,'facecolor','red','facealpha',0.3,'edgecolor','none');
    blah = isosurface(allFiles.(userDefinedallFilesName).bwcellimage, 0.3);
    patch(blah, 'facecolor', 'green', 'facealpha', 0.3, 'edgecolor', 'none');
end

%% save

fprintf('Save?\n');
preview = input('','s');

if preview == 'Y'
allFilesnames = fieldnames(allFiles);
if exist('userDefinedallFilesName','var') == 0
    warning('Working file name not found. Unable to auto-save.');
else
    allFilescheck = strcmp(userDefinedallFilesName,allFilesnames);
    if sum(allFilescheck(:)) == 1
        if size(tube2) ~= size(bwcellimage)
            warning('Dimensions of tube and cell image do not match. Usually this is a dim permutation issue.');
        end
        fprintf('Saving cell image to allFiles in %s...\n',userDefinedallFilesName);
        save(userDefinedallFilesName,'allFiles','-v7.3','-append');
    else
        warning('The working file name and save file name do not match. Skipping auto-save...');
    end
end
end

return

%% old crap

% figure;
% blah = isosurface(apicaljunctions, 0.3);
% patch(blah, 'facecolor', 'green', 'facealpha', 0.3, 'edgecolor', 'none');

% return % Use this return to display attempted cells without
% % saving the data to the folder


% This creates a folder for the jpegs that are created.  The
% folder is known as the original file name, but matlab
% acceptable.  The names of the jpegs are attempted to be named
% in order of the z-stack, therefore adding zeros in front of
% each jpeg is necessary, i.e. the number 1 is 001.  These
% jpegs will be allowed to import into ilastik.

mkdir(varname{fileNo})
for ii = 1:size(bwcellimage, 3)
    imwrite(bwcellimage(:,:,ii), sprintf('image%04d.jpg',ii));
    movefile(sprintf('image%04d.jpg',ii), varname{fileNo});
end

fprintf('\nStack has successfuly been saved to the folder %s.\nUse Pixel+Object Classification module in Ilastik to map cells. Then return to step 4\n', varname{fileNo});