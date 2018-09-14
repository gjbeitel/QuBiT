% function allFiles = s6_handle_junctions_f(allFiles,directory,dilation_radius,dilation_factor,base_int,preview,save_dir)

%% Import cell junction data using directory of image stack

% user defined parameters:
if exist('dilation_radius','var') == 0
    dilation_radius = 1;
end
if exist('dilation_factor','var') == 0
    dilation_factor = 1;
end
if 2*dilation_radius < dilation_factor
    error('Surface parameters invalid. Check dilation_radius and dilation_value');
end



fprintf('Import junctions stack directory...\n');
folder_name = uigetdir;
D = dir(folder_name);
D2 = D(not([D.isdir]));
imNames = {D2.name}';
DirSize = length(D2);

I = imread(imNames{1});
sequence = zeros([size(I) DirSize],class(I));
sequence(:,:,1) = I;

for p = 2:DirSize
    sequence(:,:,p) = imread(imNames{p});
end
%% Debug only
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
%
%
%
%
% for ii = sequence_begin:DirSize+sequence_begin-1
%     file_name = sprintf('%s/%s%03i.tif',folder_name,sequence_pattern,ii);
%     try
%         lm12(:,:,ii-sequence_begin+1) = imread(file_name);
%     catch
%         warning('Directory size exceeds image stack size. Usually this is due to hidden files in the directory.');
%     end
% end
%
% % while range(lm12(:,:,end)) == 0
% %     lm12(:,:,end) = [];
% % end
%
% base_int = input('\nInput base intensity threshold.\n You can also enter a decimal for a percentage threshold, or\n -1 to auto-threshold by the mean intensity.\n');
%
if base_int == -1
    base_int = mean(lm12(:));
elseif (0 < base_int) && (base_int < 1)
    percent_threshold = base_int*numel(lm12);
    [~,idx] = sort(lm12(:), 'descend');
    values = lm12(idx(1:ceil(percent_threshold)));
    base_int = min(values);
end

%%
lm12 = uint8(zeros(size(sequence,1),size(sequence,2),DirSize-1));
fprintf('Image sequence imported successfully.\n');

lm13 = lm12;
lm13(lm13 <= base_int) = 0;
lm13(lm13 ~= 0) = 1;


%% Create cells through inversion of cell boundaries on apical surface and export

fprintf('Mapping cells...\n');
FileList = fieldnames(allFiles);

try xscale = allFiles.(FileList{1}).xscale;
catch
    % warning('Image X scaling info not found in allFiles. Defaulted to 1.');
    xscale = 1;
end
try yscale = allFiles.(FileList{1}).yscale;
catch
    % warning('Image Y scaling info not found in allFiles. Defaulted to 1.');
    yscale = 1;
end
try zscale = allFiles.(FileList{1}).zscale;
catch
    % warning('Image Z scaling info not found in allFiles. Defaulted to 1.');
    zscale = 1;
end

jp = lm13;

tube2 = allFiles.(FileList{1}).h5bwfilled;
tube = permute(tube2, [1 3 2]); % switch XY
expdims = zeros(ceil(size(tube,1)*zscale/xscale),ceil(size(tube,2)*zscale/yscale),size(tube,3));
cubetube = expdims; % scale to 1:1:1 based on zscale, assuming xscale=yscale
cubejuns = expdims;

if (xscale ~= 1) && (yscale ~= 1) && (zscale ~= 1)
    for i=1:size(tube,3)
        cubetube(:,:,i) = imresize(tube(:,:,i),zscale/xscale);
        cubejuns(:,:,i) = imresize(jp(:,:,i),zscale/xscale);
    end
else
    cubetube = tube;
    cubejuns = jp;
end


ctouter = imdilate(cubetube,ones(dilation_radius*2+dilation_factor,dilation_radius*2+dilation_factor,dilation_radius*2+dilation_factor));
ctinner = imdilate(cubetube,ones(dilation_radius*2-dilation_factor,dilation_radius*2-dilation_factor,dilation_radius*2-dilation_factor));
ctshell = ctouter - ctinner;

% cubejuns = imdilate(cubejuns,ones(3,5,5));
% cubejuns = smooth3(cubejuns); % 170217: attempting this with initial smoothing to reduce bg noise and comp time

apicaljunctions = and(ctshell,cubejuns);

[shellx, shelly, shellz] = ind2sub(size(ctshell), find(ctshell == 1));
for i = 1:size(shellx,1)
    if apicaljunctions(shellx(i),shelly(i),shellz(i)) == 1
        apicaljunctions(shellx(i),shelly(i),shellz(i)) = 0;
    else
        apicaljunctions(shellx(i),shelly(i),shellz(i)) = 1;
    end
end

bwcellimage = permute(apicaljunctions,[1 3 2]);

allFiles.(FileList{1}).bwcellimage = bwcellimage;

%
% fprintf('Preview cell image and tube?\n');
% preview = input('','s');

if strcmp(preview,'Y')
    figure;axis equal;axis([-inf inf -inf inf -inf inf]);hold on;
    blah = isosurface(allFiles.(userDefinedallFilesName).h5bwfilled, 0.5);
    patch(blah,'facecolor','red','facealpha',0.3,'edgecolor','none');
    blah = isosurface(allFiles.(userDefinedallFilesName).bwcellimage, 0.3);
    patch(blah, 'facecolor', 'green', 'facealpha', 0.3, 'edgecolor', 'none');
end

%% save

if strcmp(overwrite,'Y')
    if size(tube2) ~= size(bwcellimage)
        warning('Dimensions of tube and cell image do not match. Usually this is a dim permutation issue.');
    end
    fprintf('Workspace and base file (%s.mat) updated.\n',FileList{1});
    save(FileList{1},'allFiles','-v7.3','-append');
else
    fprintf('Workspace updated, but base file was not overwritten.\n');
end

return

%% Debug only

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

if ~isstr(save_dir)
    mkdir(FileList{1})
else
    mkdir(save_dir)
end
for ii = 1:size(bwcellimage, 3)
    imwrite(bwcellimage(:,:,ii), sprintf('image%04d.jpg',ii));
    movefile(sprintf('image%04d.jpg',ii), FileList{1});
end

fprintf('\nStack has successfuly been saved to the folder %s.\nUse Pixel+Object Classification module in Ilastik to map cells. Then return to step 4\n', FileList{1});
% end