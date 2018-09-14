


%% ID branches

try allFiles.(userDefinedallFilesName).segmentpieces
catch
    warning('Segment length data not found in allFiles.');
end

try branchid = allFiles.(userDefinedallFilesName).branchid;
catch
    % return % remove this later
branchid = input('Branch ID array (e.g. [1 2 3 4 5]): ');
end

try branchposition = allFiles.(userDefinedallFilesName).branchposition;
catch
    % return % remove this later
branchposition = input('Position boundary array (e.g. [1 201 403 607 812]): ');
end

csdata = allFiles.(userDefinedallFilesName).csdata;


bwcellimage = allFiles.(userDefinedallFilesName).bwcellimage;



%% Normalize

celldataf = cell2mat(celldata(2:end,:));

if size(branchid,2) ~= size(branchposition,2)
    error('Array dim mismatch.');
end

k = 0;

celldataf(:,11) = zeros(size(celldataf,1),1);
for i = 1:size(celldataf,1)
    if celldataf(i,8) < branchposition(1)
        k = k + 1;
    else
        j = 1;
        while celldataf(i,8) > branchposition(j+1)
            
            j = j + 1;
            if j >= size(branchposition,2)
                break
            end
            
        end
        try celldataf(i,11) = (celldataf(i,8)-branchposition(j))/(branchposition(j+1)-branchposition(j))*(branchid(j+1)-branchid(j))+branchid(j);
        catch
            k = k + 1;
        end
    end
end


%% Apply filters

% filter params:
if exist('stage','var') == 0
    stage = 16;
end
if stage == 14
    fprintf('Using stage 14 filters...\n');
    % Filter by size, range: 1um to 4um, calculated from hist analysis
    % expected cell size, and manual cell counting/removal of cojoined
    % cells resulting from poor segmentation
    size_min = (1/.151)^2; size_max = (4/.151)^2;
    % Filter by radial dist. (max allowed values, conical A->P), range calculated from 2um^2 to 12um^2 surface area
    filter_start = 10; filter_end = 18;
else
    fprintf('Using stage 16 filters...\n');
    % Filter by size, range: 2um to 6.5um, calculated from hist analysis
    % expected cell size, and manual cell counting/removal of cojoined
    % cells resulting from poor segmentation
    size_min = (2/.151)^2; size_max = (6.5/.151)^2;
    % Filter by radial dist. (max allowed values, conical A->P). range calculated form 10um^2 to 100um^2 surface area
    filter_start = 16; filter_end = 28.9;
end

% execute filters
fprintf('Applying data filters...\n');


% Filter by rel. pos.
fprintf('Applying filters:\n-----'); 
fprintf('\nX position filter (ends)...\n');
celldataf(celldataf(:,8)==1,:) = [];
celldataf(celldataf(:,8)==max(celldataf(:,8)),:) = [];
a = length(celldata)-1;
b = length(celldataf);
fprintf('Objects filtered: %i of %i (%i%%)\n',a-b,a,round((a-b)/a*100));

% Filter by size
if exist('dilation_value','var') == 0
    dilation_value = 1;
end
fprintf('\nSize filter: [%2.1f, %2.1f]...\n',size_min*dilation_value,size_max*dilation_value);
fprintf(' (Factoring in object thickness %i)\n',dilation_value)
celldataf(celldataf(:,6)<=size_min*dilation_value,:) = [];
celldataf(celldataf(:,6)>=size_max*dilation_value,:) = [];
c = size(celldataf,1);
fprintf('Objects filtered: %i of %i (%i%%)\n',b-c,a,round((b-c)/a*100));

fprintf('\nN position filter (conical): [%2.1f, %2.1f]...\n',filter_start,filter_end);
i = 1;
while i <= size(celldataf,1)
    
    filter_value = filter_start + (filter_end-filter_start)*(celldataf(i,11)-1)/9;
    if celldataf(i,9) > filter_value
        celldataf(i,:) = [];
    else
        i = i+1;
    end
    
end
d = length(celldataf);
fprintf('Objects filtered: %i of %i (%i%%)\n',c-d,a,round((c-d)/a*100));

if size(celldataf,1) == 0
    warning('No data found after filtering!');
end

fprintf('-----\nTotal objects remaining: %i of %i (%i%%)\n',d,a,round(d/a*100));

fprintf('Creating filtered cell image...\n');
lm14 = bwlabeln(bwcellimage);
bwcellimagef = zeros(size(bwcellimage));
for i = 1:size(celldataf,1)
    bwcellimagef(lm14==celldataf(i,1)) = celldataf(i,1);
end




%% save


allFiles.(userDefinedallFilesName).celldataf = celldataf;
lm98 = celldataf;
allFiles.(userDefinedallFilesName).bwcellimagef = bwcellimagef;


save(userDefinedallFilesName,'allFiles','-append');
save([userDefinedallFilesName,'celldata'],'lm98','-append');
fprintf('Normalization data saved to %s.\n',[userDefinedallFilesName,'celldata']);

% use this return to just recalculate cell data without redoing tube calcs.

return

%% calc individual segment length

fprintf('Calculating lengths of individual TC segments...\n');
skel = allFiles.(userDefinedallFilesName).skel;
syms t

seglen = zeros(1,size(branchposition,2)-1);
for i = 1:size(branchposition,2)-1
    fprintf('TC segment: %1.1f\n',branchid(i));
    temp = branchposition(i):3:branchposition(i+1);
    if temp(end) ~= branchposition(i+1)
        temp(end) = branchposition(i+1);
    end
    skel_array = skel(temp,:);
    spline = cscvn(skel_array.');
    dspline = fnder(spline);
    
    seglens = zeros(spline.pieces,1);
    for jj = 1:spline.pieces
        fxt=spline.coefs(3*jj-2,:);
        fyt=spline.coefs(3*jj-1,:);
        fzt=spline.coefs(3*jj,:);
        fx = fxt(1)*(t^3) + fxt(2)*(t^2) + fxt(3)*t + fxt(4);
        fy = fyt(1)*(t^3) + fyt(2)*(t^2) + fyt(3)*t + fyt(4);
        fz = fzt(1)*(t^3) + fzt(2)*(t^2) + fzt(3)*t + fzt(4);
        redifx = dspline.coefs(3*jj-2,1)*t^2 + dspline.coefs(3*jj-2,2)*t + dspline.coefs(3*jj-2,3); % derivatives used for length calculation
        redify = dspline.coefs(3*jj-1,1)*t^2 + dspline.coefs(3*jj-1,2)*t + dspline.coefs(3*jj-1,3);
        redifz = dspline.coefs(3*jj,1)*t^2 + dspline.coefs(3*jj,2)*t + dspline.coefs(3*jj,3);
        
        sl = int(sqrt(redifx^2 + redify^2 + redifz^2),...
            t, 0, spline.breaks(jj+1) - spline.breaks(jj)); % distance formula
        seglens(jj)=double(sl); % get rid of syms - changes value to a double
    end
    seglen(i) = sum(seglens);
end

seglen(2,:) = branchid(1,1:end-1);
seglen([1,2],:) = seglen([2,1],:);


%% calculate csdata

fprintf('Normalizing CS data...\n');

if size(csdata,2) == 5 && strcmp(csdata(1,2),'Normalized Position')
    warning('Normalized data already exists in %s. Skipping...',userDefinedallFilesName);
    
elseif size(csdata,2) == 4
    raw_positions = cell2mat(csdata(2:end,1));
    norm_values = zeros(1,length(raw_positions));
    for i = 1:length(branchid)-1
        norm_values(raw_positions>branchposition(i) & raw_positions<branchposition(i+1)) = (raw_positions(raw_positions>branchposition(i) & raw_positions<branchposition(i+1)) - branchposition(i))/(branchposition(i+1)-branchposition(i))*(branchid(i+1)-branchid(i))+branchid(i);
    end
    csdata(:,3:5) = csdata(:,2:4);
    csdata{1,2} = 'Normalized Position';
    csdata(2:end,2) = num2cell(norm_values);
    fprintf('Complete.\n');
    
else
    warning('Data missing or invalid in %s.. Skipping...',userDefinedallFilesName);
end



%% calculate cross sectional area (10 data points per segment)

% 
% fprintf('Calculating cross-section...\n');
% segcs = zeros(10,size(branchid,2)-1);
% segper = zeros(10,size(branchid,2)-1);
% csperseg = zeros(10,size(branchid,2)-1);
% h5bwdata = double(allFiles.(userDefinedallFilesName).h5bwfilled);
% % bwcellimage = double(allFiles.(userDefinedallFilesName).bwcellimage);
% % bwcellimage2 = zeros(size(bwcellimage));
% % 
% % lm14 = bwconncomp(bwcellimage);
% % for h = lm98(2:end,1)
% % bwcellimage2(lm14.PixelIdxList{h}) = 1;
% % end
% 
% 
% for i = 1:size(branchposition,2)-1
%     fprintf('TC segment: %i\n',branchid(i));
%     temp = branchposition(i):3:branchposition(i+1);
%     if temp(end) ~= branchposition(i+1)
%         temp(end) = branchposition(i+1);
%     end
%     skel_array = skel(temp,:);
%     spline = cscvn(skel_array.');
%     
%     %     rf = round(branchid(i+1)-branchid(i));
%     %     if rf == 0
%     %         warning('This segment is too short. Skipping.');
%     %         continue
%     %     end
%     
%     range = spline.breaks(end)/11:spline.breaks(end)/11:10*spline.breaks(end)/11;
%     cv = fnval(spline,range);
%     cdv = fnval(fnder(spline),range);
%     
%     univec = zeros(3, size(cdv,2));
%     
%     for jj = 1:size(cdv, 2)
%         univec(:,jj) = cdv(:,jj)/norm(cdv(:,jj));
%     end
%     
%     pixwidth = 30; % SUBJECT TO CHANGE
%     
%     
%     
%     subplrow = -1;
%     
%     if size(cv,2) > 10
%         for bb = 10:-1:2
%             if mod(size(cv,2), bb) == 0
%                 subplrow = bb;
%                 subplcol = size(cv,2)/bb;
%                 break;
%             end
%         end
%     end
%     
%     if subplrow == -1
%         subplrow = 1;
%         subplcol = size(cv,2);
%     end
%     
%     
%     for qq = 1:size(cv,2)
%         
%         boxloc = [max([cv(1,qq) - pixwidth,1]), min([cv(1,qq) + pixwidth,floor(size(h5bwdata, 2))]), ...
%                   max([cv(2,qq) - pixwidth,1]), min([cv(2,qq) + pixwidth,floor(size(h5bwdata, 1))]), ...
%                   max([cv(3,qq) - pixwidth,1]), min([cv(3,qq) + pixwidth,floor(size(h5bwdata, 3))])];
%         
% %         boxloc = [cv(1,qq) - pixwidth, cv(1,qq) + pixwidth, ...
% %             cv(2,qq) - pixwidth, cv(2,qq) + pixwidth, ...
% %             cv(3,qq) - pixwidth, cv(3,qq) + pixwidth];
% %         
% %         % Changed from boxloc(5) <= 0
% %         if boxloc(5) < 1
% %             boxloc(5) = 1;
% %         end
% %         if boxloc(3) <  1
% %             boxloc(3) = 1;
% %         end
% %         if boxloc(1) <  1
% %             boxloc(1) = 1;
% %         end
% %         
% %         if boxloc(6) >= size(h5bwdata, 3)
% %             boxloc(6) = floor(size(h5bwdata, 3));
% %         end
% %         if boxloc(4) >= size(h5bwdata, 1)
% %             boxloc(4) = floor(size(h5bwdata, 1));
% %         end
% %         if boxloc(2) >= size(h5bwdata, 2)
% %             boxloc(2) = floor(size(h5bwdata, 2));
% %         end
%         
%         
%         h5bwdatabox = h5bwdata(floor(boxloc(3)):ceil(boxloc(4)), floor(boxloc(1)): ceil(boxloc(2)), floor(boxloc(5)):ceil(boxloc(6)));
%         
%         [xind, yind, zind] = ind2sub(size(h5bwdatabox), find(h5bwdatabox == 1));
%         xyzind = [yind + floor(boxloc(1)) - 1, xind + floor(boxloc(3)) - 1, zind + floor(boxloc(5)) - 1];
%         
%         vdist2plane = zeros(size(xyzind, 1), 1);
%         
%         for kk = 1:size(xyzind, 1)
%             vdist2plane(kk) = dot(univec(:,qq), xyzind(kk,:)' - cv(:,qq));
%         end
%         
%         v = abs(vdist2plane);
%         sortind = find(v(:) <= 1);
%         
%         vlist = zeros(size(sortind,1),3);
%         
%         for kk = 1:size(sortind, 1)
%             vlist(kk,:) = xyzind(sortind(kk), :);
%         end
%         
%         vlist = unique(vlist, 'rows');
%         
%         pproj = zeros(size(vlist, 1), 3);
%         for kk = 1:size(vlist, 1)
%             pproj(kk, :) = vlist(kk, :) - dot(vlist(kk, :) - cv(:, qq).', univec(:,qq)) * univec(:,qq).';
%         end
%         
%         while (max(pproj(:,3)) - min(pproj(:,3)) > 1)
%             w = [0 0 0];
%             while all(w == [0 0 0])
%                 rng1 = randi(size(pproj,1));
%                 w = cross(pproj(randi(size(pproj,1)), :) - pproj(rng1, :), pproj(randi(size(pproj,1)), :) - pproj(rng1, :));
%             end
%             
%             w = w/norm(w);
%             R = [null(w), w.'];
%             
%             if det(R) < 0
%                 R(:, 1:2) = R(:, 2:-1:1);
%             end
%             pproj = pproj*R;
%         end
%         
%         xuunit = [1 0 0]*R^(-1);
%         yvunit = [0 1 0]*R^(-1);
%         uunit = ((xuunit(1)*xscale)^2+(xuunit(2)*yscale)^2+(xuunit(3)*zscale)^2)^(1/2);
%         vunit = ((yvunit(1)*xscale)^2+(yvunit(2)*yscale)^2+(yvunit(3)*zscale)^2)^(1/2);
%         
%         
%         csimage = zeros(ceil(max(pproj(:,1)) - min(pproj(:,1))) + 10, ceil(max(pproj(:,2)) - min(pproj(:,2))) + 10);
%         
%         for p = 1:size(pproj,1)
%             csimage(round(pproj(p,1) - min(pproj(:,1))) + 6, round(pproj(p,2) - min(pproj(:,2))) + 6) = 1;
%         end
%         
%         [blobImage, ~] = bwlabel(csimage);
%         blobMeasurements = regionprops(blobImage, 'area');
%         allAreas = [blobMeasurements.Area];
%         
%         [~, sortIndexes] = sort(allAreas, 'descend');
%         biggestBlob = ismember(blobImage, sortIndexes(1));
%         csobj = biggestBlob > 0;
%         
%         csobj = imfill(csobj, 'holes');
%         area = regionprops(csobj, 'Area','perimeter');
%         segcs(qq,i) = area.Area * uunit * vunit;
%         segper(qq,i) = area.Perimeter;
%         
% %         % num objects per cross-section
% %         lm14 = bwconncomp(allFiles.(userDefinedallFilesName).bwcellimage);
% %         bwcellimagepostfilter = false(size(lm14.ImageSize));
% %         for ii = lm98(:,1)
% %         bwcellimagepostfilter(lm14.PixelIdxList{ii}) = true;
% %         end
% %         h5bwdatabox = bwcellimagepostfilter(floor(boxloc(3)):ceil(boxloc(4)), floor(boxloc(1)): ceil(boxloc(2)), floor(boxloc(5)):ceil(boxloc(6)));
% %         
% %         [xind, yind, zind] = ind2sub(size(h5bwdatabox), find(h5bwdatabox == 1));
% %         xyzind = [yind + floor(boxloc(1)) - 1, xind + floor(boxloc(3)) - 1, zind + floor(boxloc(5)) - 1];
% %         
% %         vdist2plane = zeros(size(xyzind, 1), 1);
% %         
% %         for kk = 1:size(xyzind, 1)
% %             vdist2plane(kk) = dot(univec(:,qq), xyzind(kk,:)' - cv(:,qq));
% %         end
% %         
% %         v = abs(vdist2plane);
% %         sortind = find(v(:) <= 1);
% %         
% %         vlist = zeros(size(sortind,1),3);
% %         
% %         for kk = 1:size(sortind, 1)
% %             vlist(kk,:) = xyzind(sortind(kk), :);
% %         end
% %         
% %         vlist = unique(vlist, 'rows');
% %         
% %         pproj = zeros(size(vlist, 1), 3);
% %         for kk = 1:size(vlist, 1)
% %             pproj(kk, :) = vlist(kk, :) - dot(vlist(kk, :) - cv(:, qq).', univec(:,qq)) * univec(:,qq).';
% %         end
% %         
% % %         while (max(pproj(:,3)) - min(pproj(:,3)) > 1)
% % %             w = [0 0 0];
% % %             while all(w == [0 0 0])
% % %                 rng1 = randi(size(pproj,1));
% % %                 w = cross(pproj(randi(size(pproj,1)), :) - pproj(rng1, :), pproj(randi(size(pproj,1)), :) - pproj(rng1, :));
% % %             end
% % %             
% % %             w = w/norm(w);
% % %             R = [null(w), w.'];
% % %             
% % %             if det(R) < 0
% % %                 R(:, 1:2) = R(:, 2:-1:1);
% % %             end
% % %             pproj = pproj*R;
% % %         end
% % %         
% % %         xuunit = [1 0 0]*R^(-1);
% % %         yvunit = [0 1 0]*R^(-1);
% % %         uunit = ((xuunit(1)*xscale)^2+(xuunit(2)*yscale)^2+(xuunit(3)*zscale)^2)^(1/2);
% % %         vunit = ((yvunit(1)*xscale)^2+(yvunit(2)*yscale)^2+(yvunit(3)*zscale)^2)^(1/2);
% %         
% %         
% %         csimage = zeros(ceil(max(pproj(:,1)) - min(pproj(:,1))) + 10, ceil(max(pproj(:,2)) - min(pproj(:,2))) + 10);
% %         
% %         for p = 1:size(pproj,1)
% %             csimage(round(pproj(p,1) - min(pproj(:,1))) + 6, round(pproj(p,2) - min(pproj(:,2))) + 6) = 1;
% %         end
% %         
% %         result = bwconncomp(csimage);
% %         csperseg(qq,i) = result.NumObjects;
%         
%     end
% end


%% save


% fprintf('Save?\n');
% preview = input('','s');
% 
% if preview == 'Y'
    
fprintf('Saving all data...\n');

allFiles.(userDefinedallFilesName).branchid = branchid;
allFiles.(userDefinedallFilesName).branchposition = branchposition;
allFiles.(userDefinedallFilesName).celldata = lm98;
allFiles.(userDefinedallFilesName).csdata = csdata;


% allFiles.(userDefinedallFilesName).seglen = seglen;
% allFiles.(userDefinedallFilesName).segcs = segcs;
% allFiles.(userDefinedallFilesName).segper = segper;
% allFiles.(userDefinedallFilesName).csperseg = csperseg;
save(userDefinedallFilesName,'allFiles','-append');


fprintf('Done. Updated allFiles in %s\n',userDefinedallFilesName);

% end
