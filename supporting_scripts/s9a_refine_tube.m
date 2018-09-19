% function s9a_refine_tube(varargin)

% narginchk(0,1)

%% voronoi diagram v2


% import .h5 pixel class unrolled tube
% if strcmp(userDefinedallFilesName,FileName(1:length(userDefinedallFilesName)))
%     fprintf('Using previously imported data as:\n    %s...\n',userDefinedallFilesName);
% else
%     fprintf('Import the .h5 pixel class data for %s.\n',userDefinedallFilesName);
%
%     [FileName,PathName] = uigetfile('*.h5', 'Select the .h5 file.', 'MultiSelect', 'off');
%     ilastik_data = h5read(fullfile(PathName,FileName),'/exported_data');
%     if ~strcmp(userDefinedallFilesName,FileName(1:length(userDefinedallFilesName)))
%         fprintf('Warning: Current file name is:\n    %s\nbut .h5 file appears to be for a different dataset:\n    %s\nContinue anyway? (Y/N)\n',userDefinedallFilesName,FileName);
%         preview = input('','s');
%         if preview ~= 'Y'
%             return
%         end
%     else
%         fprintf('%s imported successfully.\n',FileName);
%     end

fprintf('Reading %s...\n',[userDefinedallFilesName,'unrobj.h5']);

% ilastik_data = h5read([userDefinedallFilesName,'unrobj.h5'],'/exported_data');
% 
% if ~ismatrix(ilastik_data)
%     ilastik_data = squeeze(ilastik_data);
%     if ~ismatrix(ilastik_data)
%         error('Data dims > 2');
%     else
%         warning('Extra dims exported from Ilastik.');
%     end
% end
% % basic export param check
% if size(ilastik_data,2) < size(ilastik_data,1)
%     ilastik_data = permute(ilastik_data,[2 1]);
% end

ilastik_data = squeeze(h5read([userDefinedallFilesName,'unrobj.h5'],'/exported_data'));


    if ~ismatrix(ilastik_data)
        error('Data dims > 2');
    else
        warning('Extra dims exported from Ilastik.');
    end

% basic export param check
if size(ilastik_data,2) < size(ilastik_data,1)
    ilastik_data = permute(ilastik_data,[2 1]);
end


%%

% required variables: skelnums, S
S = allFiles.(userDefinedallFilesName).S;
skelnums = allFiles.(userDefinedallFilesName).skelnums;
load([userDefinedallFilesName,'unrdata'],'all_indices');
branchposition = allFiles.(userDefinedallFilesName).branchposition;
branchid = allFiles.(userDefinedallFilesName).branchid;
segcs = allFiles.(userDefinedallFilesName).segcs;



%%

% ilastik_data_exp = imdilate(ilastik_data,strel('disk',5));

ilastik_regionprops = regionprops(bwconncomp(ilastik_data));
% ilastik_regionprops = regionprops(ilastik_data);
ilastik_centroids = zeros(size(ilastik_regionprops,1),2);
for i = 1:size(ilastik_regionprops,1)
    ilastik_centroids(i,:) = ilastik_regionprops(i).Centroid;
end

ilastik_centroids(ilastik_centroids(:,1)<branchposition(1) | ilastik_centroids(:,1)>branchposition(end),:)=[];

% end



%%




% generate voronoi

fprintf('Generating Voronoi diagram...\n');
% [ilastik_v, ilastik_c] = voronoin(ilastik_centroids);
% plot initial voronoi diagram
ilastik_xlim = [branchposition(1) branchposition(end)];
ilastik_ylim = [0 size(ilastik_data,1)];
ilastik_ext = [ilastik_xlim(1) ilastik_xlim(2) ilastik_xlim(2) ilastik_xlim(1) ilastik_xlim(1); ...
    ilastik_ylim(2) ilastik_ylim(2) ilastik_ylim(1) ilastik_ylim(1) ilastik_ylim(2)]';
[ilastik_v, ilastik_c] = VoronoiLimit(ilastik_centroids(:,1),ilastik_centroids(:,2),'bs_ext',ilastik_ext,'figure','off');
% plot
% if nargin == 1
    figure();
    axis equal;hold on;
    set(gca,'FontSize',24);
    title('Voronoi Diagram','FontSize',32);
    axis([branchposition(1)-5 branchposition(end)+5 -5 size(ilastik_data,1)+5]);
    xlabel('Relative A/P position (pix)','FontSize',24);
    ylabel('Circumferential position (deg/4)','FontSize',24);

    if ilastik_v(1) == Inf
        for i = 1:length(ilastik_c)
            if all(ilastik_c{i}~=1)
                patch(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2),i); % use color i.
            end
        end
    else
        for i = 1:length(ilastik_c)
            patch(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2),i); % use color i.
        end
    end
    colorbar('ticks',[1 length(ilastik_c)],'ticklabels',{'Anterior','Posterior'},'location','southoutside');
% end



% calculate polygon areas
fprintf('Calculating polygon areas...\n');
ilastik_polyareas = nan(size(ilastik_c));
if ilastik_v(1) == Inf
    for i = 1:length(ilastik_c)
        if all(ilastik_c{i}~=1)
            ilastik_polyareas(i) = polyarea(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2));
        end
    end
else
    for i = 1:length(ilastik_c)
        ilastik_polyareas(i) = polyarea(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2));
    end
end
% if nargin == 1
    % plot
    figure();
    axis equal;hold on;
    set(gca,'FontSize',24);
    title('Voronoi Diagram (Cell Areas)','FontSize',32);
    axis([branchposition(1)-5 branchposition(end)+5 -5 size(ilastik_data,1)+5]);
    xlabel('Relative A/P position (pix)','FontSize',24);
    ylabel('Circumferential position (deg/4)','FontSize',24);
    
    if ilastik_v(1) == Inf
        for i = 1:length(ilastik_c)
            if all(ilastik_c{i}~=1)
                patch(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2),ilastik_polyareas(i)); % use color i.
            end
        end
    else
        for i = 1:length(ilastik_c)
            patch(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2),ilastik_polyareas(i)); % use color i.
        end
    end
    colorbar;
% end



% % remove nans for data analysis
% ilastik_y = ilastik_polyareas; % temp data storage
% ilastik_x = ilastik_centroids(:,1);
%
% ilastik_isoutlier_inds = isoutlier(ilastik_y); % keep track of excluded datapoints
%
% ilastik_y(ilastik_isoutlier_inds) = NaN; % define NaNs as outliers and remove them
% ilastik_x(isnan(ilastik_y)) = [];
% ilastik_y(isnan(ilastik_y)) = [];
% % % check
% % figure();
% % scatter(ilastik_x,ilastik_y);



% filter (size, orientation?)
fprintf('Calculating cells per cross-section...\n');


% calc cells per cross (this version which uses unique objects id'd from
% ilastik only requires the object classification and is much faster than
% the VoronoiLimit method but doesn't take empty space into account
% i.e. junctions are not infinitely thin lines.)

% % method 1: window width of 5
%
% ilastik_cpc = zeros(length(ilastik_data),1);
% for i = 3:length(ilastik_data)-2
%     ilastik_cpc(i) = length(unique(ilastik_data(round(size(ilastik_data,1)/3):round(size(ilastik_data,1)*2/3),i-2:i+2)));
% end
% figure();
% plot(ilastik_cpc(2:end-2));

% % method 2: window width of 1
%
% ilastik_cpc = zeros(length(ilastik_data),1);
% for i = 1:length(ilastik_data)
%     ilastik_cpc(i) = length(unique(ilastik_data(round(size(ilastik_data,1)/3):round(size(ilastik_data,1)*2/3),i)));
% end
% figure();hold on;
% scatter(1:length(ilastik_cpc),ilastik_cpc);
% y = sgolayfilt(ilastik_cpc/3,5,299);
% plot(y,'.-');

% method 3: object expansion before cpc

ilastik_data_exp = imdilate(bwlabel(ilastik_data),strel('disk',5));
% ilastik_data_unique = unique(ilastik_data_exp(:));
% ilastik_data_unique(1) = [];
% for i = 1:length(ilastik_data_unique)
%     ilastik_data_exp(ilastik_data_exp==ilastik_data_unique(i))=i;
% end
ilastik_cpc = zeros(length(ilastik_data_exp),1);
for i = 1:length(ilastik_data_exp)
    ilastik_cpc(i) = sum(unique(ilastik_data_exp(:,i))~=0)/3;
end

% if nargin == 1
    % plot
    figure();
    hold on;
    set(gca,'FontSize',24);
    title('Cells per cross-section','FontSize',32);
    axis([branchposition(1)-5 branchposition(end)+5 0 inf]);
    xlabel('Relative A/P position (pix)','FontSize',24);
    ylabel('Number of cells','FontSize',24);
    
    scatter(1:length(ilastik_cpc),ilastik_cpc);
    y = sgolayfilt(ilastik_cpc,5,299);
    plot(y,'.-');
    
    for i = 1:length(branchposition)
        plot([branchposition(i) branchposition(i)],[0 10],'r-');
%         plot([branchposition(i+1) branchposition(i+1)],[0 10],'g-');
    end
    
    % figure();
    % plot(ilastik_cpc(2:end-2));
% end



% calculate regionprops 2d
ilastik_n_position = zeros(branchposition(end)-branchposition(1)+1,1);
branchposition = branchposition - branchposition(1) + 1;
ilastik_n_position(1) = branchid(1);
for i = 1:length(branchid)-1
    ilastik_n_position(branchposition(i)+1:branchposition(i+1)) = (branchid(i)+(branchid(i+1)-branchid(i))/(branchposition(i+1)-branchposition(i)):(branchid(i+1)-branchid(i))/(branchposition(i+1)-branchposition(i)):branchid(i+1))';
end
branchposition = allFiles.(userDefinedallFilesName).branchposition;

ilastik_centroids(:,3) = ilastik_n_position(round(ilastik_centroids(:,1))-branchposition(1));

% normalization data
xx = [];
for i = 1:length(branchid)-1
    xx = [xx,branchid(i):(branchid(i+1)-branchid(i))/10:branchid(i+1)];
    xx(end) = [];
end
[pix_a,pix_b,~] = linreg(xx,segcs(:));
clear xx


ilastik_rp = regionprops(ilastik_data_exp,'Image','Eccentricity','Orientation','Centroid','MajorAxisLength','MinorAxisLength');
ilastik_rp_2d = cell(length(ilastik_rp)+1,7);
ilastik_rp_2d{1,1} = 'Cell ID';
ilastik_rp_2d{1,2} = 'Centroid*';
ilastik_rp_2d{1,3} = 'Eccentricity*';
ilastik_rp_2d{1,4} = 'Aspect Ratio*';
ilastik_rp_2d{1,5} = 'Rel. size (pix)';
ilastik_rp_2d{1,6} = 'Orientation (degrees)';
ilastik_rp_2d{1,7} = 'Normalized Position';
for i = 1:length(ilastik_rp)
    if ~isnan(ilastik_rp(i).Centroid(2))
        try
            ilastik_rp_2d{i+1,7} = ilastik_n_position(round(ilastik_rp(i).Centroid(1))-branchposition(1));
        catch
            ilastik_rp_2d{i+1,7} = NaN;
        end
    end
    ilastik_y_pix_scale = pix_a*ilastik_rp_2d{i+1,7}+pix_b;
    ilastik_rp_2d{i+1,1} = i;
    ilastik_rp_2d{i+1,2} = ilastik_rp(i).Centroid;
    ilastik_rp_2d{i+1,3} = ilastik_rp(i).Eccentricity;
    ilastik_rp_2d{i+1,4} = ilastik_rp(i).MajorAxisLength/ilastik_rp(i).MinorAxisLength;
    ilastik_rp_2d{i+1,5} = bwarea(ilastik_rp(i).Image)*ilastik_y_pix_scale;
    ilastik_rp_2d{i+1,6} = atand(sind(ilastik_rp(i).Orientation)*ilastik_y_pix_scale/(cosd(ilastik_rp(i).Orientation)*ilastik_y_pix_scale));
    
end
ilastik_rp_2d = ilastik_rp_2d(all(cellfun(@(x)any(~isnan(x)),ilastik_rp_2d),2),:); % remove rows that have nan in them


% % calc cells per cross (this version which uses VoronoiLimit to predefine
% % the unrolled image area to strips (width 3) before asking number of cell
% % objects takes a really long time to compute).
%
% ilastik_cpc = zeros(length(ilastik_data),1);
% for i = 2:5:length(ilastik_data)-1
%     fprintf('Iteration: %i of %i\n',i-1,length(ilastik_data)-2);
%     ilastik_xlim = [i-1, i+1];
%     ilastik_ylim = [round(size(ilastik_data,1)/3) round(size(ilastik_data,1)*2/3)];
%     ilastik_ext = [ilastik_xlim(1) ilastik_xlim(2) ilastik_xlim(2) ilastik_xlim(1) ilastik_xlim(1); ...
%                    ilastik_ylim(2) ilastik_ylim(2) ilastik_ylim(1) ilastik_ylim(1) ilastik_ylim(2)]';
%     ilastik_int = [];
%     [~, ilastik_vl_c] = VoronoiLimit(ilastik_centroids(:,1),ilastik_centroids(:,2),'bs_ext',ilastik_ext,'bs_int',ilastik_int,'figure','off');
%     ilastik_cpc(i) = length(ilastik_vl_c);
% end
%
% figure();
% plot(ilastik_cpc(2:end-1))





% calc cell connectivity
fprintf('Calculating object connectivity...\n');

ilastik_cc = zeros(length(ilastik_centroids),1);
for i = 1:length(ilastik_cc)
    ilastik_cc(i) = max([4,min([8,length(ilastik_c{i})])]);
end


% if nargin == 1
    % plot voronoi diagram
    figure();
    axis equal;hold on;
    set(gca,'FontSize',24);
    title('Voronoi Diagram (Cell Connectivity)','FontSize',32);
    axis([branchposition(1)-5 branchposition(end)+5 -5 size(ilastik_data,1)+5]);
    xlabel('Relative A/P position (pix)','FontSize',24);
    ylabel('Circumferential position (deg/4)','FontSize',24);
    
    for i = 1:length(ilastik_c)
        patch(ilastik_v(ilastik_c{i},1),ilastik_v(ilastik_c{i},2),ilastik_cc(i)); % use color i.
    end
    colorbar;
    
%     hold on;
%     for i = 1:length(ilastik_unrolled_node_position)
%     t = linspace(0,2*pi);
%     r = 5;
%     x = ilastik_unrolled_node_position(i,1)+r*cos(t);
%     y = ilastik_unrolled_node_position(i,2)+r*sin(t)+90;
%     patch(x,y,[1 1 1]);
%     end


% end
% % plot connectivity (edge number = vertex number) based on axial position
% figure();
% scatter(ilastik_centroids(:,1),ilastik_cc);




% calc positions of branch points
fprintf('Calculating branch positions...\n');

% include holes?


% find indices of branch points
% stores co-ords of beginning (:,:,1) and ending (:,:,2) of branch segments
% removes co-ords of primary skeleton segments (stores as NaN)
ilastik_S_inds = zeros(length(S),3,2);
for i = 1:length(S)
    ilastik_S_inds(i,:,1) = S{i}(1,:);
    ilastik_S_inds(i,:,2) = S{i}(end,:);
end
ilastik_S_inds(skelnums,:,:) = NaN;

% finds indices branch points id'd by skelnums
% i.e. this includes nodes for TCs, DBs, and any other stray nodes in
% between. need to filter again in the next for loop
ilastik_S_index_positions = zeros(length(skelnums),1);
ilastik_S_dist = 0;
for i = 1:length(skelnums)
    ilastik_S_dist = length(S{skelnums(i)}) + ilastik_S_dist - 1;
    ilastik_S_index_positions(i) = ilastik_S_dist;
end

% finds indices of branch points id'd by branchid
% i.e. this only uses the locations of TCs that we manually input using
% branchid
ilastik_S_node_positions = nan(length(branchposition),1);
for i = 1:length(branchposition)
    if branchid(i) - round(branchid(i)) < 10^-5 % this ignores all non-branch indices (non-whole numbers in branchid)
        j = find(ismember(ilastik_S_index_positions, branchposition(i)));
        if branchposition(i) == 1
            ilastik_S_node_positions(i) = NaN; % the start of DT will ignored for now bc we usually don't include branches there
        elseif ~isempty(j)
            ilastik_S_node_positions(i) = skelnums(j);
        end
    end
end

% find circumferential vector that closest matches each TC branch vector,
% i.e. calculate y position of hole in unrolled tube image
ilastik_unrolled_node_position = zeros(length(ilastik_S_node_positions),2);
for i = 1:length(ilastik_S_node_positions)
    % find branch id# that corresponds to TC node location
    if isnan(ilastik_S_node_positions(i))
        continue
    end
    ilastik_intersect_ind = find(ismember(squeeze(ilastik_S_inds(:,:,2)),S{ilastik_S_node_positions(i)}(end,:),'rows'));
    if isempty(ilastik_intersect_ind)
        warning('Failed to find branch point at the end of node %i. Skipping...',ilastik_S_node_positions(i));
        continue
    end
    if numel(ilastik_intersect_ind) >= 2 % accounts for branches occurring before DT ROI
        continue
    end
    % determine branch vector (closest to TC node)
    ilastik_intersect_v = S{ilastik_intersect_ind}(end-1,:) - S{ilastik_intersect_ind}(end,:);
    % find array of circumferential vectors
    ilastik_intersect_array = squeeze(all_indices(branchposition(i),:,:,:));
    ilastik_cross_array = nan(size(ilastik_intersect_array,1), 3);
    for j = 1:size(ilastik_intersect_array,1)
        % ilastik_cross_array(j,:) = (ilastik_intersect_array(j,:,2) - ilastik_intersect_array(j,:,1))';
        ilastik_cross_array(j,:) = (ilastik_intersect_array(j,:,min([size(ilastik_intersect_array,3) 10])) - ilastik_intersect_array(j,:,1))';
    end
    % calculate direction of TC vector w.r.t. normalized circumferential angle
    ilastik_intcross_v_delta = nan(size(ilastik_intersect_array,1),1);
    for j = 1:size(ilastik_intersect_array,1)
        ilastik_intcross_v_delta(j) = dot(ilastik_intersect_v,ilastik_cross_array(j,:));
    end
    % store unrolled TC location data
    ilastik_unrolled_node_position(i,1) = ilastik_S_index_positions(skelnums==ilastik_S_node_positions(i));
    % ilastik_unrolled_node_position(i,1) = branchposition(i); % backup code incase above doesn't work
    ilastik_unrolled_node_position(i,2) = median(find(ismember(ilastik_intcross_v_delta,max(ilastik_intcross_v_delta))));
end
ilastik_unrolled_node_position(ilastik_unrolled_node_position(:,1)==branchposition(end),:) = [];

% % create holes in unrolled image
fprintf('Generating branch points as holes in unrolled image...\n');
% unr_max_int_values_t = false(size(unr_max_int_values));
% unr_max_int_values_w_holes = unr_max_int_values;
% for i = find(any(ilastik_unrolled_node_position,2))
%     unr_max_int_values_t(round(ilastik_unrolled_node_position(i,2)),round(ilastik_unrolled_node_position(i,1))) = true;
% end
% unr_max_int_values_t = imdilate(unr_max_int_values_t,strel('disk',7));
% unr_max_int_values_w_holes(unr_max_int_values_t) = 0;
%
% figure();imshow(unr_max_int_values_w_holes); % show plot with holes
% hold on;
% plot(ilastik_unrolled_node_position(:,1),ilastik_unrolled_node_position(:,2),'r*'); % plot nodes with a red asterisk
% figure();imshow(unr_max_int_values); % show plot with no holes
% hold on;
% plot(ilastik_unrolled_node_position(:,1),ilastik_unrolled_node_position(:,2),'r*'); % plot nodes with a red asterisk



% plot voronoi with branch holes

ilastik_xlim = [branchposition(1) branchposition(end)];
ilastik_ylim = [0 round(size(ilastik_data,1)/3)]; % round(size(ilastik_data,1)*2/3)];
ilastik_ext = [ilastik_xlim(1) ilastik_xlim(2) ilastik_xlim(2) ilastik_xlim(1) ilastik_xlim(1); ...
    ilastik_ylim(2) ilastik_ylim(2) ilastik_ylim(1) ilastik_ylim(1) ilastik_ylim(2)]';
ilastik_int = cell(length(ilastik_unrolled_node_position),1);
ilastik_t = linspace(0,2*pi)';
ilastik_tube_r = 11;
for i = find(~ismember(ilastik_unrolled_node_position,[0 0],'rows'))'
    ilastik_int{i} = [(cos(ilastik_t)*ilastik_tube_r)+ilastik_unrolled_node_position(i,1) (sin(ilastik_t)*ilastik_tube_r)+ilastik_unrolled_node_position(i,2)];
end
ilastik_int = ilastik_int(~cellfun('isempty',ilastik_int));
[ilastik_vl_v, ilastik_vl_c] = VoronoiLimit(ilastik_centroids(:,1),ilastik_centroids(:,2),'bs_ext',ilastik_ext,'bs_int',ilastik_int,'figure','off');

% figure();axis equal;axis([-5 size(ilastik_data,2)+5 -5 size(ilastik_data,1)+5]);hold on;

% if nargin == 1
    % plot
    figure();
    axis equal;hold on;
    set(gca,'FontSize',24);
    title('Voronoi Diagram (with transverse branch points as holes)','FontSize',32);
    axis([branchposition(1)-5 branchposition(end)+5 -5 round(size(ilastik_data,1)/3)+5]);
    xlabel('Relative A/P position (pix)','FontSize',24);
    ylabel('Circumferential position (deg/4)','FontSize',24);
    
    if ilastik_vl_v(1) == Inf
        for i = 1:length(ilastik_vl_c)
            if all(ilastik_vl_c{i}~=1)
                patch(ilastik_vl_v(ilastik_vl_c{i},1),ilastik_vl_v(ilastik_vl_c{i},2),i); % use color i.
            end
        end
    else
        for i = 1:length(ilastik_vl_c)
            patch(ilastik_vl_v(ilastik_vl_c{i},1),ilastik_vl_v(ilastik_vl_c{i},2),i); % use color i.
        end
    end
    colorbar('ticks',[1 length(ilastik_vl_c)],'ticklabels',{'Anterior','Posterior'},'location','southoutside');
% end


return

%% save

save([userDefinedallFilesName,'unrdata'],'ilastik_polyareas','ilastik_data','ilastik_rp_2d','ilastik_cpc','ilastik_cc','ilastik_n_position','ilastik_centroids','-append')
fprintf('Data saved to %s.\n',[userDefinedallFilesName,'unrdata']);

% if strcmp(userDefinedallFilesName,FileName(1:length(userDefinedallFilesName)))
%     save([userDefinedallFilesName,'unrdata'],'ilastik_polyareas','ilastik_data','ilastik_rp_2d','ilastik_cpc','ilastik_cc','ilastik_n_position','ilastik_centroids','-append')
%     fprintf('Data saved to %s.\n',[userDefinedallFilesName,'unrdata']);
% else
%     warning('Data not saved due to allFiles and .h5 name mismatch:\n    %s\n    %s\n',userDefinedallFilesName,FileName);
% end

return





%%

% filter w.r.t. x
ilastik_xlim = [400 600];
ilastik_ylim = [round(size(ilastik_data,1)/3) round(size(ilastik_data,1)*2/3)];
[ilastik_vl_v, ilastik_vl_c] = VoronoiLimit(ilastik_centroids(:,1),ilastik_centroids(:,2),'bs_ext', ...
    [ilastik_xlim(1) ilastik_xlim(2) ilastik_xlim(2) ilastik_xlim(1) ilastik_xlim(1); ...
    ilastik_ylim(2) ilastik_ylim(2) ilastik_ylim(1) ilastik_ylim(1) ilastik_ylim(2)]', ...
    'figure','off');

figure();axis equal;axis([ilastik_xlim(1)-5 ilastik_xlim(2)+5 ilastik_ylim(1)-5 ilastik_ylim(2)+5]);hold on;
if ilastik_vl_v(1) == Inf
    for i = 1:length(ilastik_vl_c)
        if all(ilastik_vl_c{i}~=1)
            patch(ilastik_vl_v(ilastik_vl_c{i},1),ilastik_vl_v(ilastik_vl_c{i},2),i); % use color i.
        end
    end
else
    for i = 1:length(ilastik_vl_c)
        patch(ilastik_vl_v(ilastik_vl_c{i},1),ilastik_vl_v(ilastik_vl_c{i},2),i); % use color i.
    end
end


%% use same search indices to replot data
% surface data should be under bwcellimage.
% need to load scurve and bwcellimage from allFiles

scurve = allFiles.(userDefinedallFilesName).scurve;
bwcellimage = allFiles.(userDefinedallFilesName).bwcellimage; % use this to choose the bwcellimage file you need to unroll

values = fnval(scurve, 1:scurve.pieces);
der = fnder(scurve);
tangent_vectors = fnval(der, 1:scurve.pieces);

% base_ref_vector = [0;0;0];
% [~,direction] = min(size(bwcellimage));
% base_ref_vector(direction) = 1;% set this

unr_ave_int_values = zeros(ceil(360/angle_delta)+1,floor(scurve.pieces/x_delta)); % allocate size
unr_max_int_values = zeros(ceil(360/angle_delta)+1,floor(scurve.pieces/x_delta));

for w = 1:floor(scurve.pieces/x_delta)
    W = w*x_delta;
    fprintf('Progress: %i of %i\n',W,scurve.pieces);
    
    int_averages = zeros(ceil(360/angle_delta)+1,1); % allocate size
    int_maxes = zeros(ceil(360/angle_delta)+1,1);
    
    
    for angle_iteration = 0:360/angle_delta % num of iterations?
        search_indices = squeeze(all_indices(w,angle_iteration+1,:,:));
        
        search_indices(:,search_indices(3,:)>size(bwcellimage,3))=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(2,:)>size(bwcellimage,2))=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(1,:)>size(bwcellimage,1))=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(3,:)<1)=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(2,:)<1)=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(1,:)<1)=[]; % remove out-of-bounds coords
        % if size(search_indices,2) ~= 0
        try
            ind = sub2ind(size(bwcellimage),search_indices(1,:),search_indices(2,:),search_indices(3,:));
            int_averages(angle_iteration+1) = mean(bwcellimage(ind(~isnan(ind)))); % pull int values and average them
            int_maxes(angle_iteration+1) = max(bwcellimage(ind(~isnan(ind)))); % or pull max int values
        % else
        catch
            % need to make an exception for inds that are outside of the
            % image, usually this means that the ROI is hitting the edge of
            % the image (need to take a bigger stack).
            warning('Radial indices are out of bounds at position [%i,%i] of %i. Skipping...',W,angle_iteration*angle_delta,scurve.pieces);
        end
    end
    
    unr_ave_int_values(:,W) = int_averages;
    unr_max_int_values(:,W) = int_maxes;
    
    % base_ref_vector = proj_ref_vector; % use this to dynamically change the ref vector on each cross section
end

% unrolled_plot = mat2gray(unr_max_int_values,[0,max(unr_max_int_values(:))]);
figure();imshow(unr_max_int_values);

% to save as userDefinedallFilesName for 2d segmentation in ilastik, use the generic line:
imwrite(repmat(unr_max_int_values,3,1),[userDefinedallFilesName,'.jpg'],'jpg');


%% 171120 - comments on pixel class
% to fix the "hashes", we can just use ilastik pixel class. the initial
% training is a little time consuming but afterward it should be easy to
% batch process and export that data.

% steps:
% 1. save unr_max_int_values as a tiff file using
%imwrite(unr_max_int_values,'test.tiff','tiff');
% this saves it in the matlab folder as test.tiff.
% 2. run the pixel class module in ilastik and teach it to recognize cells
% from non-cells. export the binary segmentation result as an h5 file.
% 3. use
%testimage = h5read('test.h5','/exported_data');
% to import back into matlab and then
%im2bw(squeeze(testimage(1,:,:)),.6);
% and change the .6 to the threshold that makes the most cells.
% 4. there's an ilastik save file with some teaching params already set on
% the desktop which means we can skip that step in future runs.


%% voronoi diagram
a = unr_max_int_values;
d = imdilate(imerode(a,strel('disk',1,4)),strel('disk',1,4)); % smoothing?)
c = repmat(d,3,1);
e = regionprops(bwconncomp(c));
f = zeros(size(e,1),2);
for i = 1:size(e,1)
    f(i,:) = e(i).Centroid;
end
voronoi(f(:,1),f(:,2));axis equal;


%% voronoi diagram with colors based on edge number

a = unr_max_int_values;
d = imdilate(imerode(a,strel('disk',1,4)),strel('disk',1,4)); % smoothing?)
c = repmat(d,5,1);
e = regionprops(bwconncomp(c));
f = zeros(size(e,1),2);
for i = 1:size(e,1)
    f(i,:) = e(i).Centroid;
end
[v,c] = voronoin(f);
% color = [[1 0 0];[226/255, 87/255, 30/255];[1 .5 0];[1 1 0];[0 1 0];[150/255 191/255 51/255];[0 0 1];[75/255 0 130/255];[139/255 0 1]];
figure();axis equal;axis([0 size(unr_max_int_values,2) size(unr_max_int_values,1)-5 3*size(unr_max_int_values,1)+5]);hold on;
for i = 1:length(c)
    if all(c{i}~=1)   % If at least one of the indices is 1,
        % then it is an open region and we can't
        % patch that.
        patch(v(c{i},1),v(c{i},2),i); % use color i.
    end
end


%%

a = unr_max_int_values;
%a(:,sum(a)==0)=[]; % clear empty columns
% a(sum(a,2)==0,:)=[]; % clear empty rows
b = mat2gray(a,[0,max(a(:))]);
%figure();imshow(b);
c = repmat(b,3,1);
%figure();imshow(c);
d = imdilate(imerode(c,strel('disk',1,4)),strel('disk',1,4)); % smoothing?
figure();imshow(d);

%% generate watershed

e = watershed(imcomplement(d));
f = imcomplement(logical(e));
figure();imshow(f);