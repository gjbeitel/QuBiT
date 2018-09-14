%% pull data

scurve = allFiles.(userDefinedallFilesName).scurve;
h5bwfilled = allFiles.(userDefinedallFilesName).h5bwfilled;
% image source
bwcellimage = bwlabeln(allFiles.(userDefinedallFilesName).bwcellimage);
load([userDefinedallFilesName,'unrdata'],'all_indices');

% need to select image
% 
% if exist('lm12','var')==1
%     fprintf('Cell image found.\n');
% else
%     fprintf('\nSelect the cell (objects) .h5 file.\n');
%     FileName = uigetfile('*.h5', 'Select the cell (objects) .h5 file.', 'MultiSelect', 'off');
%     try
%         lm12 = h5read(FileName,'/exported_data');
%     catch
%         error('Invalid file or unsupported data type.');
%     end
%     if numel(size(squeeze(lm12))) > 3
%         error('Objects .h5 file has more than one data channel.\nBe sure that "Objects" was exported correctly from the Pixel+Object Classification module in Ilastik.');
%     end
% end
% image = squeeze(lm12);



%% unroll tube

if ~exist('all_indices','var')
values = fnval(scurve, 1:scurve.pieces);
der = fnder(scurve);
tangent_vectors = fnval(der, 1:scurve.pieces);

base_ref_vector = [0;0;0];
[~,direction] = min(size(bwcellimage));
base_ref_vector(direction) = 1;% set this

syms scan_x scan_y scan_z

% vars:
radial_search_dist = 27; % max pixel dist of how far to search from the center
angle_delta = 4; % search angle per cross-section
x_delta = 1;
%

% unr_ave_int_values = zeros(ceil(360/angle_delta)+1,floor(scurve.pieces/x_delta)); % allocate size
unr_image = zeros(ceil(360/angle_delta)+1,floor(scurve.pieces/x_delta));
all_indices = nan(floor(scurve.pieces/x_delta),floor(360/angle_delta),3,40);
% time = zeros(1,scurve.pieces);

for w = 1:floor(scurve.pieces/x_delta) % call from values and tangent_vectors

    W = w*x_delta;
    
    fprintf('Progress: %i of %i\n',W,scurve.pieces);
%    tic
    
    spline_location = values(:,W);
    spline_tangent = tangent_vectors(:,W);
    
    unit_norm = spline_tangent/norm(spline_tangent);
    proj_ref_vector = base_ref_vector - dot(base_ref_vector, unit_norm)*unit_norm;
    unit_proj = proj_ref_vector/norm(proj_ref_vector);
    
    int_averages = zeros(ceil(360/angle_delta)+1,1); % allocate size
    ints = zeros(ceil(360/angle_delta)+1,1);
    
    
    for angle_iteration = 0:360/angle_delta % num of iterations?
        % fprintf('  Subprogress: %i of 72\n',angle_iteration);
        
        angle = angle_iteration * angle_delta;
        
        %eqn1 = unit_proj(2)*scan_z - unit_proj(3)*scan_y == unit_norm(1)*sind(angle);
        %eqn2 = unit_proj(1)*scan_z - unit_proj(3)*scan_x == unit_norm(2)*sind(angle);
        eqn3 = unit_proj(1)*scan_y - unit_proj(2)*scan_x == unit_norm(3)*sind(angle);
        
        eqn1 = unit_proj(1)*scan_x + unit_proj(2)*scan_y + unit_proj(3)*scan_z == cosd(angle);
        eqn2 = unit_norm(1)*scan_x + unit_norm(2)*scan_y + unit_norm(3)*scan_z == 0;
        
        
        [A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [scan_x, scan_y, scan_z]);
        
        solutions = linsolve(A,B);
        
        % pull intensity values
        search_vector = solutions/norm(solutions); % define search direction
        search_indices = round(repmat(spline_location,1,radial_search_dist) + ...
            times(repmat(search_vector,1,radial_search_dist),repmat(1:radial_search_dist,3,1))); % specify search coords
        
        % this line variable depending on the shape of bwcellimage
        search_indices = search_indices([2 1 3],:);
        % 
        
        all_indices(w,angle_iteration+1,1:size(search_indices,1),1:size(search_indices,2)) = search_indices;
        
        search_indices(:,search_indices(3,:)>size(bwcellimage,3))=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(2,:)>size(bwcellimage,2))=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(1,:)>size(bwcellimage,1))=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(3,:)<1)=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(2,:)<1)=[]; % remove out-of-bounds coords
        search_indices(:,search_indices(1,:)<1)=[]; % remove out-of-bounds coords
        % if size(search_indices,2) ~= 0
        try
            ind = sub2ind(size(bwcellimage),search_indices(1,:),search_indices(2,:),search_indices(3,:));
            % int_averages(angle_iteration+1) = mean(bwcellimage(ind(~isnan(ind)))); % pull int values and average them
            ints(angle_iteration+1) = mode(bwcellimage(ind(~isnan(ind)))); % or pull max int values
        % else
        catch
            % need to make an exception for inds that are outside of the
            % image, usually this means that the ROI is hitting the edge of
            % the image (need to take a bigger stack).
            
            % warning('Radial indices are out of bounds at position [%i,%i] of %i. Skipping...',W,angle_iteration*angle_delta,scurve.pieces);
        end
        
    end
    
    % unr_ave_int_values(:,W) = int_averages;
    unr_image(:,W) = ints;
    
    base_ref_vector = proj_ref_vector; % use this to dynamically change the ref vector on each cross section
    
%    time(W) = toc;
    
end

else
%% unroll using pre-calculated indices
% (saves a lot of time from calculating cross-sections again)


unr_image = zeros(size(all_indices,2),size(all_indices,1));


for i = 1:size(all_indices,1)
    for j = 1:size(all_indices,2)
        ints = zeros(1,40);
        for k = 1:40
            % if any(all_indices(i,j,:,k))
            try
                ints(k) = bwcellimage(all_indices(i,j,1,k),all_indices(i,j,2,k),all_indices(i,j,3,k));
            catch
            end
            % end
        end
% debug use only
%         if sum(ints)==0
%             fprintf('  Skipped position [%i %i] due to out of bounds.\n',i,j);
%         end
        unr_image(j,i) = mode(nonzeros(ints));
        
    end
end

end

unr_image(isnan(unr_image))=0;
% plot
% figure();imshow(unr_image,[]);

%% 
% % calculate pixel-indexed orientations
% 
% aa = imdilate(imerode(repmat(unr_image,3,1),strel('disk',1,4)),strel('disk',1,4));
% cc = regionprops(aa,'filledimage');
% bb = regionprops(cc.FilledImage,'centroid','orientation','area','majoraxislength','minoraxislength','equivdiameter','eccentricity');
% unr_regionprops = zeros(size(bb,1),8);
% 
% for i = 1:size(bb,1)
%     unr_regionprops(i,[1 2]) = bb(i).Centroid;
%     unr_regionprops(i,3) = bb(i).Orientation;
%     unr_regionprops(i,4) = bb(i).Area;
%     unr_regionprops(i,5) = bb(i).MajorAxisLength;
%     unr_regionprops(i,6) = bb(i).MinorAxisLength;
%     unr_regionprops(i,7) = bb(i).EquivDiameter;
%     unr_regionprops(i,8) = bb(i).Eccentricity;
% end
% 
% unr_or_values = zeros(size(unr_regionprops,1),4);
% % y vector component corrected for circumferential distortion
% unr_or_values(:,1) = unr_regionprops(:,5).*sind(unr_regionprops(:,3)).*(48.9794+10.7863*unr_regionprops(:,1)/max(unr_regionprops(:,1))*7+3)*.151/90;
% % x vector component assuming length = .151
% unr_or_values(:,2) = unr_regionprops(:,5).*cosd(unr_regionprops(:,3)).*.151;
% % orientation
% unr_or_values(:,3) = atand(unr_or_values(:,1)./unr_or_values(:,2));
% % length?
% unr_or_values(:,4) = sqrt(unr_or_values(:,1).^2+unr_or_values(:,2).^2);
% 
% % figure();imshow(imfill(aa));


% hold on; quiver(unr_max_regionprops(:,1)-unr_or_values(:,2)*2,unr_max_regionprops(:,2)-unr_or_values(:,1)*2,unr_or_values(:,2),unr_or_values(:,1),0.8);
% quiver(unr_max_regionprops(:,1)-unr_or_values(:,2)*2,unr_max_regionprops(:,2)-unr_or_values(:,1)*2,cosd(unr_max_regionprops(:,3)),sind(unr_max_regionprops(:,3)),0.4);

%%
% % calculate normalized x positions (requires branchid and branchposition
% % data from allFiles)
% 
% branchid = allFiles.(userDefinedallFilesName).branchid;
% branchposition = allFiles.(userDefinedallFilesName).branchposition;
% 
% x_pos_adj = nan(1,branchposition(1)-1);
% for i = 1:numel(branchid)-1
%     x_pos_adj = [x_pos_adj, branchid(i):(branchid(i+1)-branchid(i))/(branchposition(i+1)-branchposition(i)):branchid(i+1)];
%     x_pos_adj(end) = [];
% end
% x_pos_adj(end+1) = branchid(end);
% 
% 
% % calculate data
% 
% clear unr_calc_data
% 
% unr_calc_data(1,:) = x_pos_adj(floor(unr_regionprops(unr_regionprops(:,1)<length(x_pos_adj))));
% unr_calc_data(2,:) = unr_or_values(floor(unr_regionprops(:,1))<length(x_pos_adj),3)';
% 
% 
% % plot orientations of unrolled objects
% numofFiles = 1;
% lower_bound = 1;
% step_size = 0.1;
% upper_bound = 10;
% 
%     stats_x_position = lower_bound:step_size:upper_bound;
%     stats_unror_value = zeros(numofFiles,size(stats_x_position,2));
%     stats_unror_sterr = zeros(numofFiles,size(stats_x_position,2));
%     for j = 1:size(stats_x_position,2)
%         %data source
%         aa = unr_calc_data(:, unr_calc_data(1,:)>stats_x_position(j)-0.5);
%         bb = aa(:, aa(1,:)<stats_x_position(j)+0.5);
%         %data location
%         stats_unror_value(1,j) = mean(bb(2,:));
%         stats_unror_sterr(1,j) = std(bb(2,:))/numel(bb(2,:))^(1/2);
%     end


%% save

imwrite(logical(unr_image),[userDefinedallFilesName,'unrimage'],'tiff');
    save([userDefinedallFilesName,'unrdata'],'unr_image','all_indices');


