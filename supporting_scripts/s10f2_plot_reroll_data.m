%% import/aggregate 2d cell data for a particular genotype
% this should come from the unrdata, e.g.
% <genotype><replicate_number>unrdata.mat

% aggregate 2d cell data

seg_start = input('Segment start: ');
seg_end = input('Segment end: ');

lower_bound = seg_start+.25;
upper_bound = seg_end-.25;
step_size = .1;

[file_array, ~, ~] = uigetfile('*.mat', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
end


f30 = figure('name','Cell Size');
f31 = figure('name','Aspect Ratio');
f32 = figure('name','Orientation');

obj_array = gobjects(numofFiles,3,4);

    stats_x_position = lower_bound:step_size:upper_bound;
    stats_sa_value = zeros(numofFiles,size(stats_x_position,2));
    stats_sa_sterr = zeros(numofFiles,size(stats_x_position,2));
    stats_ar_value = zeros(numofFiles,size(stats_x_position,2));
    stats_ar_sterr = zeros(numofFiles,size(stats_x_position,2));
    stats_or_value = zeros(numofFiles,size(stats_x_position,2));
    stats_or_sterr = zeros(numofFiles,size(stats_x_position,2));

% file_array = input('Files array {string1 ,string2...}: ');


for i = 1:numofFiles
    
clear reroll_agg
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    
    if ~exist('reroll_agg','var')
        warning('Data not found in %s. Skipping...',varname);
        continue
    end
    
    for j = 1:length(stats_x_position)
        aa = reroll_agg(reroll_agg(:,1)<stats_x_position(j)+.2 & reroll_agg(:,1)>stats_x_position(j)-.2,:);
        
        % surface area
        stats_sa_value(i,j) = mean(abs(aa(:,2)))*.151*.151;
        stats_sa_sterr(i,j) = std(abs(aa(:,2))*.151*.151)/numel(aa(:,2))^(1/2);
        % aspect ratio
        stats_ar_value(i,j) = mean(abs(aa(:,3)));
        stats_ar_sterr(i,j) = std(abs(aa(:,3)))/numel(aa(:,3))^(1/2);
        % orientation
        stats_or_value(i,j) = mean(abs(aa(:,4)));
        stats_or_sterr(i,j) = std(abs(aa(:,4)))/numel(aa(:,4))^(1/2);
        
    end
    % plot all
    
    figure(f30);hold on;
    x_data = stats_x_position;
    y_data = smooth(stats_sa_value(i,:))';
    y_sterr = smooth(stats_sa_sterr(i,:))';
    obj_array(i,2,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,2,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,2,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,2,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f31);hold on;
    x_data = stats_x_position;
    y_data = smooth(stats_ar_value(i,:))';
    y_sterr = smooth(stats_ar_sterr(i,:))';
    obj_array(i,3,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,3,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,3,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,3,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f32);hold on;
    x_data = stats_x_position;
    y_data = smooth(stats_or_value(i,:))';
    y_sterr = smooth(stats_or_sterr(i,:))';
    obj_array(i,1,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,1,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,1,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,1,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
end

if numofFiles == 2
[~,p_or] = ttest(stats_or_value(1,:),stats_or_value(2,:),'alpha',.01);
[~,p_ar] = ttest(stats_ar_value(1,:),stats_ar_value(2,:),'alpha',.01);
[~,p_sa] = ttest(stats_sa_value(1,:),stats_sa_value(2,:),'alpha',.01);
end

legend_array = cell(1,numofFiles);
for i = 1:numofFiles
    try legend_array{i} = file_array{i}(1:end-10);
    catch
        legend_array{i} = file_array(1:end-10);
    end
end

figure(f30)
set(gca,'FontSize',24);
title('Cell apical surface area','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Area (\mum^2)','FontSize',24);
axis([seg_start seg_end 0 inf]);
legend(obj_array(:,1,2),legend_array,'location','best');

figure(f31)
set(gca,'FontSize',24);
title('Cell apical aspect ratio','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Aspect Ratio','FontSize',24);
axis([seg_start seg_end 1 3]);
legend(obj_array(:,2,2),legend_array,'location','best');

figure(f32)
set(gca,'FontSize',24);
title('Cell orientation','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Angle w.r.t. local DT','FontSize',24);
axis([seg_start seg_end 0 90]);
legend(obj_array(:,3,2),legend_array,'location','best');

% figure(f7)
% set(gca,'FontSize',12);
% title('Cell orientation','FontSize',20);
% xlabel('Normalized position (TC)','FontSize',16);
% ylabel('Cell orientation relative to DT axis (^{\circ})','FontSize',16);
% axis([seg_start seg_end 0 90]); 
% legend(obj_array(:,4,1),legend_array,'location','best');