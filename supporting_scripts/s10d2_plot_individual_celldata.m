%% plot moving window average of (absolute) orientation, cell surface area, and aspect ratio vs. position

% this should come from the aggregated celldata, e.g.
% <genotype>celldata.mat


% seg_start = 3;
% seg_end = 10;
step_size = 10;

seg_start = input('Segment start: ');
seg_end = input('Segment end: ');

lower_bound = seg_start*100+25;
upper_bound = seg_end*100-25;


[file_array, ~, ~] = uigetfile('*.mat', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
end

% figure('color',[0 0 0]);set(gca,'color',[0 0 0]);hold on
% dark bg figure base code

f4 = figure('name','Cell surface area vs position');
f5 = figure('name','Cell aspect ratio vs position');
f3 = figure('name','Cell orientation vs position');
% f7 = figure('name','Cell orientation vs position');

obj_array = gobjects(numofFiles,4,4);

    stats_x_position = lower_bound:step_size:upper_bound;
    stats_or_value = zeros(numofFiles,size(stats_x_position,2));
    stats_or_sterr = zeros(numofFiles,size(stats_x_position,2));
    stats_sa_value = zeros(numofFiles,size(stats_x_position,2));
    stats_sa_sterr = zeros(numofFiles,size(stats_x_position,2));
    stats_ar_value = zeros(numofFiles,size(stats_x_position,2));
    stats_ar_sterr = zeros(numofFiles,size(stats_x_position,2));

for i = 1:numofFiles
    
    clear lm98;
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    
    if exist('lm98','var') == 0
        warning('Data missing. Skipping...');
        continue;
    end
    if size(lm98,2) < 11
        warning('Data missing. Skipping...');
        continue;
    end

    for j = 1:size(stats_x_position,2)
        aa = lm98(lm98(:,11)>(stats_x_position(j)-20)/100,:);
        bb = aa(aa(:,11)<(stats_x_position(j)+20)/100,:);
        % orientation
        stats_or_value(i,j) = mean(abs(bb(:,7)));
        stats_or_sterr(i,j) = std(abs(bb(:,7)))/numel(bb(:,7))^(1/2);
        % surface area
        stats_sa_value(i,j) = mean(abs(bb(:,6)))*.151*.151;
        stats_sa_sterr(i,j) = std(abs(bb(:,6))*.151*.151)/numel(bb(:,6))^(1/2);
        % aspect ratio
        cc = bb(:,5);
        cc(cc<10^-5) = [];
        cc = 1./cc;
        stats_ar_value(i,j) = mean(cc);
        stats_ar_sterr(i,j) = std(cc)/numel(cc)^(1/2);
    end
    
    % plot all
    figure(f3);hold on;
    x_data = stats_x_position/100;
    y_data = smooth(stats_or_value(i,:))';
    y_sterr = smooth(stats_or_sterr(i,:))';
    obj_array(i,1,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,1,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,1,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,1,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f4);hold on;
    x_data = stats_x_position/100;
    y_data = smooth(stats_sa_value(i,:))';
    y_sterr = smooth(stats_sa_sterr(i,:))';
    obj_array(i,2,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,2,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,2,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,2,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f5);hold on;
    x_data = stats_x_position/100;
    y_data = smooth(stats_ar_value(i,:))';
    y_sterr = smooth(stats_ar_sterr(i,:))';
    obj_array(i,3,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,3,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,3,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,3,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    % figure(f7);hold on;
    % obj_array(i,4,1) = scatter(lm98(:,18),abs(lm98(:,13)),'markeredgecolor',[r(i)/2 g(i)/2 b(i)/2],'markerfacecolor',[r(i) g(i) b(i)]);
    
end

if numofFiles > 1
[~,p_or] = ttest(stats_or_value(1,:),stats_or_value(2,:),'alpha',.01);
[~,p_ar] = ttest(stats_ar_value(1,:),stats_ar_value(2,:),'alpha',.01);
[~,p_sa] = ttest(stats_sa_value(1,:),stats_sa_value(2,:),'alpha',.01);
end

legend_array = cell(1,numofFiles);
for i = 1:numofFiles
    try legend_array{i} = file_array{i}(1:end-12);
    catch
        legend_array{i} = file_array(1:end-12);
    end
end

figure(f4)
set(gca,'FontSize',24);
title('Cell apical surface area','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Cell apical surface area (\mum^2)','FontSize',24);
axis([seg_start seg_end 0 inf]);
legend(obj_array(:,1,2),legend_array,'location','best');

figure(f5)
set(gca,'FontSize',24);
title('Cell apical aspect ratio','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Cell apical aspect ratio (major:minor axis)','FontSize',24);
axis([seg_start seg_end 1 3]);
legend(obj_array(:,2,2),legend_array,'location','best');

figure(f3)
set(gca,'FontSize',24);
title('Cell orientation','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Cell orientation relative to DT axis (^{\circ})','FontSize',24);
axis([seg_start seg_end 0 90]);
legend(obj_array(:,3,2),legend_array,'location','best');

% figure(f7)
% set(gca,'FontSize',12);
% title('Cell orientation','FontSize',20);
% xlabel('Normalized position (TC)','FontSize',16);
% ylabel('Cell orientation relative to DT axis (^{\circ})','FontSize',16);
% axis([seg_start seg_end 0 90]); 
% legend(obj_array(:,4,1),legend_array,'location','best');