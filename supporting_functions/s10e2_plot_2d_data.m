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


f11 = figure('name','Cell connectivity vs position');
f12 = figure('name','Cells per cross-section vs position');

obj_array = gobjects(numofFiles,2,6);

stats_x_position = lower_bound:step_size:upper_bound;
stats_cc_value = zeros(numofFiles,size(stats_x_position,2),5);
stats_cpc_value = zeros(numofFiles,size(stats_x_position,2));
stats_cpc_sterr = zeros(numofFiles,size(stats_x_position,2));

% file_array = input('Files array {string1 ,string2...}: ');


for i = 1:numofFiles
    
clear cc cpc;
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    
    if ~exist('cc','var') || ~exist('cpc','var')
        warning('Data not found in %s. Skipping...',varname);
        continue
    end
    
    for j = 1:length(stats_x_position)
        % cc
        aa = cc(cc(:,1)<stats_x_position(j)+.5 & cc(:,1)>stats_x_position(j)-.5,:);
        % stats_cc_value(i,j) = mean(aa(:,2));
        % stats_cc_sterr(i,j) = std(aa(:,2))/length(aa)^(1/2);
          
         stats_cc_value(i,j,1) = sum(aa(:,2)<=4);
        stats_cc_value(i,j,2) = sum(aa(:,2)==5);
        stats_cc_value(i,j,3) = sum(aa(:,2)==6);
        stats_cc_value(i,j,4) = sum(aa(:,2)==7);
         stats_cc_value(i,j,5) = sum(aa(:,2)>=8);
          
        stats_cc_value(i,j,:) = stats_cc_value(i,j,:)/sum(stats_cc_value(i,j,:));
        % cpc
        aa = cpc(cpc(:,1)<stats_x_position(j)+.25 & cpc(:,1)>stats_x_position(j)-.25,:);
        stats_cpc_value(i,j) = mean(aa(:,2));
        stats_cpc_sterr(i,j) = std(aa(:,2))/(length(aa)/3)^(1/2);
        % stats_cpc_sterr(i,j) = std(aa(:,2));
    end
    
    figure(f11);hold on;
%     x_data = stats_x_position;
%     y_data = stats_cc_value(i,:);
%     y_sterr = stats_cc_sterr(i,:);
%     obj_array(i,1,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
%     obj_array(i,1,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
%     obj_array(i,1,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
%     obj_array(i,1,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,1,2:6) = bar(stats_x_position,squeeze(stats_cc_value(i,:,:)),'stacked','FaceAlpha',.5);

    
    
    figure(f12);hold on;
    x_data = stats_x_position;
    y_data = smooth(stats_cpc_value(i,:))';
    y_sterr = smooth(stats_cpc_sterr(i,:))';
    obj_array(i,2,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,2,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,2,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,2,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
end

if i > 1
    % [~,p_cc] = ttest2(stats_cc_value(1,:),stats_cc_value(2,:),'alpha',.01);
    [~,p_cpc] = ttest2(stats_cpc_value(1,:),stats_cpc_value(2,:),'alpha',.01);
end

legend_array = cell(1,numofFiles);
for i = 1:numofFiles
    try legend_array{i} = file_array{i}(1:end-11);
    catch
        legend_array{i} = file_array(1:end-11);
    end
end

figure(f11)
set(gca,'FontSize',24);
title('Cell connectivity vs position','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Proportion of cells','FontSize',24);
axis([seg_start seg_end -inf inf]);
% legend(obj_array(:,1,2),legend_array,'location','best');
legend(obj_array(1,1,2:6),{'4-gon','5-gon','6-gon','7-gon','8-gon'},'location','best');

figure(f12)
set(gca,'FontSize',24);
title('Cells per cross-section','FontSize',32);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Cell number','FontSize',24);
axis([seg_start seg_end -inf inf]);
legend(obj_array(:,2,2),legend_array,'location','best');
