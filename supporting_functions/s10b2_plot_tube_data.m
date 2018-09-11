%% plot cross section and segment length and tube perim
% this should come from the aggregated tubedata, e.g.
% <genotype>tubedata.mat

% seg_start = 1;
% seg_end = 10;

seg_start = input('Segment start: ');
seg_end = input('Segment end: ');
step_size = input('Interval length: ');

[file_array, ~, ~] = uigetfile('*.*', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
    
end

f1 = figure('name','Length data');
f2 = figure('name','Cross-section data');
f6 = figure('name','Surface Area data');
f22 = figure('name','Cross-section v2');
f26 = figure('name','Perim v2');
f27 = figure('name','Radius');
f28 = figure('name','Laplacian');

obj_array = gobjects(numofFiles,7,4);

for i = 1:numofFiles
    
    clear length_data cross_data perim_data
    
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)))
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    
    if exist('length_data','var') == 0 || exist('perim_data','var') == 0 || exist('cross_data','var') == 0
        warning('Data missing. Skipping...')
        continue
    end

    figure(f1);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:seg_end-0.5;
    y_data = mean(length_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan');
    y_sterr = std(length_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan')./sqrt(sum(isfinite(length_data(2:end,seg_start:seg_end-1)),1));
    obj_array(i,1,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,1,2) = plot(x_data,y_data,'LineStyle','none','Marker','o','MarkerSize',9,'MarkerFaceColor',[r(i) g(i) b(i)],'Color',[r(i) g(i) b(i)]);
    obj_array(i,1,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,1,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f2);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:seg_end-0.5;
    y_data = mean(cross_data(2:end,seg_start:seg_end-1)*.151^2,1,'omitnan');
    y_sterr = std(cross_data(2:end,seg_start:seg_end-1)*.151^2,1,'omitnan')./sqrt(sum(isfinite(length_data(2:end,seg_start:seg_end-1)),1));
    obj_array(i,2,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,2,2) = plot(x_data,y_data,'LineStyle','none','Marker','o','MarkerSize',9,'MarkerFaceColor',[r(i) g(i) b(i)],'Color',[r(i) g(i) b(i)]);
    obj_array(i,2,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,2,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f6);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:seg_end-0.5;
    y_data = mean(perim_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan').*mean(length_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan');
    y_sterr = (mean(perim_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan').*mean(length_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan')).* sqrt(...
         (std(perim_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan')/mean(perim_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan')).^2 ...
        +(std(length_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan')/mean(length_data(2:end,seg_start:seg_end-1)*.151,1,'omitnan')).^2);
    % propagation of error (error = v1*v2 * sqrt((e1/v1)^2 + (e2/v2)^2)
    obj_array(i,3,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,3,2) = plot(x_data,y_data,'LineStyle','none','Marker','o','MarkerSize',9,'MarkerFaceColor',[r(i) g(i) b(i)],'Color',[r(i) g(i) b(i)]);
    obj_array(i,3,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,3,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
        figure(f22);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:step_size:seg_end-0.5;
    y_data = nan(1,length(x_data)-1);
    y_std = nan(1,length(x_data)-1);
    for j = 1:length(x_data)-1
    y_data(j) = .151^2*mean(cs_aggregate(cs_aggregate(:,1)>x_data(j) & cs_aggregate(:,1)<x_data(j+1),2));
    y_std(j) = std(y_data,'omitnan')/length(cs_aggregate(cs_aggregate(:,1)>x_data(j) & cs_aggregate(:,1)<x_data(j+1),2))^(1/2);
    end
    y_data = smooth(y_data)';
    y_std = smooth(y_std)';
    x_data = x_data(1:end-1)+step_size/2;
    % propagation of error (error = v1*v2 * sqrt((e1/v1)^2 + (e2/v2)^2)
    obj_array(i,4,1) = fill([x_data';flipud(x_data')],[(y_data+y_std)';flipud((y_data-y_std)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,4,2) = plot(x_data,y_data,'LineWidth',2,'Marker','none','Color',[r(i) g(i) b(i)]);
    obj_array(i,4,3) = plot(x_data,y_data+y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,4,4) = plot(x_data,y_data-y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    
        figure(f26);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:step_size:seg_end-0.5;
    y_data = nan(1,length(x_data)-1);
    y_std = nan(1,length(x_data)-1);
    for j = 1:length(x_data)-1
    y_data(j) = mean(cs_aggregate(cs_aggregate(:,1)>x_data(j) & cs_aggregate(:,1)<x_data(j+1),3));
    y_std(j) = std(y_data,'omitnan');
    end
    x_data = x_data(1:end-1)+step_size/2;
    % propagation of error (error = v1*v2 * sqrt((e1/v1)^2 + (e2/v2)^2)
    obj_array(i,5,1) = fill([x_data';flipud(x_data')],[(y_data+y_std)';flipud((y_data-y_std)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,5,2) = plot(x_data,y_data,'LineStyle','none','Marker','o','MarkerSize',9,'MarkerFaceColor',[r(i) g(i) b(i)],'Color',[r(i) g(i) b(i)]);
    obj_array(i,5,3) = plot(x_data,y_data+y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,5,4) = plot(x_data,y_data-y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    figure(f27);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:step_size:seg_end-0.5;
    y_data = nan(1,length(x_data)-1);
    y_std = nan(1,length(x_data)-1);
    for j = 1:length(x_data)-1
    y_data(j) = (mean(cs_aggregate(cs_aggregate(:,1)>x_data(j) & cs_aggregate(:,1)<x_data(j+1),3)))/(2*pi)*.151;
    y_std(j) = std(y_data,'omitnan');
    end
    y_data = smooth(y_data)';
    y_std = smooth(y_std)';
    x_data = x_data(1:end-1)+step_size/2;
    % propagation of error (error = v1*v2 * sqrt((e1/v1)^2 + (e2/v2)^2)
    obj_array(i,6,1) = fill([x_data';flipud(x_data')],[(y_data+y_std)';flipud((y_data-y_std)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,6,2) = plot(x_data,y_data,'LineWidth',2,'Marker','none','Color',[r(i) g(i) b(i)]);
    obj_array(i,6,3) = plot(x_data,y_data+y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,6,4) = plot(x_data,y_data-y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
    
        figure(f28);hold on;
    % x_data = 1.5:9.5;
    x_data = seg_start+0.5:step_size:seg_end-0.5;
    y_data = nan(1,length(x_data)-1);
    % y_std = nan(1,length(x_data)-1);
    for j = 1:length(x_data)-1
    y_data(j) = (mean(cs_aggregate(cs_aggregate(:,1)>x_data(j) & cs_aggregate(:,1)<x_data(j+1),3)))/(2*pi)*.151;
    % y_std(j) = std(y_data,'omitnan');
    end
    y_data = smooth(del2(y_data))';
    % y_std = smooth(y_std)';
    x_data = x_data(1:end-1)+step_size/2;
    % propagation of error (error = v1*v2 * sqrt((e1/v1)^2 + (e2/v2)^2)
    % obj_array(i,7,1) = fill([x_data';flipud(x_data')],[(y_data+y_std)';flipud((y_data-y_std)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,7,2) = plot(x_data,y_data,'LineWidth',2,'Marker','none','Color',[r(i) g(i) b(i)]);
    % obj_array(i,7,3) = plot(x_data,y_data+y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % obj_array(i,7,4) = plot(x_data,y_data-y_std,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    % errorbar(x_data,y_data,y_sterr,'LineStyle','none','Color',[r(i) g(i) b(i)]);
    
end


% legend_array = cell(1,size(file_array,2)*2);
% for i = 1:size(file_array,2)
%     try legend_array{2*i-1} = file_array{i}(1:end-12);
%     catch
%         legend_array{2*i-1} = file_array(1:end-12);
%     end
%     legend_array{2*i} = '';
% end

legend_array = cell(1,numofFiles);
for i = 1:numofFiles
    try legend_array{i} = file_array{i}(1:end-12);
    catch
        legend_array{i} = file_array(1:end-12);
    end
end

figure(f1);legend(obj_array(:,1,2),legend_array,'location','best');
% figure(f1);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('Dorsal trunk segment length','FontSize',32);
axis([seg_start seg_end 0 inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Segment length (\mum)','FontSize',24);

figure(f2);legend(obj_array(:,2,2),legend_array,'location','best');
%figure(f2);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('Dorsal trunk cross-sectional area','FontSize',32);
axis([seg_start seg_end 0 inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Cross-sectional area (\mum^2)','FontSize',24);

figure(f6);legend(obj_array(:,3,2),legend_array,'location','best');
% figure(f6);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('Dorsal trunk surface area','FontSize',32);
axis([seg_start seg_end 0 inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Surface area (\mum^2)','FontSize',24);

figure(f22);legend(obj_array(:,4,2),legend_array,'location','best');
% figure(f6);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('DT cross-sectional area','FontSize',32);
axis([seg_start seg_end 0 inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Cross-sectional area (\mum^2)','FontSize',24);

figure(f26);legend(obj_array(:,5,2),legend_array,'location','best');
% figure(f6);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('CS Perimeter v2','FontSize',32);
axis([seg_start seg_end 0 inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Perimeter (pix)','FontSize',24);

figure(f27);legend(obj_array(:,6,2),legend_array,'location','best');
% figure(f6);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('Tube Radius','FontSize',32);
axis([seg_start seg_end 0 inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('Radius (\mum)','FontSize',24);

figure(f28);legend(obj_array(:,7,2),legend_array,'location','best');
% figure(f6);legend(legend_array,'location','best');
set(gca,'FontSize',24);
title('Laplacian','FontSize',32);
axis([seg_start seg_end -inf inf]);
xlabel('Normalized position (TC)','FontSize',24);
ylabel('(\mum/S^2)','FontSize',24);