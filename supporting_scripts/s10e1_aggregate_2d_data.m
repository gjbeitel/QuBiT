%% import/aggregate 2d cell data for a particular genotype
% this should come from the celldata, e.g.
% <genotype><replicate_number>unrdata.mat

% aggregate 2d cell data



[file_array, ~, ~] = uigetfile('*.mat', 'Select .mat file(s).', 'MultiSelect', 'on');

if iscell(file_array)
    numofFiles = length(file_array);
elseif ischar(file_array)
    numofFiles = 1;
elseif (file_array == 0)
    error('No data selected.');
end

fprintf('\nSave aggregate data as? (Leave blank to skip saving)\n');
save_name = input('','s');

% file_array = input('Files array {string1 ,string2...}: ');
ilastik_rp_2d_aggregate = [];
cpc = [];
cc = [];


f13 = figure('name','Cells per cross-section vs position (individual tube data)');
obj_array = gobjects(numofFiles,4);

for i = 1:numofFiles
    clear ilastlik_rp_2d
    clear ilastik_cc ilastik_cpc ilastik_centroids ilastik_n_position
    if numofFiles > 1
        fprintf('Reading %s... (file %i of %i)\n',char(file_array(i)),i,numofFiles);
        load(char(file_array(i)));
        varname = char(file_array(i));
    else
        fprintf('Reading %s... (file %i of %i)\n',file_array,i,numofFiles);
        load(file_array);
        varname = file_array;
    end
    if ~exist('ilastik_cc','var') || ~exist('ilastik_cpc','var') || ~exist('ilastik_centroids','var') || ~exist('ilastik_rp_2d','var') || size(ilastik_rp_2d,2) ~= 7
        warning('Data not found in %s. Skipping...',varname);
        continue
    end

    ilastik_rp_2d_aggregate = [ilastik_rp_2d_aggregate;cell2mat(ilastik_rp_2d(2:end,:))];
    cpc = [cpc;[ilastik_n_position,ilastik_cpc(end-length(ilastik_n_position)+1:end)]]; % this is fine
    cc = [cc;[ilastik_centroids(ilastik_centroids(:,2)>90 & ilastik_centroids(:,2)<180,3),ilastik_cc(ilastik_centroids(:,2)>90 & ilastik_centroids(:,2)<180)]];
    
    
    
    
    % 180322

lower_bound = ilastik_n_position(1);step_size = 0.1;upper_bound = ilastik_n_position(end);
stats_x_position = lower_bound:step_size:upper_bound;
stats_cpc_value = zeros(1,size(stats_x_position,2));
stats_cpc_sterr = zeros(1,size(stats_x_position,2));

    for j = 1:length(stats_x_position)
        % cpc
        aa = cpc(cpc(:,1)<stats_x_position(j)+.25 & cpc(:,1)>stats_x_position(j)-.25,:);
        stats_cpc_value(1,j) = mean(aa(:,2));
        stats_cpc_sterr(1,j) = std(aa(:,2))/length(aa)^(1/2);
    end
    
    figure(f13);hold on;
    x_data = stats_x_position;
    y_data = stats_cpc_value;
    y_sterr = stats_cpc_sterr;
    obj_array(i,1) = fill([x_data';flipud(x_data')],[(y_data+y_sterr)';flipud((y_data-y_sterr)')],[(r(i)+2)/3 (g(i)+2)/3 (b(i)+2)/3],'linestyle','none','facealpha',0.3);
    obj_array(i,2) = plot(x_data,y_data,'LineWidth',2,'Color',[r(i) g(i) b(i)]);
    obj_array(i,3) = plot(x_data,y_data+y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    obj_array(i,4) = plot(x_data,y_data-y_sterr,'Color',[(r(i)+1)/2 (g(i)+1)/2 (b(i)+1)/2]);
    
end

legend_array = cell(1,numofFiles);
for i = 1:numofFiles
    try legend_array{i} = file_array{i}(1:end-11);
    catch
        legend_array{i} = file_array(1:end-11);
    end
end


figure(f13)
set(gca,'FontSize',12);
title('Cells per cross-section','FontSize',16);
xlabel('Normalized position (TC)','FontSize',12);
ylabel('Cell number','FontSize',12);
% axis([1 10 -inf inf]);
legend(obj_array(isgraphics(obj_array(:,2)),2),legend_array(isgraphics(obj_array(:,2))),'location','best');


if ~isempty(save_name)
save(save_name,'cc','cpc','ilastik_rp_2d_aggregate');    
end

fprintf('Done.\n');
