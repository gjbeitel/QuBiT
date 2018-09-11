function s2_display_tubes_f(allFiles)

%% Display tube structure
if ~isstruct(allFiles)
    error('The input argument must be a struct of the H5 data.')
end

FileList = fieldnames(allFiles);

colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
    [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
    [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1]];



% Show the iso-surface of the vessels
figure('name', FileList{1});

FV = isosurface(allFiles.(FileList{1}).h5bwfilled, 0.5);
patch(FV,'facecolor',[0.4 0.6 0.9],'facealpha',0.3,'edgecolor','none');
view(3)
camlight
axis equal


% Display the skeleton
hold on;
% legendInfo = cell(length(allFiles.(varname{fileNo}).S),1);
for i=1:length(allFiles.(FileList{1}).S)
    L=allFiles.(FileList{1}).S{i};
    % plot3(L(:,1),L(:,2),L(:,3),'-','Color',rand(1,3)); % Plot skeleton
    plot3(L(:,1),L(:,2),L(:,3),'-','color',colorspec(mod(i,size(colorspec,1))+1,:),'linewidth',4);
    plot3(L(1,1),L(1,2),L(1,3), '*', 'MarkerSize', 12); % Plot nodes
    
    % Label the respective nodes with segment numbers.
    % Exception included in case of < 10 points in a segment to avoid
    % "ind exceeds dims" error when displaying segment #'s.
    text(L(min([ceil(size(L,1)/2),10]),1), L(min([ceil(size(L,1)/2),10]),2), L(min([ceil(size(L,1)/2),10]),3),num2str(i),'FontSize',20);
    % legendInfo{i} = ['X = ' num2str(i)];
end

% legend(legendInfo)


return
%% Debug use only

% Show the iso-surface of the vessels

FV = isosurface(allFiles.(varname{fileNo}).h5bwfilled, 0.5);
figure('color',[0 0 0]);set(gca,'color',[0 0 0]);axis equal;axis([-inf inf -inf inf -inf inf]);
patch(FV,'facecolor',[1 1 1],'facealpha',0.3,'edgecolor','none');

% Display the skeleton
hold on;
for i=1:length(allFiles.(varname{fileNo}).S)
    L=allFiles.(varname{fileNo}).S{i};
    plot3(L(:,1),L(:,2),L(:,3),'-','color',colorspec(mod(i,size(colorspec,1))+1,:));
    %plot3(L(1,1),L(1,2),L(1,3), '*', 'MarkerSize', 12); % Plot nodes
    %text(L(min([ceil(size(L,1)/2),10]),1), L(min([ceil(size(L,1)/2),10]),2), L(min([ceil(size(L,1)/2),10]),3),num2str(i));
end

end

