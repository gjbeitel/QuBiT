

% 1: Cell ID
% 2-4: Centroid co-ords
% 5-6: Eigenvalues
% 7-8: Eccentricity
% 9: Aspect ratio
% 10: Rel. position
% 11: Rel. size
% 12: Cell dist. from centerline
% 13: Cell orientation (angle)
% 14-16: Cell orientation (vector)
% 17: Circumferential position


%% Plot filtered vectors as 2d map
%  180714: Probably better just to do unrolling first instead of this bs
% 
% fprintf('Plot vectors...\n');
% lm97 = lm98;
% 
% for i = 1:size(lm97,1)
%     if lm97(i,7)>90
%         lm97(i,7) = 180-lm97(i,7);
%     end
% end
% 
% figure();axis equal;axis([-inf inf -180 180]);hold on
% colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
%     [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
%     [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1];[0.2 0.5 0.8];[0.2 0.8 0.5];...
%     [0.8 0.2 0.5];[0.8 0.5 0.2];[0.5 0.8 0.2];[0.5 0.2 0.8]];
% for i = 1:size(lm97,1)
%     quiver(lm97(i,10),... % X position (rel. position centerline)
%         ...% max(lm97(:,6))*rand(),... % Y position (rel. position cir.)
%         lm97(i,17)*180/pi,... % Y position (rel. position cir.)
%         cosd(lm97(i,7)),... % unit vector direction U (don't change)
%         sind(lm97(i,7)),... % unit vector direction V (don't change)
%         40*lm97(i,9),... % Vector size
%         'color',colorspec(mod(i,size(colorspec,1))+1,:)) % Other properties
% end

%% Plot filtered cells
% This might take a while...

fprintf('Plot filtered cells 3D...\n');
% lm14 = bwconncomp(allFiles.(userDefinedallFilesName).bwcellimage);
 lm14 = bwlabeln(allFiles.(userDefinedallFilesName).bwcellimage);
% bwcellimagef = allFiles.(userDefinedallFilesName).bwcellimagef;

figure();axis equal;axis([-inf inf -inf inf -inf inf]);hold on;
fprintf('Plotting cells...\n');
colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
    [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
    [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1];[0.2 0.5 0.8];[0.2 0.8 0.5];...
    [0.8 0.2 0.5];[0.8 0.5 0.2];[0.5 0.8 0.2];[0.5 0.2 0.8]];

for kk = 1:length(lm98)
% for kk = 1:5
    fprintf('Plotting object with index %3.0f (%3.0f of %3.0f)...\n',lm98(kk,1),kk,length(lm98));
    % lm15 = lm14.PixelIdxList{lm98(kk,1)};
    % lm16 = false(size(allFiles.(userDefinedallFilesName).bwcellimage));
    % lm16(lm15) = true;
%     [I,J,K] = ind2sub(size(lm14),find(lm14==lm98(kk,1)));
%       lm17 = regionprops3([J,I,K],'IsPixList');
%       scatter3(lm17.Centroid(1),lm17.Centroid(2),lm17.Centroid(3),50,'markerfacecolor',colorspec(mod(kk,size(colorspec,1))+1,:),'markeredgecolor','none')
%       quiver3(lm17.Centroid(1)-lm17.FirstAxis(1)*8,lm17.Centroid(2)-lm17.FirstAxis(2)*8,lm17.Centroid(3)-lm17.FirstAxis(3)*8,lm17.FirstAxis(1)*16,lm17.FirstAxis(2)*16,lm17.FirstAxis(3)*16,0,...
%           'color',colorspec(mod(kk,size(colorspec,1))+1,:),'linewidth',3,'showarrowhead','off');
    blah = isosurface(lm14==lm98(kk,1),0.5);
    patch(blah, 'facecolor', colorspec(mod(kk,size(colorspec,1))+1,:), 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
end


return

%% plot cell orientation fields for a single genotype
%  180714: Scrapped. Use unroll tube, but this is here as a reference

% 
% 
% figure('name','Cell orientation vector field');axis equal;axis([-inf inf -180 180]);hold on
% colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
%     [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
%     [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1];[0.2 0.5 0.8];[0.2 0.8 0.5];...
%     [0.8 0.2 0.5];[0.8 0.5 0.2];[0.5 0.8 0.2];[0.5 0.2 0.8]];
% for i = 1:size(lm97,1)
%     quiver(lm97(i,18)*100,... % X position (rel. position/branch)
%         ...%max(lm97(:,17))*rand(),... % Y position (rel. position cir.)
%         lm97(i,17)*180/pi,... % Y position (rel. position cir.)
%         cosd(lm97(i,13)),... % unit vector direction U (don't change)
%         sind(lm97(i,13)),... % unit vector direction V (don't change)
%         20*lm97(i,9),... % Vector size
%         'color',colorspec(mod(i,size(colorspec,1))+1,:)) % Other properties
% end

%% plot normalized cell vectors (rainbow colors based on angle)
%  180714: Scrapped. Use unroll tube, but this is here as a reference

% 
% figure('name','Cell orientation vector field','color',[0 0 0]);set(gca,'color',[0 0 0]);axis equal;axis([-inf inf -180 180]);hold on
% for i = 1:size(lm97,1)
%     if abs(lm97(i,13)) < 10
%         color = [1 0 0];
%     elseif abs(lm97(i,13)) < 20
%         color = [226/255, 87/255, 30/255];
%     elseif abs(lm97(i,13)) < 30
%         color = [1 .5 0];
%     elseif abs(lm97(i,13)) < 40
%         color = [1 1 0];
%     elseif abs(lm97(i,13)) < 50
%         color = [0 1 0];
%     elseif abs(lm97(i,13)) < 60
%         color = [150/255 191/255 51/255];
%     elseif abs(lm97(i,13)) < 70
%         color = [0 0 1];
%     elseif abs(lm97(i,13)) < 80
%         color = [75/255 0 130/255];
%     else
%         color = [139/255 0 1];
%     end
%     quiver(lm97(i,18)*100,... % X position (rel. position/branch)
%         ...% y_values(i),...
%         lm97(i,17)*180/pi,... % Y position (rel. position cir.)
%         cosd(lm97(i,13)),... % unit vector direction U (don't change)
%         sind(lm97(i,13)),... % unit vector direction V (don't change)
%         20*lm97(i,9),... % Vector size
%         'color',color,'linewidth',2) % Other properties
% end
% 
% return

%% Plot filtered cells
% lm14 = bwlabeln(lm13);
% figure();axis equal;axis([-inf inf -inf inf -inf inf]);hold on;
% colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
%     [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
%     [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1];[0.2 0.5 0.8];[0.2 0.8 0.5];...
%     [0.8 0.2 0.5];[0.8 0.5 0.2];[0.5 0.8 0.2];[0.5 0.2 0.8]];
% fprintf('Plotting filtered cells...\n');
% for i = 1:size(lm98,1)
%     kk = lm98(i,1);
%     fprintf('Plotting object %4.0f of %4.0f...\n',kk,max(lm14(:)));
%     lm14(lm14~=kk) = 0;
%     blah = isosurface(lm14,0.3);
%     patch(blah, 'facecolor', colorspec(mod(kk,size(colorspec,1))+1,:), 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
%
%     text(lm98(i,2),lm98(i,3),lm98(i,4),num2str(kk),...
%         'color',colorspec(mod(i,size(colorspec,1))+1,:)); % This labels objects with their IDs.
%     quiver3(lm98(i,2)-5*lm98(i,14)/2,lm98(i,3)-5*lm98(i,15)/2,lm98(i,4)-5*lm98(i,16)/2,...
%         10*lm98(i,14),10*lm98(i,15),10*lm98(i,16)); % Plot primary vector of cell orientation
%
%     lm14 = bwlabeln(lm13);
% end



%% 180713: Various ts / debugging tools below


%% Plot subset of cells that meet conditions

bwcellimage = allFiles.(userDefinedallFilesName).bwcellimage;
lm14 = bwlabeln(bwcellimage);
lm15 = lm14;


figure();hold on;
fprintf('Plotting cells...\n');

for kk = lm98(:,1).'
    
    
    lm15(lm15~=kk) = 0;
    lm15 = logical(lm15);
    
    if sum(lm15(:)) < (1/.151)^2
        continue
    end
    fprintf('Plotting object %i of %i...\n',kk,max(lm98(:,1)));
    if sum(lm15(:)) < (2/.151)^2
        blah = isosurface(lm15,0.3);
        patch(blah, 'facecolor', [1 0 0], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (3/.151)^2
        blah = isosurface(lm15,0.3);
        patch(blah, 'facecolor', [226/255 87/255 30/255], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (4/.151)^2
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [1 .5 0], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (5/.151)^2
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [1 1 0], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (6/.151)^2
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [0 1 0], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (7/.151)^2
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [0 0 1], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (8/.151)^2
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [150/255 191/255 51/255], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    elseif sum(lm15(:)) < (9/.151)^2
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [75/255 0 130/255], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    else
        blah = isosurface(lm14,0.3);
        patch(blah, 'facecolor', [139/255 0 1], 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
    end
    
    
    lm15 = lm14;
end

%% check filtered cells

figure();hold on;
bwcellimage = allFiles.(userDefinedallFilesName).bwcellimage;
lm14 = bwconncomp(bwcellimage);
lm15 = false(lm14.ImageSize);
fprintf('Plotting objects...\n');
for i = 1:lm14.NumObjects
    if size(lm14.PixelIdxList{i},1) > 30
        lm15(lm14.PixelIdxList{i}) = true;
        object = isosurface(lm15,0.3);
        if ismember(i,lm98(:,1))
            fprintf('  True Object  %4.0f of %4.0f...\n',i,lm14.NumObjects);
            patch(object, 'facecolor', [226/255 87/255 30/255], 'facealpha', 0.3, 'edgecolor', 'none');
        else
            fprintf('  False Object %4.0f of %4.0f...\n',i,lm14.NumObjects);
            patch(object, 'facecolor', [.8 .8 .8], 'facealpha', 0.3, 'edgecolor', 'none');
        end
        lm15(lm14.PixelIdxList{i}) = false;
    end
end

%% testing using small subset of cells

% lm14 = bwlabeln(lm13);
% figure();hold on;
% colorspec = [[0.5 0.5 1];[0.5 1 0.5];[1 0.5 0.5];[1 0.5 0];[1 0 0.5];...
%     [0 1 0.5];[0 0.5 1];[0.5 1 0];[0.5 0 1];[1 1 0.5];[1 0.5 1];[0.5 1 1];...
%     [1 1 0];[0 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1]];
% figure(1);
% for kk = 80:100
%     fprintf('Plotting object %4.0f of %4.0f...\n',kk,max(lm14(:)))
%     lm14(lm14~=kk) = 0;
%     blah = isosurface(lm14,0.3);
%     patch(blah, 'facecolor', colorspec(mod(kk,size(colorspec,1))+1,:), 'facealpha', 0.3, 'edgecolor', 'none') %, 'DisplayName', textbox);
%     lm14 = bwlabeln(lm13);
% end

%% skeleton plotting
% % skel, scurve, and scurvep should have been imported above.
% 
% syms t
% hold on;
% fprintf('Plotting skeleton...\n');
% for jj = 1:scurve.pieces
%     if mod(jj,10) == 0
%         fprintf('Working on segment %4.0f of %4.0f...\n',jj,scurve.pieces);
%     end
%     fxt=scurve.coefs(3*jj-2,:);
%     fyt=scurve.coefs(3*jj-1,:);
%     fzt=scurve.coefs(3*jj,:);
%     fx = fxt(1)*(t^3) + fxt(2)*(t^2) + fxt(3)*t + fxt(4);
%     fy = fyt(1)*(t^3) + fyt(2)*(t^2) + fyt(3)*t + fyt(4);
%     fz = fzt(1)*(t^3) + fzt(2)*(t^2) + fzt(3)*t + fzt(4);
%     
%     ezplot3(fx,fy,fz,[0,(scurve.breaks(jj+1) - scurve.breaks(jj))]);
% end
% 
% % fbreaks = [fx,fy,fz];


if strcmp(overwrite,'Y')
    fprintf('Workspace and base file (%s.mat) updated.\n',FileList{1});
    save(FileList{1},'allFiles','-v7.3','-append');
else
    fprintf('Workspace updated, but base file was not overwritten.\n');
end

%% Random stuff (don't actually run this)
return

%% aspect ratio vs. cell size
lm97(lm97(:,9)<10^-5,:)=[];
lm97(lm97(:,9)>1/1.2,:)=[];
figure('name','Aspect ratio vs cell size');

% use this format for future linreg
x_data = sqrt(lm97(:,11))*.151;
y_data = 1./lm97(:,9);
scatter(x_data,y_data); % view data
xx = [ones(length(x_data),1),x_data];
bb = xx\y_data;
hold on;
plot(xx(:,2),xx*bb); % plot line
rsq = 1 - sum((y_data - xx*bb).^2)/sum((y_data - mean(y_data)).^2) % gof value

set(gca,'FontSize',12);
title('Aspect Ratio vs. Cell Size','FontSize',20);
xlabel('Cell size (\mum)','FontSize',16);
ylabel('Cell apical aspect ratio','FontSize',16);
% aspect ratio vs. cell orientation

figure('name','Aspect ratio vs cell orientation');

scatter(abs(lm97(:,13)),1./lm97(:,9));
xx = [ones(length(lm97(:,13)),1),abs(lm97(:,13))];
yy = 1./lm97(:,9);
bb = xx\yy;
hold on;
plot(xx(:,2),xx*bb);
rsq = 1 - sum((yy - xx*bb).^2)/sum((yy - mean(yy)).^2)

set(gca,'FontSize',12);
title('Aspect Ratio vs. Cell Orientation','FontSize',20);
xlabel('Cell orientation (^{\circ})','FontSize',16);
ylabel('Cell apical aspect ratio','FontSize',16);
% cell orientation vs cell size
figure('name','Cell orientation vs cell size');

scatter(sqrt(lm97(:,11))*.151,abs(lm97(:,13)));
xx = [ones(length(lm97(:,11)),1),sqrt(lm97(:,11))*.151];
yy = abs(lm97(:,13));
bb = xx\yy;
hold on;
plot(xx(:,2),xx*bb);
rsq = 1 - sum((yy - xx*bb).^2)/sum((yy - mean(yy)).^2)

set(gca,'FontSize',12);
title('Cell orientation vs. cell size','FontSize',20);
xlabel('Cell size (\mum)','FontSize',16);
ylabel('Cell orientation (^{\circ})','FontSize',16);

%% more random stuff

% Send angle data to array
dataset1 = lm97(:,10).';

% Specify dataset2 (usually the mock DN one) and run the monte carlo
figure();plot(array_of_subsets/10);

% Compare boxplots of data distributions
a = [dataset1,dataset2];
group = [zeros(size(dataset1)),ones(size(dataset2))];
figure();boxplot(a,group);

% Rose angle histogram plot
figure();
rose(lm97(:,10)*pi/180);

% Compass plot of all orientations normalized to eccentricity
u = (cosd(lm97(:,10))).*lm97(:,6);
v = (sind(lm97(:,10))).*lm97(:,6);
figure();
compass(u,v);

% Histogram normalized to eccentricity?

edges=linspace(min(lm97(:,10)),max(lm97(:,10)),10);
[~, bin]=histc(lm97(:,10),edges);
count=accumarray(bin,lm97(:,6));

bar(count)



