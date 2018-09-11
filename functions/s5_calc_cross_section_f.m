function allFiles = s5_calc_cross_section_f(allFiles)

if ~isstruct(allFiles)
    error('The input argument must be a struct of the H5 data.')
end

FileList = fieldnames(allFiles);

%% Calculate cross section

% [selected, ok] = listdlg('PromptString', 'Select file(s):', 'SelectionMode', 'multiple', 'ListString', varname);
%
% if (isempty(selected))
%     error('No files were selected.');
% else
%     fileNo = selected;
% end
%
% for ii = fileNo
%     cmd2005 = '?';
%     while(~strcmp(cmd2005, 'Y') && ~strcmp(cmd2005,'N'))
%         fprintf('\nProvide length per pixel for file %i?\n', ii);
%         cmd2005 = input('Y/N:  ', 's');
%         if isempty(cmd2005)
%             cmd2005 = '?';
%         end
%     end
%
%     switch cmd2005
%         case 'Y'
%             prompt = {'Length per pixel (in microns) for x dimension:','Length per pixel (in microns) for y dimension:', 'Length per pixel (in microns) for z dimension:'};
%             dlg_title = 'Input';
%             num_lines = 1;
%             default = {'1','1','1'};
%             answer = inputdlg(prompt, dlg_title, num_lines, default);
%
%             if ~isempty(answer)
%                 xscale = str2double(answer{1});
%                 yscale = str2double(answer{2});
%                 zscale = str2double(answer{3});
%             else
%                 xscale = 1;
%                 yscale = 1;
%                 zscale = 1;
%             end
%         case 'N'
%             fprintf('\nLength per pixel for file %i has defaulted to 1 by 1 by 1.\n', ii);
%             xscale = 1;
%             yscale = 1;
%             zscale = 1;
%     end
% end
%
% cmd2006 = zeros(size(1:fileNo, 2));
% cmd2007 = zeros(size(1:fileNo, 2));
%
% for ii = fileNo
%     % no_of_cs_samples = 0;
%     while (no_of_cs_samples == 0)
%         fprintf('\nHow many evenly spaced data points would you like to test for cross-sectional area for file %i?\n', ii);
%         no_of_cs_samples = input('Input:  ');
%         if (~isa(no_of_cs_samples, 'numeric'))
%             no_of_cs_samples = 0;
%             fprintf('\nNot a valid input.\n');
%         elseif (no_of_cs_samples <= 0)
%             fprintf('\nNot a valid input.\n');
%         end
%     end
% end

% for ii = fileNo

no_of_cs_samples = 100;
csdata = cell(no_of_cs_samples+1,4);
csdata{1,1} = 'Position';
csdata{1,2} = 'Image';
csdata{1,3} = 'Area';
csdata{1,4} = 'Perimeter';

%tic
fprintf('\nUsing data from: %s...\n', FileList{1});
allFiles.(FileList{1}).scurve = cscvn((allFiles.(FileList{1}).skel).');
% syms t;
allFiles.(FileList{1}).scurvep = fnder(allFiles.(FileList{1}).scurve);
% scurve = allFiles.(FileList{1}).scurve;

% PROMPT:  How many cross-sections do you want?
%                         cmd2006 = 0;
%                         while (cmd2006 == 0)
%                             fprintf('\nHow many evenly spaced data points would you like to test for cross-sectional area for file %i?\n', ii);
%                             cmd2006 = input('Input:  ');
%                             if (~isa(cmd2006, 'numeric'))
%                                 cmd2006 = 0;
%                                 fprintf('\nNot a valid input.\n');
%                             elseif (cmd2006 <= 0)
%                                 fprintf('\nNot a valid input.\n');
%                             end
%                         end
%
%                         cmd2007 = '?';
%
%                         while (~strcmp(cmd2007, 'Y') && ~strcmp(cmd2007,'N'))
%                             fprintf('\nWould you like the perimeter and area presented above the cross-section in the figure for file %i?\n', ii);
%                             cmd2007 = input('Y/N:  ', 's');
%                             if isempty(cmd2007)
%                                 cmd2007 = '?';
%                             end
%                         end

% fprintf('\nCalculating cross-section(s) for file %i...\n', ii);

% allFiles.(FileList{1}).NumberOfCrossSections = no_of_cs_samples;
% allFiles.(FileList{1}).NumberOfCrossSections = cmd2006;
% allFiles.(FileList{1}).csarea = [];
% allFiles.(FileList{1}).perimeter = [];
% range = allFiles.(FileList{1}).scurve.breaks(end)/(cmd2006 + 1):allFiles.(FileList{1}).scurve.breaks(end)/(cmd2006 + 1):cmd2006*allFiles.(FileList{1}).scurve.breaks(end)/(cmd2006 + 1);  % Number of data points + 1 = Number of segments
range = allFiles.(FileList{1}).scurve.breaks(end)/(no_of_cs_samples + 1):allFiles.(FileList{1}).scurve.breaks(end)/(no_of_cs_samples + 1):no_of_cs_samples*allFiles.(FileList{1}).scurve.breaks(end)/(no_of_cs_samples + 1);  % Number of data points + 1 = Number of segments

cv = fnval(allFiles.(FileList{1}).scurve, range);
cdv = fnval(fnder(allFiles.(FileList{1}).scurve), range);

univec = zeros(3, size(cdv,2));

for jj = 1:size(cdv, 2)
    univec(:,jj) = cdv(:,jj)/norm(cdv(:,jj));
end

% make "box" - boxloc(1) is xmin, (2) is xmax, (3) is ymin, (4) is
% ymax, (5) is zmin, (6) is zmax

pixwidth = 30; % SUBJECT TO CHANGE
h5bwdata = double(allFiles.(FileList{1}).h5bwfilled);

% NOTE:  VERY IMPORTANT TO UNDERSTANDING THIS CODE.
% x refers to columns and y refers to rows, due to translation from plots to
% matrices, the appropriate adjustments must be made.

% subplrow = -1;
%
% if size(cv,2) > 10
%     for bb = 10:-1:2
%         if mod(size(cv,2), bb) == 0
%             subplrow = bb;
%             subplcol = size(cv,2)/bb;
%             break;
%         end
%     end
% end
%
% if subplrow == -1
%     subplrow = 1;
%     subplcol = size(cv,2);
% end


for qq = 1:size(cv,2)
    fprintf('Working on object %3.0f of %3.0f (position %4.0f of %4.0f)...\n',qq,size(cv,2),range(qq),allFiles.(FileList{1}).scurve.breaks(end));
    csdata{qq+1,1} = range(qq);
    boxloc = [cv(1,qq) - pixwidth, cv(1,qq) + pixwidth, cv(2,qq) - pixwidth, cv(2,qq) + pixwidth, cv(3,qq) - pixwidth, cv(3,qq) + pixwidth];
    
    boxloc(boxloc<1)=1;
    if boxloc(2) > size(allFiles.(FileList{1}).h5bwfilled,2)
        boxloc(2) = size(allFiles.(FileList{1}).h5bwfilled,2);
    end
    if boxloc(4) > size(allFiles.(FileList{1}).h5bwfilled,1)
        boxloc(4) = size(allFiles.(FileList{1}).h5bwfilled,1);
    end
    if boxloc(6) > size(allFiles.(FileList{1}).h5bwfilled,3)
        boxloc(6) = size(allFiles.(FileList{1}).h5bwfilled,3);
    end
    
    % try
    h5bwdatabox = h5bwdata(floor(boxloc(3)):ceil(boxloc(4)), floor(boxloc(1)): ceil(boxloc(2)), floor(boxloc(5)):ceil(boxloc(6)));
    
    [xind, yind, zind] = ind2sub(size(h5bwdatabox), find(h5bwdatabox == 1));
    xyzind = [yind + floor(boxloc(1)) - 1, xind + floor(boxloc(3)) - 1, zind + floor(boxloc(5)) - 1];
    
    vdist2plane = zeros(size(xyzind, 1), 1);
    
    for kk = 1:size(xyzind, 1)
        vdist2plane(kk) = dot(univec(:,qq), xyzind(kk,:)' - cv(:,qq));
    end
    
    v = abs(vdist2plane);
    sortind = find(v(:) <= 1);
    
    if isempty(sortind)
        warning('Invalid coords for this cross-section. Skipping...');
        continue
    end
    
    vlist = zeros(size(sortind,1),3);
    
    for kk = 1:size(sortind, 1)
        vlist(kk,:) = xyzind(sortind(kk), :);
    end
    
    vlist = unique(vlist, 'rows');
    
    pproj = zeros(size(vlist, 1), 3);
    for kk = 1:size(vlist, 1)
        pproj(kk, :) = vlist(kk, :) - dot(vlist(kk, :) - cv(:, qq).', univec(:,qq)) * univec(:,qq).';
    end
    
    while (max(pproj(:,3)) - min(pproj(:,3)) > 1)
        w = [0 0 0];
        while all(w == [0 0 0])
            rng1 = randi(size(pproj,1));
            w = cross(pproj(randi(size(pproj,1)), :) - pproj(rng1, :), pproj(randi(size(pproj,1)), :) - pproj(rng1, :));
        end
        
        w = w/norm(w);
        R = [null(w), w.'];
        
        if det(R) < 0
            R(:, 1:2) = R(:, 2:-1:1);
        end
        pproj = pproj*R;
    end
    
    xuunit = [1 0 0]*R^(-1);
    yvunit = [0 1 0]*R^(-1);
    uunit = ((xuunit(1)*xscale)^2+(xuunit(2)*yscale)^2+(xuunit(3)*zscale)^2)^(1/2);
    vunit = ((yvunit(1)*xscale)^2+(yvunit(2)*yscale)^2+(yvunit(3)*zscale)^2)^(1/2);
    
    
    csimage = zeros(ceil(max(pproj(:,1)) - min(pproj(:,1))) + 10, ceil(max(pproj(:,2)) - min(pproj(:,2))) + 10);
    
    for p = 1:size(pproj,1)
        csimage(round(pproj(p,1) - min(pproj(:,1))) + 6, round(pproj(p,2) - min(pproj(:,2))) + 6) = 1;
    end
    
    [blobImage, numberOfBlombs] = bwlabel(csimage);
    blobMeasurements = regionprops(blobImage, 'area');
    allAreas = [blobMeasurements.Area];
    
    [sortedAreas, sortIndexes] = sort(allAreas, 'descend');
    biggestBlob = ismember(blobImage, sortIndexes(1));
    csobj = biggestBlob > 0;
    
    csobj = imfill(csobj, 'holes');
    csobj_rp = regionprops(csobj,'Area','Perimeter');
    % allFiles.(FileList{1}).csarea(qq) = area.Area * uunit * vunit;
    
    csdata{qq+1,2} = csobj;
    csdata{qq+1,3} = csobj_rp.Area * uunit * vunit;
    csdata{qq+1,4} = csobj_rp.Perimeter * uunit * vunit;
    %
    %         ok1 = 0 ;
    %         u = abs(uunit);
    %         v = abs(vunit);
    %         z = sqrt(uunit^2 + vunit^2);
    %
    %         while ok1 ~= 1
    %
    %             objotl = bwperim(csobj);
    %
    %             allFiles.(FileList{1}).perimeter(qq) = 0;
    %             method = 0;
    %             [row, col] = find(objotl == 1);
    %             current = [row(1), col(1)];
    %             num = 0;
    %             if objotl(current(1) + 1, current(2)) == 1 % down
    %                 current = [current(1) + 1, current(2)];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                 method = 1;
    %                 num = num + 1;
    %             elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                 current = [current(1) - 1, current(2)];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                 method = 2;
    %                 num = num + 1;
    %             elseif objotl(current(1), current(2) + 1) == 1 % right
    %                 current = [current(1), current(2) + 1];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                 method = 3;
    %                 num = num + 1;
    %             elseif objotl(current(1), current(2) - 1) == 1 % left
    %                 current = [current(1), current(2) - 1];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                 method = 4;
    %                 num = num + 1;
    %             elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                 current = [current(1) + 1, current(2) + 1];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                 method = 5;
    %                 num = num + 1;
    %             elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                 current = [current(1) + 1, current(2) - 1];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                 method = 6;
    %                 num = num + 1;
    %             elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                 current = [current(1) - 1, current(2) + 1];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                 method = 7;
    %                 num = num + 1;
    %             elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                 current = [current(1) - 1, current(2) - 1];
    %                 allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                 method = 8;
    %                 num = num + 1;
    %             else
    %                 csobj = bwmorph(csobj, 'spur');
    %                 continue;
    %             end
    %             ok2 = 1;
    %             while (num ~= (size(row, 1)))
    %                 switch method
    %                     case 1
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 2
    %                         if objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 3
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 4
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 5
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 6
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 7
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) + 1) == 1 % right/down
    %                             current = [current(1) + 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 5;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                     case 8
    %                         if objotl(current(1) + 1, current(2)) == 1 % down
    %                             current = [current(1) + 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 1;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2)) == 1 % up
    %                             current = [current(1) - 1, current(2)];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + v;
    %                             method = 2;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) + 1) == 1 % right
    %                             current = [current(1), current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 3;
    %                             num = num + 1;
    %                         elseif objotl(current(1), current(2) - 1) == 1 % left
    %                             current = [current(1), current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + u;
    %                             method = 4;
    %                             num = num + 1;
    %                         elseif objotl(current(1) + 1, current(2) - 1) == 1 % left/down
    %                             current = [current(1) + 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 6;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) + 1) == 1 % right/up
    %                             current = [current(1) - 1, current(2) + 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 7;
    %                             num = num + 1;
    %                         elseif objotl(current(1) - 1, current(2) - 1) == 1 % left/up
    %                             current = [current(1) - 1, current(2) - 1];
    %                             allFiles.(FileList{1}).perimeter(qq) = allFiles.(FileList{1}).perimeter(qq) + z;
    %                             method = 8;
    %                             num = num + 1;
    %                         else
    %                             csobj = bwmorph(csobj, 'spur');
    %                             ok2 = 0;
    %                             break;
    %                         end
    %                 end
    %             end
    %             if ok2 == 0
    %                 continue;
    %             end
    %             ok1 = 1;
    %         end
    %
    %
    % subplot(subplrow,subplcol, qq);
    % imshow(csobj);
    % title({['Area:  ', num2str(allFiles.(FileList{1}).csarea(qq))]; ['Perimeter:  ', num2str(allFiles.(FileList{1}).perimeter(qq))]});
    
    % figure();imshow(csobj)
    
    % textbox = {['Area:  ', num2str(allFiles.(FileList{1}).csarea(qq))]; ['Perimeter:  ', num2str(allFiles.(FileList{1}).perimeter(qq))]};
    
    
    
    % scatter3(vlist(:,1), vlist(:,2), vlist(:,3)); text(cv(1,qq), cv(2,qq), cv(3,qq) + 15, textbox);
    
    
    % catch
    %    warning('Cross section area out of bounds. Skipping...');
    % end
end


allFiles.(FileList{1}).csdata = csdata;


return

%% Debug only

%% save


if strcmp(overwrite,'Y')
    fprintf('Workspace and base file (%s.mat) updated.\n',FileList{1});
    save(FileList{1},'allFiles','-v7.3','-append');
else
    fprintf('Workspace updated, but base file was not overwritten.\n');
end



end
