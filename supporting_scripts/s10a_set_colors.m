%% use matlab set colors

a = get(gca,'ColorOrder');
r = a(:,1);
g = a(:,2);
b = a(:,3);

fprintf('1:Blue 2:Red 3:Yellow 4:Purple 5:Green 6:LightBlue 7:Maroon\n');

return

%% reroll rand values

 r = rand(1,30);
 g = rand(1,30);
 b = rand(1,30);