function plot(T)

treeplot(T.connect);
nbnodes = size(T.connect,2);
set(gca,'xtick',[],'ytick',[])
[x,y] = treelayout(T.connect);
x = x';
y = y';

name1 = cellstr(num2str( (1:nbnodes)'));
name2= cell(nbnodes,1);
name3= cell(nbnodes,1);
for i=1:nbnodes
    name2{i} = num2str( T.Cvar{i} ) ;
    name2{i} = strcat('   [',name2{i},']');
%     if ~isempty(T.rank)
%     if ~isRleaf(T,i) && length(T.var{i})~=1
%         name3{i} = num2str( T.rank(i) ) ;
%         name3{i} = strcat('rank=[',name3{i},']   .');
%     else
%         name3{i}=' ';
%     end
%     end
end


text(x(:,1), y(:,1), name1, 'HorizontalAlignment','center','EdgeColor','red','BackgroundColor',[.7 .9 .7])
text(x(:,1), y(:,1), name2, 'HorizontalAlignment','left')
% if ~isempty(T.rank)
% text(x(:,1), y(:,1), name3, 'HorizontalAlignment','right')
% end
title({['Arbre ' inputname(1)]},'FontSize',20,'FontName','Times New Roman');