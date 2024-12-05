function plot(L)

T=getparam(L,'tree');
connect=getconnect(T);
var    =getvar(T);
maxorder =getparam(L,'maxorder');
tol =getparam(L,'tol');
treeplot(connect);
nbnodes = size(connect,2);
set(gca,'xtick',[],'ytick',[])
[x,y] = treelayout(connect);
x = x';
y = y';

name1 = cellstr(num2str( (1:nbnodes)'));
name2= cell(nbnodes,1);
name3= cell(nbnodes,1);
for i=1:nbnodes
    
    if var(i)~=0
        name2{i} = strcat('   [',num2str( var(i) ),']');
    else
        name2{i} = strcat('   m=',num2str( maxorder(i) ),' tol=',num2str( tol(i),'%1.1e\n' )  );
    end
    % if ~isempty(T.rank)
    % if ~isRleaf(T,i) && length(T.var{i})~=1
    %     name3{i} = num2str( T.rank(i) ) ;
    %     name3{i} = strcat('rank=[',name3{i},']   .');
    % else
    %     name3{i}=' ';
    % end
    % end
end


text(x(:,1), y(:,1), name1, 'HorizontalAlignment','center','EdgeColor','red','BackgroundColor',[.7 .9 .7])
text(x(:,1), y(:,1), name2, 'HorizontalAlignment','left')
% if ~isempty(T.rank)
%     text(x(:,1), y(:,1), name3, 'HorizontalAlignment','right')
% end
title({['Arbre ' inputname(1)]},'FontSize',20,'FontName','Times New Roman');