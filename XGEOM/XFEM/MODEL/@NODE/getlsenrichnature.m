function u = getlsenrichnature(node,num)
if nargin==1
u = node.lsenrichnature;
else
rep = getrepnature(node,num);
temp = find(rep);
u = cell(size(num,1),size(num,2));
u(temp)=node.lsenrichnature(rep(temp));
if numel(num)==1
u = u{1};    
end
end

