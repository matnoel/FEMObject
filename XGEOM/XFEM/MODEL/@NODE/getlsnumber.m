function l = getlsnumber(node,num)

if nargin==1
    l = node.lsnumber;
else
    l = zeros(size(num));
    [rep,j] = ismember(num,node.lsenrichnode);
    l(rep) = node.lsnumber(j(rep));
end
