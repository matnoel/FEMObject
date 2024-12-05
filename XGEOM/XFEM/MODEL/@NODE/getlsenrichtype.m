function l = getlsenrichtype(node,num);

if nargin==1
    l = node.lsenrichtype;
else
    l = zeros(size(num));
    [rep,j] = ismember(num,node.lsenrichnode);
    l(rep) = node.lsenrichtype(j(rep));
end
