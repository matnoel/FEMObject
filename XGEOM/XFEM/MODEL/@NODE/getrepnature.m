function l = getrepnature(node,num)

if nargin==1
    l = node.repnature;
else
    l = zeros(size(num));
    [rep,j] = ismember(num,node.lsenrichnode);
    l(rep) = node.repnature(j(rep));
end

