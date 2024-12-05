function g = getnbddlpergroup(node,num)

if nargin == 1
    num = 1:length(node.groupddl);
end

g = zeros(1,length(num));
for k=1:length(num)
    g(k) = node.groupddl{num(k)}.nbddl;
end
