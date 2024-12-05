function [dim,isranddim] = getranddim(x)


isranddim = zeros(1,getnbgroups(x));
for i=1:getnbgroups(x)
isranddim(i) = numel(x.phi{i})>1;
end
dim = find(isranddim);