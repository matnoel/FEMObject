function [dim,isdetdim] = getdetdim(x)

isdetdim = zeros(1,getnbdim(x));
for i=1:getnbdim(x)
isdetdim(i) = numel(x.phi{i})==1;
end
dim = find(isdetdim);