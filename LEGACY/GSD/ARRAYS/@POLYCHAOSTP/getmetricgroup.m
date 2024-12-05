function m = getmetricgroup(PC,k)
% function m = getmetricgroup(PC,k)

if nargin==1
    m = cellfun(@getmetric,PC.PCgroups,'UniformOutput',false);
else  
m = getmetric(PC.PCgroups{k});
end