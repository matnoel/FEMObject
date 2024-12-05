function m = getmassegroup(PC,k)
%  function m = getmassegroup(PC,k)
if nargin==1
    m = cellfun(@getmasse,PC.PCgroups,'UniformOutput',false);
else  
m = getmasse(PC.PCgroups{k});
end