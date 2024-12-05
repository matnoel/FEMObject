function PC=removeindices(PC,ind,varargin)
% function PC=removeindices(PC,ind)
% Remove indices from the multi-index set of POLYCHAOS PC
% 
% function PC=removeindices(PC,ind,'sort')
% Remove indices from the multi-index set of POLYCHAOS PC and sort indices
% 
% function PC=removeindices(PC,ind,'update')
% Remove indices from the multi-index set of POLYCHAOS PC and update PC
% 
% See also POLYCHAOS/setindices, POLYCHAOS/getindices, POLYCHAOS/addindices

PC.indices = setdiff(PC.indices,ind,'rows','stable');
if ischarin('sort',varargin)
    PC.indices = sortrows(PC.indices,size(PC.indices,2):-1:1);
end
if ischarin('update',varargin)
    PC = update(PC);
end