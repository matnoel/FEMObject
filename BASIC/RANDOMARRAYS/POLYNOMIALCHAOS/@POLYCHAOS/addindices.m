function PC=addindices(PC,ind,varargin)
% function PC=addindices(PC,ind)
% Add indices to the multi-index set of POLYCHAOS PC with no repetitions
% 
% function PC=addindices(PC,ind,'sort')
% Add indices to the multi-index set of POLYCHAOS PC with no repetitions and sort indices
% 
% function PC=addindices(PC,ind,'update')
% Add indices to the multi-index set of POLYCHAOS PC with no repetitions and update PC
% 
% See also POLYCHAOS/setindices, POLYCHAOS/getindices, POLYCHAOS/removeindices

PC.indices = unique([PC.indices;ind],'rows','stable');
if ischarin('sort',varargin)
    PC.indices = sortrows(PC.indices,size(PC.indices,2):-1:1);
end
if ischarin('update',varargin)
    PC = update(PC);
end