function PC=setindices(PC,ind,varargin)
% function PC=setindices(PC,ind)
% Set the multi-index set of POLYCHAOS PC
% 
% function PC=setindices(PC,ind,'sort')
% Set the multi-index set of POLYCHAOS PC and sort indices
% 
% function PC=setindices(PC,ind,'update')
% Set the multi-index set of POLYCHAOS PC and update PC
% 
% See also POLYCHAOS/getindices, POLYCHAOS/addindices, POLYCHAOS/removeindices

PC.indices=ind;
if ischarin('sort',varargin)
    PC.indices = sortrows(PC.indices,size(PC.indices,2):-1:1);
end
if ischarin('update',varargin)
    PC = update(PC);
end