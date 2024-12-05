function M = sortnodecoord(M,varargin)
% function M = sortnodecoord(M,varargin)

[M.node,oldnumber] = sortnodecoord(M.node,varargin{:});

for j=1:M.nbgroupelem
   M.groupelem{j}=changenodenumber(M.groupelem{j},oldnumber,getnumber(M.node)) ;
end

