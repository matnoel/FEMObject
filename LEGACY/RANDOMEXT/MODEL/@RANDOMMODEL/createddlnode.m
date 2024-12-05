function M = createddlnode(M,varargin)

if nargin>1 & isa(varargin{1},'NODE')% ----- COPIE DES DDL SUR varargin{1}
 
[M.node,M.nbddl] = copy_groupddl(M.node,varargin{1});

else % --------------- CREATION DE NOUVEAUX DDL --------------
    
for p=1:M.nbgroupelem
M.groupelem{p}=actualise_ddl(M.groupelem{p},varargin{:});  
end

[M.node,M.nbddl] = create_groupddl(M.node,calcgroupddl(M));

M.BCOND = BCOND(M.node);

end

for i=1:M.nbgroupelem
M.groupelem{i}=calcnumddl(M.groupelem{i},M.node);
end
