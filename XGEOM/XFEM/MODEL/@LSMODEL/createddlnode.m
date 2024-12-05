function M = createddlnode(M,varargin)
%function M = createddlnode(M,varargin)

if nargin>1 && (isa(varargin{1},'NODE') || isa(varargin{1},'MODEL')) 
    % ----- COPIE DES DDL SUR varargin{1}
 
M = copy_groupddl(M,varargin{1});   


else % --------------- CREATION DE NOUVEAUX DDL --------------
    
M = actualise_ddl(M,varargin{:});
M = create_groupddl(M);
M = setbcond(M,BCOND(getnode(M)));
end


M = calcnumddl(M);



M = applyfunctiontofaces(M,@copy_groupddl,getnode(M));

