function M = createddlnode(M,varargin)
%function M = createddlnode(M)
% creation de la structure des ddl d'un modele
%
%function M = createddlnode(M,A)
% si A est un NODE ou MODEL, copie des ddl
% si A est DDL, on associe directement ces ddl aux noeuds

if nargin>1 && (isa(varargin{1},'NODE') || isa(varargin{1},'MODEL') )
    % ----- COPIE DES DDL SUR varargin{1}
    
    M = copy_groupddl(M,varargin{1});
    M = calcnumddl(M);
    %M = applyfunctiontofaces(M,@createddlnode,DDL('u'));
    M = applyfunctiontofaces(M,@copy_groupddl,getnode(M));
    M = applyfunctiontofaces(M,@calcnumddl);
    
elseif nargin>1 % --------------- CREATION DE NOUVEAUX DDL --------------
    M = actualise_ddl(M,varargin{:});
    M = create_groupddl(M);
    M = setbcond(M,BCOND(getnode(M)));
    M = calcnumddl(M);
    try
        M = applyfunctiontofaces(M,@createddlnode,getnode(M));
    end
else
    M = actualise_ddl(M,varargin{:});
    M = create_groupddl(M);
    M = setbcond(M,BCOND(getnode(M)));
    M = calcnumddl(M);
    M = applyfunctiontofaces(M,@copy_groupddl,getnode(M));
    M = applyfunctiontofaces(M,@calcnumddl);
end


