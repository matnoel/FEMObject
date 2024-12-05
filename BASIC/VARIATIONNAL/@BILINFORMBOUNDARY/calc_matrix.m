function b = calc_matrix(a,S,varargin)
% function b = calc_matrix(a,S,varargin)

if isa(S,'SEPMODEL')
    b = calc_matrix(S,a,varargin{:});
    return
end

if nargin==2
    varargin = cell(1,2);
    varargin(:) = {[]};
end

rep = find(~isemptyarg(varargin(1:getn(a))));
for i=1:length(rep)
    if isfloat(varargin{rep(i)}) && numel(varargin{rep(i)})>1
        varargin{rep(i)} = unfreevector(S,varargin{rep(i)});
    elseif isfloat(varargin{rep(i)})
        varargin{rep(i)} = varargin{rep(i)}*ones(getnbddl(S),1);
    end
end

if ~isempty(getk(a)) && ~isempty(getpk(a))
    try
        a = setk(a,unfreevector(S,getk(a))) ;
    end
end


if isboundary(a)
    Sini = S;
    S = create_boundary(Sini);
    if ~isempty(getboundaryobject(a))
        S = intersect(S,getboundaryobject(a));
    end
    S = createddlnode(S,getnode(Sini));
end


if isa(getselgroup(a),'char')
    sel = 1:getnbgroupelem(S);
else
    sel = getselgroup(a);
end

% L'utilisateur specifie un champs par element :
if ischarin('CHE',varargin)
    % On va plutot utiliser la fonction unassembled_matrixelem au lieu
    % de assembled_matrixelem !
    dim2change = getcharin('CHE',varargin);
    p = getp(a);
    if p(dim2change)==1
        error('On ne peut deriver un champ par element !');
    end
    
    % Donnee du type be(e,:,:,:)
    be = cell(1,getnbgroupelem(S));
    for p=sel
        be{p} = eval_elem(a,getgroupelem(S,p),getnode(S),varargin{:});
    end
    b = unassembled_matrixelem(S,be,'selgroup',sel);
    
    % Deuxieme etape : Permuter les dimensions 1 et dim2change
    % be(:,:,e,:)
    Perm = 1:b.ordre;
    Perm(1) = dim2change+1;
    Perm(dim2change+1) = 1;
    b = permute(b,Perm);
    % Derniere etape : Multiplication par le vecteur 1 :
    % be(1,:,e,:)
    ONE = one(S);
    b = SPTENSORtimesVECTOR(b,ONE,1);
    
    % L'utilisateur ne specifie pas de CHE
else % ANCIEN CODE (suite) :
    be = cell(1,getnbgroupelem(getcharin('parent',varargin,S)));
    selParent = [] ;
    for p=sel
        groupelem = getgroupelem(S,p) ;
        pParent = getparam(groupelem,'parentgroupelem') ;
        selParent = union(selParent,pParent) ;
        be{pParent} = eval_elem(a,groupelem,getnode(S),varargin{:});
    end
    b = assemble_matrixelem(getcharin('parent',varargin,S),be,...
        'selgroup',selParent);
    b = getfact(a)*b;
    if isboundary(a)
        P = calc_P(Sini,S);
        b = P'*b*P;
    end
    if isfree(a)
        b = freematrix(S,b);
    end
end

end
