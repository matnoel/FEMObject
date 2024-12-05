function M = addcl(M,D,ddls,value)
% function M = addcl(M,D,ddls,value)
% ajout de conditions aux limites au MODEL M
% D : MODEL ou objet GEOMETRIQUE
% si D=[] , on impose sur tout le bord de M
% si D double : liste de noeuds
% ddls : liste des ddls a bloquer {'UX','U','R',...} pour la mecanique,
%        'T' pour la thermique ou autre ... selon le choix du MODEL
% ddls = 'all' par defaut
% value : valeur imposee
% value = 0 par defaut

if M.nbddl==0
    error('ddls non definis dans le modele -> finaliser le modele')
end

if isempty(D)
    numnode = getnumber(getnode(create_boundary(M)));
elseif isa(D,'GEOMOBJECT')
    [~,numnode] = intersect(M,D);
elseif isa(D,'MODEL')
    [~,numnode] = intersect(getnode(M),getnode(D));
elseif isa(D,'double')
    numnode = D;
else
    error('bad argument')
end

if nargin<4
    value = 0;
end
if nargin<3
    ddls = getddlname(getddl(M.node,numnode));
elseif isa(ddls,'char')
    ddls = {ddls};
end

bloque = [];
nbddlcompo = 0;
ddlsto = DDL();
for l=1:length(ddls)
    switch ddls{l}
        case {'T'}
            ddl = DDL(DDLSCAL('T'));
        case {'U'}
            ddl = DDL(DDLVECT('U',M.syscoord,'TRANS'));
        case {'R'}
            ddl = DDL(DDLVECT('R',M.syscoord,'ROTA'));
        otherwise
            ddl = DDL(ddls{l});
    end
    ddlsto = addddl(ddlsto,ddl);
    ddlbloque = findddl(M,ddl,numnode);
    nbddlcompo = nbddlcompo+length(ddl);
    bloque = union(bloque,ddlbloque);
end
if isa(value,'double')
    if size(value,1)==1
        value  = repmat(value,nbddlcompo,1);
    end
    value  = repmat(value,length(numnode),1);
elseif isa(value,'function_handle') || isa(value,'inline')
    fun = fcnchk(value);
    x = double(getcoord(M.node,numnode));
    value = fun(x);
    value = value';
    value = value(:);
else
    error('wrong format for values')
end

M.BCOND = addbc(M.BCOND,bloque,value,ddlsto);
