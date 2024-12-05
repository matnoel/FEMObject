function P = calc_P_edges(M,D1,D2,ddls)
% function M = calc_P_edges(M,D1,D2,ddls)

% D1 et D2 : MODEL ou objet GEOMETRIQUE
% si D1 et D2 double : liste de noeuds
% ddls : liste des ddls a bloquer {'UX','U','R',...} pour la mecanique,
%        'T' pour la thermique ou autre ... selon le choix du MODEL
% ddls = 'all' par defaut

% Calcul de la matrice permettant de passer des ddl de M sur D1 aux ddl de
% M sur D2. Si v1 est un vecteur de ddl defini sur D1
% v2 = P*v1 est le vecteur de ddl correspondant defini sur D2
  
if M.nbddl==0
    error('ddls non definis dans le modele -> finaliser le modele')
end

if isa(D1,'GEOMOBJECT')
    [~,numnode1] = intersect(M,D1);
elseif isa(D1,'MODEL')
    [~,numnode1] = intersect(getnode(M),getnode(D1));
elseif isa(D1,'double')
    numnode1 = D1;
elseif isinteger(D1)
    numnode1 = double(D1) ;
else
    error('bad argument')
end
if isa(D2,'GEOMOBJECT')
    [~,numnode2] = intersect(M,D2);
elseif isa(D2,'MODEL')
    [~,numnode2] = intersect(getnode(M),getnode(D2));
elseif isa(D2,'double')
    numnode2 = D2;
elseif isinteger(D2)
    numnode2 = double(D2) ;
else
    error('bad argument')
end

if nargin<4
    ddls1 = getddlname(getddl(M.node,numnode1));
    ddls2 = getddlname(getddl(M.node,numnode2));
    ddls = intersect(ddls1,ddls2);
    if isempty(ddls)
        error('aucun ddl commun aux deux listes de noeuds')
    end
elseif isa(ddls,'char')
    ddls = {ddls};
end

xnode1 = double(getcoord(getnode(M),numnode1));
xnode2 = double(getcoord(getnode(M),numnode2));

for i=1:size(xnode1,2)
    [xnode1,I1] = sortrows(xnode1,i);
    numnode1 = numnode1(I1);
    [xnode2,I2] = sortrows(xnode2,i);
    numnode2 = numnode2(I2);
end

ddl1 = [];
ddl2 = [];
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
    numddl1 = findddl(M,ddl,numnode1);
    numddl2 = findddl(M,ddl,numnode2);
    ddl1 = union(ddl1,numddl1);
    ddl2 = union(ddl2,numddl2);
end

P = sparse(getnbddl(M),getnbddl(M));
P(ddl2,ddl1) = speye(length(ddl2));

