function M2 = transfercl(M1,M2,ddls)
% function M2 = transfercl(M1,M2,ddls)
% transfert des conditions aux limites d'un MODEL M1 vers un MODEL M2
% ddls : liste des ddls a transferer {'UX','U','R',...} pour la mecanique,
%        'T' pour la thermique ou autre ... selon le choix du MODEL
% ddls = 'all' par defaut

if M1.nbddl==0
    error('ddl non definis dans le premier modele -> finaliser le premier modele')
end
if M2.nbddl==0
    error('ddl non definis dans le second modele -> finaliser le second modele')
end

if nargin<3
    ddls = getddlname(getddl(M1.node,getnumber(M1.node)));
elseif isa(ddls,'char')
    ddls={ddls};
end

BC = getbc(M1.BCOND);

for i=1:length(BC)
    [numnode,ddl] = findnodeinnumddl(M1,BC{i}.ddlbloque);
    ddl = unique(ddl);
    if ismember(ddl,ddls)
        [~,~,repnode] = intersect(getnode(M1,numnode),getnode(M2));
        if ~isempty(repnode)
            numnode = getnumber(M2.node,repnode);
            value = BC{i}.value;
            M2 = addcl(M2,numnode,ddl,value);
        end
    end
end

