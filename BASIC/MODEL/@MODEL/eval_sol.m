function v = eval_sol(S,u,P,ddl)
% function v = eval_sol(S,u,P,ddl)
% P : POINT ou on veut evaluer la solution
%   : si P est un MODEL , on prend les noeuds
% u : solution d'une resolution EF (double ou PCMATRIX)

if isa(P,'MODEL')
    P = POINT(P.node);
elseif ~isa(P,'POINT')
    error('P must be a POINT or a MODEL')
end

u = unfreevector(S,u);

if isa(u,'PCRADIALMATRIX')
    ud = double(getV(u));
else
    ud = double(u);
end

[Pinter,repP,numnode] = intersect(POINT(P),POINT(S.node));
numnode = getnumber(S.node,numnode);

repP2 = setdiff(1:numel(P),repP);

v = zeros(length(DDL(ddl)),size(ud,2),numel(P));
con = 0;
if numel(repP)>0
    rep = findddl(S,ddl,numnode);
    if ~isempty(rep)
        vtemp = full(ud(rep,:));
        vtemp = reshape(vtemp,[length(DDL(ddl)),length(numnode),size(vtemp,2)]);
        v(:,:,repP) = permute(vtemp,[1,3,2]);
        con = con + length(repP);
    end
end

if numel(repP2)>0
    for p=1:S.nbgroupelem
        [vtemp,repP] = eval_sol(S.groupelem{p},S.node,ud,P(repP2(:)),ddl);
        v(:,:,repP2(repP)) = vtemp ;
        con = con + length(repP);
    end
end

if isa(u,'PCMATRIX')
    v = permute(v,[1,3,2]);
    v = PCMATRIX(v,[size(v,1),size(v,2)],getPC(u));
elseif isa(u,'PCRADIALMATRIX')
    v = permute(v,[1,3,2]);
    v = MULTIMATRIX(v,[size(v,1),size(v,2)],[size(v,3),1]);
    v = PCRADIALMATRIX(v,size(v),getL(u));
elseif ~isa(u,'double')
    error('input is not a double')
end

if con<numel(P)
    warning('la solution n''est pas definie au point demande')
    v=[];
end
