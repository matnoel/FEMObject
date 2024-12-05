function M=calc_connec(M)


numelemglob=getnumelem(M);
inda=[];
indb=[];

for i=1:M.nbgroupelem
    numelem=getnumber(M.groupelem{i});
    co = getconnec(M.groupelem{i});
    [a,co] = ismember(co,getnumber(M.node)); % co : connectivite avec noeud numerotes de 1 a nbnode
    
    [numelem,numelemt]=intersect(numelemglob,numelem);
    e=repmat(numelemt(:),1,size(co,2));
    tempinda = co+M.nbnode*(e-1);
    inda=[inda;tempinda(:)];
    tempindb = e+M.nbelem*(co-1);
    indb=[indb;tempindb(:)];
end
[ia,ja]=ind2sub([M.nbnode,M.nbelem],inda);
[ib,jb]=ind2sub([M.nbelem,M.nbnode],indb);
A=sparse(ia,ja,1,M.nbnode,M.nbelem);
B=sparse(ib,jb,1,M.nbelem,M.nbnode);
%A(co+M.nbnode*(e-1))=1;
%B(e+M.nbelem*(co-1))=1;
C=B*A;
D=A*B;

M.connec.elem2node = A;
M.connec.node2elem = B;
M.connec.elem2elem = C;
M.connec.node2node = D;

% A  : taille nbnode*nbelem
% donne les elements connectes a un noeud
% ->ligne i contient des 1 aux colonnes
% correspondant aux elements connectes au noeud i

% B  : taille nbelem*nbnode
% donne les noeuds connectes a un element
% ->ligne i contient des 1 aux colonnes
% correspondant aux noeuds connectes a l'element i

% C  : taille nbelem*nbelem
% donne les elements connectes a un element et le nb de noeuds communs
% ->Cij est le nombre de noeuds communs entre les elements i et j
% Cii est naturellement la connectivite de l'element i

% D  : taille nbnode*nbnode
% donne les noeuds connectes a un noeud via des elements
% Cij est le nombre d'elements qui connectent le noeud i au noeud j
% Cii est donc le nombre d'elements connectes au noeud i
